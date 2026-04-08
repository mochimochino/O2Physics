// ============================================================
//  Step3_Unfolding.C (Bayesian Unfolding & Correlation Only)
// ============================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>
#include <vector>

// RooUnfold headers
#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ROOUNFOLD_ROOT / include)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include <TSystem.h>

// ----------------------------------------------------------
// Calculate Average Global Correlation
// ----------------------------------------------------------
double GetAverageGlobalCorrelation(TMatrixD covMat)
{
  int nBins = covMat.GetNrows();
  TMatrixD covInv = covMat;
  Double_t det = 0;
  covInv.Invert(&det);

  double sumRho = 0.0;
  int validBins = 0;

  for (int i = 0; i < nBins; ++i) {
    double v_ii = covMat(i, i);
    double vinv_ii = covInv(i, i);

    if (v_ii > 0 && vinv_ii > 0) {
      double term = 1.0 - 1.0 / (v_ii * vinv_ii);
      if (term < 0)
        term = 0.0;
      double rho_i = std::sqrt(term);
      sumRho += rho_i;
      validBins++;
    }
  }
  if (validBins > 0)
    return sumRho / validBins;
  return 1.0;
}

void Step3_Unfolding()
{
  // Settings
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0408/Resolution/";
  const TString inFile = inDir + "Step2_Response_for_Unfolding.root";
  const TString outDir = inDir;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  gStyle->SetPalette(kLightTemperature);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open " << inFile << std::endl;
    return;
  }

  TH1D* hGen = (TH1D*)fIn->Get("hMCGen_Rebinned");
  TH1D* hReco = (TH1D*)fIn->Get("hMCReco_Rebinned");
  TH2D* hMat = (TH2D*)fIn->Get("hResponseMatrix");
  TH1D* hData = (TH1D*)fIn->Get("hDataReco_Rebinned");

  int nBins = hGen->GetNbinsX();

  std::cout << "[Info] Constructing RooUnfoldResponse..." << std::endl;
  RooUnfoldResponse response(hReco, hGen, hMat);

  // ==========================================================
  // 1. Global Correlation Scan (Bayesian Unfolding Optimization)
  // ==========================================================
  std::cout << "\n[Info] Running Global Correlation Scan for Bayes..." << std::endl;

  TGraph* grScanBayes = new TGraph();
  grScanBayes->SetTitle("Bayesian Unfolding Optimization;Number of Iterations;Average Global Correlation #LT#rho#GT");
  grScanBayes->SetMarkerStyle(20);
  grScanBayes->SetMarkerColor(kRed);
  grScanBayes->SetLineColor(kRed);
  grScanBayes->SetLineWidth(2);

  int bestBayesIter = 1;
  double minRhoBayes = 1.0;

  for (int iter = 1; iter <= 15; ++iter) {
    RooUnfoldBayes unfoldBayes(&response, hData, iter);
    TMatrixD covMat = unfoldBayes.Ereco(RooUnfold::kCovariance);
    double avgRho = GetAverageGlobalCorrelation(covMat);

    grScanBayes->SetPoint(grScanBayes->GetN(), iter, avgRho);

    if (avgRho < minRhoBayes) {
      minRhoBayes = avgRho;
      bestBayesIter = iter;
    }
  }

  std::cout << " -> Best Bayes Iteration : " << bestBayesIter << " (Avg Rho = " << minRhoBayes << ")" << std::endl;

  TCanvas* cScan = new TCanvas("cScan", "Global Correlation Scan", 600, 500);
  cScan->SetGrid();
  grScanBayes->Draw("ALP");
  cScan->SaveAs(outDir + "Step3_GlobalCorrelationScan_Bayes.png");

  // ==========================================================
  // 2. Unfolding and Closure Test with Optimal Parameters
  // ==========================================================
  std::cout << "\n[Info] Running Unfolding with Optimal Iteration (" << bestBayesIter << ")..." << std::endl;

  RooUnfoldBayes unfoldBayesOpt(&response, hData, bestBayesIter);
  TH1D* hUnfoldedBayes = (TH1D*)unfoldBayesOpt.Hreco();
  hUnfoldedBayes->SetName("hUnfoldedBayes");
  hUnfoldedBayes->SetLineColor(kRed);
  hUnfoldedBayes->SetMarkerColor(kRed);
  hUnfoldedBayes->SetMarkerStyle(20);

  hGen->SetLineColor(kBlack);
  hGen->SetLineWidth(2);
  TCanvas* c1 = new TCanvas("c1", "Closure Test", 800, 800);
  c1->Divide(1, 2);

  TPad* pad1 = (TPad*)c1->GetPad(1);
  pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  hGen->SetTitle(Form("Closure Test: Unfolded vs True MC (Opt. Iteration = %d)", bestBayesIter));
  hGen->GetXaxis()->SetLabelSize(0);
  hGen->GetXaxis()->SetTitleSize(0);
  hGen->Draw("HIST");
  hUnfoldedBayes->Draw("PE SAME");

  TLegend* leg = new TLegend(0.55, 0.7, 0.85, 0.85);
  leg->AddEntry(hGen, "True (MC Gen)", "l");
  leg->AddEntry(hUnfoldedBayes, Form("Bayes (iter=%d)", bestBayesIter), "pe");
  leg->Draw();

  c1->cd();
  TPad* pad2 = (TPad*)c1->GetPad(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1D* hRatioBayes = (TH1D*)hUnfoldedBayes->Clone("hRatioBayes");
  hRatioBayes->Divide(hGen);
  hRatioBayes->SetTitle("");
  hRatioBayes->GetYaxis()->SetTitle("Unfolded / True");
  hRatioBayes->GetYaxis()->SetRangeUser(0.5, 1.5);
  hRatioBayes->GetYaxis()->SetNdivisions(505);
  hRatioBayes->GetYaxis()->SetLabelSize(0.1);
  hRatioBayes->GetYaxis()->SetTitleSize(0.12);
  hRatioBayes->GetYaxis()->SetTitleOffset(0.4);
  hRatioBayes->GetXaxis()->SetTitle("p_{T}^{2} (GeV^{2}/c^{2})");
  hRatioBayes->GetXaxis()->SetLabelSize(0.1);
  hRatioBayes->GetXaxis()->SetTitleSize(0.12);
  hRatioBayes->Draw("PE");

  TLine* line = new TLine(hGen->GetXaxis()->GetXmin(), 1.0, hGen->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw("SAME");
  c1->SaveAs(outDir + "Step3_ClosureTest_BayesOnly.png");

  // ==========================================================
  // 3. 2D Correlation Matrix
  // ==========================================================
  std::cout << "[Info] Drawing 2D Correlation Matrix for optimal iteration..." << std::endl;

  TMatrixD covMatOpt = unfoldBayesOpt.Ereco(RooUnfold::kCovariance);
  TH2D* hCorr = new TH2D("hCorr", Form("Correlation Matrix (Bayes iter=%d);Bin i;Bin j", bestBayesIter),
                         nBins, 0.5, nBins + 0.5, nBins, 0.5, nBins + 0.5);

  for (int i = 0; i < nBins; ++i) {
    for (int j = 0; j < nBins; ++j) {
      double var_i = covMatOpt(i, i);
      double var_j = covMatOpt(j, j);
      if (var_i > 0 && var_j > 0) {
        // ρ_ij = V_ij / sqrt(V_ii * V_jj)
        double corr = covMatOpt(i, j) / std::sqrt(var_i * var_j);
        hCorr->SetBinContent(i + 1, j + 1, corr);
      } else {
        hCorr->SetBinContent(i + 1, j + 1, 0.0);
      }
    }
  }

  TCanvas* cCorr = new TCanvas("cCorr", "Correlation Matrix", 800, 600);
  cCorr->SetRightMargin(0.15);
  hCorr->GetZaxis()->SetRangeUser(-1.0, 1.0);
  gStyle->SetPaintTextFormat(".2f");
  hCorr->Draw("COLZ TEXT");
  cCorr->SaveAs(outDir + "Step3_CorrelationMatrix_BayesOpt.png");

  std::cout << "[Done] Unfolding optimization and validation complete." << std::endl;
}
