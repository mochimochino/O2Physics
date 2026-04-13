// ============================================================
//  Step3_1_Unfolding.C (Bayesian Unfolding 4-Pad Comparison)
// ============================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include <iostream>
#include <vector>

// RooUnfold headers
#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ROOUNFOLD_ROOT / include)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

void Step3_1_Unfolding()
{
  // --- Settings ---
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/IncoherentJpsi/50/";
  const TString inFile = inDir + "Step2_Response_for_Unfolding.root";
  const TString outDir = inDir;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open " << inFile << std::endl;
    return;
  }

  TH1D* hGen = (TH1D*)fIn->Get("hMCGen_Rebinned");     // Truth (MC Gen)
  TH1D* hReco = (TH1D*)fIn->Get("hMCReco_Rebinned");   // MC Reco
  TH2D* hMat = (TH2D*)fIn->Get("hResponseMatrix");     // Response Matrix
  TH1D* hData = (TH1D*)fIn->Get("hDataReco_Rebinned"); // Measured (Data)

  std::cout << "[Info] Constructing RooUnfoldResponse..." << std::endl;
  RooUnfoldResponse response(hReco, hGen, hMat);

  std::vector<int> iters = {1, 3, 5, 7}; // iteration

  // 2x2
  TCanvas* c1 = new TCanvas("c1", "Bayes Iteration Comparison", 1200, 1000);
  c1->Divide(2, 2);

  for (size_t i = 0; i < iters.size(); ++i) {
    int iter = iters[i];

    RooUnfoldBayes unfoldBayes(&response, hData, iter);
    TH1D* hUnfolded = (TH1D*)unfoldBayes.Hreco()->Clone(Form("hUnfolded_iter%d", iter));
    hUnfolded->SetLineColor(kRed + 1);
    hUnfolded->SetMarkerColor(kRed + 1);
    hUnfolded->SetMarkerStyle(20);

    c1->cd(i + 1);

    // ==========================================================
    // Top Pad: Spectra
    // ==========================================================
    TPad* padTop = new TPad(Form("padTop_%zu", i), "padTop", 0.0, 0.3, 1.0, 1.0);
    padTop->SetBottomMargin(0.02);
    padTop->SetLeftMargin(0.12);
    padTop->SetRightMargin(0.05);
    padTop->SetLogy();
    padTop->Draw();
    padTop->cd();

    TH1D* hGenDraw = (TH1D*)hGen->Clone(Form("hGenDraw_%zu", i));
    hGenDraw->SetTitle(Form("Bayes Unfolding (Iter = %d)", iter));
    hGenDraw->SetLineColor(kBlack);
    hGenDraw->SetLineWidth(2);
    hGenDraw->GetXaxis()->SetLabelSize(0);
    hGenDraw->GetXaxis()->SetTitleSize(0);
    hGenDraw->GetYaxis()->SetTitle("Counts");
    hGenDraw->GetYaxis()->SetTitleSize(0.05);
    hGenDraw->GetYaxis()->SetLabelSize(0.045);

    double maxY = hGenDraw->GetMaximum();
    double minY = hGenDraw->GetMinimum();
    if (hData->GetMaximum() > maxY)
      maxY = hData->GetMaximum();
    if (hUnfolded->GetMaximum() > maxY)
      maxY = hUnfolded->GetMaximum();
    hGenDraw->SetMaximum(maxY * 10.0);
    hGenDraw->SetMinimum(minY * 0.1);

    hGenDraw->Draw("HIST");

    TH1D* hDataDraw = (TH1D*)hData->Clone(Form("hDataDraw_%zu", i));
    hDataDraw->SetLineColor(kBlue + 1);
    hDataDraw->SetMarkerColor(kBlue + 1);
    hDataDraw->SetMarkerStyle(24);
    // hDataDraw->SetLineWidth(2);
    hDataDraw->Draw("PE SAME");

    hUnfolded->Draw("PE SAME");

    TLegend* leg = new TLegend(0.65, 0.65, 0.90, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hGenDraw, "Truth", "l");
    leg->AddEntry(hDataDraw, "Measured", "l");
    leg->AddEntry(hUnfolded, Form("Unfolded (iter=%d)", iter), "pe");
    leg->Draw();

    // ==========================================================
    // Bottom Pad: Ratio to Truth
    // ==========================================================
    c1->cd(i + 1);
    TPad* padBot = new TPad(Form("padBot_%zu", i), "padBot", 0.0, 0.0, 1.0, 0.3);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.35);
    padBot->SetLeftMargin(0.12);
    padBot->SetRightMargin(0.05);
    padBot->Draw();
    padBot->cd();

    TH1D* hRatio = (TH1D*)hUnfolded->Clone(Form("hRatio_%zu", i));
    hRatio->Divide(hGen); // Unfolded / Truth
    hRatio->SetTitle("");

    hRatio->GetYaxis()->SetTitle("Unfolded/Truth");
    hRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetTitleOffset(0.5);

    hRatio->GetXaxis()->SetTitle("p_{T}^{2} (GeV^{2}/c^{2})");
    hRatio->GetXaxis()->SetLabelSize(0.12);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->Draw("PE");

    // y = 1.0
    TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->Draw("SAME");
  }

  c1->SaveAs(outDir + "Step3_1_4Pads_Comparison.png");

  std::cout << "[Done] Step3-1: 4-Pad Iteration comparison plot created successfully." << std::endl;
}
