#include "Math/MinimizerOptions.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>

// ==========================================
// Function definitions
// ==========================================
Double_t pT_func(Double_t* x, Double_t* par)
{
  double pt = x[0];
  double N = par[0];
  double b = par[1]; // b_pd
  double n = par[2]; // n_pd

  if (n <= 0)
    return 0; // Prevent division by zero or negative powers

  return N * pt * TMath::Power(1.0 + (b / n) * pt * pt, -n);
}

// ==========================================
// Main macro
// ==========================================
void pT_fit()
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFontSize(0.04);

  TString filePath = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/pT/AnalysisResults.root";
  TString dirName = "my-upc-mass-jpsi/";
  int rebinFactor = 2;
  double drawMin = 0.0;
  double drawMax = 2.5;
  double fitMin = 0.01;
  double fitMax = 2.5;

  // 1. ファイルとヒストグラムの読み込み (フェイルセーフのみ残す)
  TFile* f = TFile::Open(filePath, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open file: " << filePath << std::endl;
    return;
  }

  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dirName + "PtMuonUnlike"));
  TH1* hLike = dynamic_cast<TH1*>(f->Get(dirName + "PtMuonLike"));
  if (!hUnlike || !hLike) {
    std::cerr << "ERROR: Cannot get histograms from: " << dirName << std::endl;
    return;
  }

  // Rebin
  hUnlike->Rebin(rebinFactor);
  hLike->Rebin(rebinFactor);

  // Signal Extraction (Unlike - Like)
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->Add(hLike, -1.0);

  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);
  hSignal->SetTitle(Form("Signal p_{T} (Unlike - Like); p_{T} [GeV/c]; Counts / (%.3f GeV/c)", binWidth));
  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);

  // pT Fit
  double maxVal = hSignal->GetMaximum();
  TF1* fitFunc = new TF1("fitFunc", pT_func, fitMin, fitMax, 3);
  fitFunc->SetParameters(maxVal * 10.0, 40.0, 2.0);
  fitFunc->SetParLimits(0, 0.0, maxVal * 1000.0);
  fitFunc->SetParLimits(1, 1.0, 10000.0);
  fitFunc->SetParLimits(2, 0.1, 20.0);

  std::cout << "\n=== pT Fit ===" << std::endl;
  hSignal->Fit("fitFunc", "RME0");

  double chi2 = fitFunc->GetChisquare();
  int ndf = fitFunc->GetNDF();
  double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
  double paramB = fitFunc->GetParameter(1);
  double paramB_err = fitFunc->GetParError(1);
  double paramn = fitFunc->GetParameter(2);
  double paramn_err = fitFunc->GetParError(2);

  // 5. 描画設定
  TCanvas* c = new TCanvas("c", "pT Fit Canvas", 800, 800);

  TPad* p1 = new TPad("p1", "p1", 0, 0.3, 1, 1.0);
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(0.12);
  p1->SetRightMargin(0.05);
  p1->SetTopMargin(0.05);
  p1->Draw();

  TPad* p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.3);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.30);
  p2->SetLeftMargin(0.12);
  p2->SetRightMargin(0.05);
  p2->Draw();

  // --- Upper Pad (Main Plot) ---
  p1->cd();
  gPad->SetGrid();

  TH1* hS = (TH1*)hSignal->Clone("hS");
  hS->GetXaxis()->SetLabelSize(0);
  hS->GetXaxis()->SetTitleSize(0);
  hS->Draw("E");

  fitFunc->SetLineColor(kBlue);
  fitFunc->SetLineWidth(2);
  fitFunc->Draw("SAME");

  TLegend* leg = new TLegend(0.55, 0.75, 0.90, 0.90);
  leg->SetBorderSize(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(hS, "Data (Unlike - Like)", "lep");
  leg->AddEntry(fitFunc, "Fit: p_{T}(1 + #frac{b}{n}p_{T}^{2})^{-n}", "l");
  leg->Draw();

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextFont(42);
  latex->DrawLatex(0.55, 0.65, Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2, ndf, chi2ndf));
  latex->DrawLatex(0.55, 0.58, Form("b_{pd} = %.2f #pm %.2f (GeV/c)^{-2}", paramB, paramB_err));
  latex->DrawLatex(0.55, 0.51, Form("n_{pd} = %.2f #pm %.2f", paramn, paramn_err));

  // --- Lower Pad (Residuals) ---
  p2->cd();
  gPad->SetGridy();
  TH1D* hR = (TH1D*)hS->Clone("hR");
  hR->SetTitle("");
  for (int i = 1; i <= hR->GetNbinsX(); i++) {
    double xc = hR->GetBinCenter(i);
    if (xc >= drawMin && xc <= drawMax) {
      double d = hS->GetBinContent(i);
      double fv = fitFunc->Eval(xc);
      if (fv > 0) {
        hR->SetBinContent(i, (d - fv) / fv);
        hR->SetBinError(i, hS->GetBinError(i) / fv);
      } else {
        hR->SetBinContent(i, 0);
        hR->SetBinError(i, 0);
      }
    } else {
      hR->SetBinContent(i, 0);
      hR->SetBinError(i, 0);
    }
  }

  hR->SetMarkerStyle(20);
  hR->GetYaxis()->SetTitle("(Data-Fit)/Fit");
  hR->GetYaxis()->SetRangeUser(-1.5, 1.5);
  hR->GetYaxis()->SetNdivisions(505);
  hR->GetYaxis()->SetLabelSize(0.10);
  hR->GetYaxis()->SetTitleSize(0.12);
  hR->GetYaxis()->SetTitleOffset(0.4);

  hR->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hR->GetXaxis()->SetLabelSize(0.10);
  hR->GetXaxis()->SetTitleSize(0.12);
  hR->GetXaxis()->SetTitleOffset(0.9);
  hR->Draw("EP");

  TF1* l0 = new TF1("l0", "0", drawMin, drawMax);
  l0->SetLineColor(kRed);
  l0->SetLineStyle(2);
  l0->Draw("SAME");

  c->SaveAs("pT_fit_result.png");
  f->Close();
}
