#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include <iostream>
#include <vector>

void CheckPurityStability(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0414/Coherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0414/test/Logic1_PurityStability.root",
  const std::vector<double>& binEdges = {0.0, 0.010, 0.065}
) {
  gStyle->SetOptStat(0);
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) return;
  TH2D* h2 = (TH2D*)fIn->Get("hResponseMatrixPairPt2");
  if (!h2) return;

  int nBins = binEdges.size() - 1;
  const double* bins = &binEdges[0];

  TH1D* hPurity = new TH1D("hPurity", "Purity;p_{T}^{2} (GeV^{2}/c^{2});Fraction", nBins, bins);
  TH1D* hStability = new TH1D("hStability", "Stability;p_{T}^{2} (GeV^{2}/c^{2});Fraction", nBins, bins);

  for (int i = 1; i <= nBins; ++i) {
    double xLo = binEdges[i-1];
    double xHi = binEdges[i];
    int binLo = h2->GetXaxis()->FindBin(xLo + 1e-6);
    int binHi = h2->GetXaxis()->FindBin(xHi - 1e-6);

    double diag = h2->Integral(binLo, binHi, binLo, binHi);
    double allReco = h2->Integral(1, h2->GetNbinsX(), binLo, binHi);
    double allGen = h2->Integral(binLo, binHi, 1, h2->GetNbinsY());

    if (allReco > 0) hPurity->SetBinContent(i, diag / allReco);
    if (allGen > 0) hStability->SetBinContent(i, diag / allGen);
  }

  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  hPurity->Write();
  hStability->Write();
  fOut->Close();
  std::cout << "[Info] Saved ROOT: " << outputFile << std::endl;

  TString outBase = outputFile;
  outBase.ReplaceAll(".root", "");

  auto drawAndSave = [&](TH1D* h, const TString& title, const TString& suffix, int color) {
    TCanvas* c = new TCanvas("c_" + suffix, title, 800, 600);
    c->SetGridy();
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(20);
    h->SetMinimum(0.0);
    h->SetMaximum(1.1);
    h->Draw("P E");

    TLine* line = new TLine(binEdges.front(), 0.5, binEdges.back(), 0.5);
    line->SetLineStyle(2);
    line->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.15, 0.84, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.15, 0.79, suffix == "Purity" ? "Diag / Reco_in_bin" : "Diag / Gen_in_bin");

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  drawAndSave(hPurity, "Purity", "Purity", kBlue+1);
  drawAndSave(hStability, "Stability", "Stability", kRed+1);

  std::cout << "[Info] Saved PNGs: " << outBase << "_Purity/Stability.png" << std::endl;
  fIn->Close();
}
