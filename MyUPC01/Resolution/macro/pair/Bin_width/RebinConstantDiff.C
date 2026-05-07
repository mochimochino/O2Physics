#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <vector>

// -------- Make2DHistFromBinnedHists との連携 --------
#include "../Make2DHistFromBinnedHists.C"

void RebinConstantDiff(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Incoherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Incoherent/constantdiff/15bin/Step2_Rebinned.root",
  double initialWidth = 0.016,
  double deltaW = 0.012,
  double xMin = 0.0,
  double xMax = 1.2,
  bool respectOverflow = true)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  std::vector<double> bins;
  bins.push_back(xMin);
  double currentX = xMin;
  double currentWidth = initialWidth;
  while (currentX < xMax) {
    currentX += currentWidth;
    if (currentX > xMax && !respectOverflow)
      break;
    bins.push_back(currentX);
    currentWidth += deltaW;
  }
  if (bins.back() != xMax && respectOverflow)
    bins.back() = xMax;

  int nBins = bins.size() - 1;
  double* binArray = &bins[0];

  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie())
    return;
  TH1D* hGen = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hReco = (TH1D*)fIn->Get("hRecoPt2");
  if (!hGen || !hReco)
    return;

  TH1D* hGenRebin = (TH1D*)hGen->Rebin(nBins, "hGenPt2_rebin", binArray);
  TH1D* hRecoRebin = (TH1D*)hReco->Rebin(nBins, "hRecoPt2_rebin", binArray);

  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  hGenRebin->Write();
  hRecoRebin->Write();
  fOut->Close();
  std::cout << "[Info] Saved ROOT: " << outputFile << std::endl;

  TString outBase = outputFile;
  outBase.ReplaceAll(".root", "");

  auto drawAndSave = [&](TH1D* h, const TString& title, const TString& suffix) {
    TCanvas* c = new TCanvas("c_" + suffix, title, 800, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
    c->SetGrid();
    h->SetLineColor(suffix == "Gen" ? kGreen + 2 : kOrange + 7);
    h->SetLineWidth(2);
    h->Draw("HIST E");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.03);
    // tex.DrawLatex(0.55, 0.84, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.55, 0.79, "Method: Constant Bin Difference");
    tex.DrawLatex(0.55, 0.74, Form("w0: %.3f, dw: %.3f", initialWidth, deltaW));

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  drawAndSave(hGenRebin, "Gen p_{T}^{2} (Const Diff)", "Gen");
  drawAndSave(hRecoRebin, "Reco p_{T}^{2} (Const Diff)", "Reco");

  std::cout << "[Info] Saved PNGs: " << outBase << "_Gen/Reco.png" << std::endl;
  fIn->Close();

  // 2D hist
  const TString matrixFile = inputFile;
  const TString matrixName = "hResponseMatrixPairPt2";
  const TString outFile2D = outputFile;

  Make2DHistFromBinnedHists(
    matrixFile,
    matrixName,
    outputFile,
    "hGenPt2_rebin",
    "hRecoPt2_rebin",
    outFile2D);
}
