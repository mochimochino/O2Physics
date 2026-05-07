#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// -------- Make2DHistFromBinnedHists との連携 --------
#include "../Make2DHistFromBinnedHists.C"

void RebinFlatStats(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Sq/flatstats/4bin/Step2_Rebinned.root",
  int nTargetBins = 4,
  double xMin = 0.0,
  double xMax = 0.065)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie())
    return;
  TH1D* hGen = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hReco = (TH1D*)fIn->Get("hRecoPt2");
  if (!hGen || !hReco)
    return;

  int binMin = hGen->FindBin(xMin);
  int binMax = hGen->FindBin(xMax);

  std::vector<double> cumSum;
  cumSum.push_back(0.0);
  for (int i = binMin; i <= binMax; ++i) {
    cumSum.push_back(cumSum.back() + hGen->GetBinContent(i));
  }

  double totalIntegral = cumSum.back();
  double step = totalIntegral / nTargetBins;

  std::vector<double> bins;
  bins.push_back(xMin);

  for (int k = 1; k < nTargetBins; ++k) {
    double target = k * step;
    int bestIdx = 1;
    double minDiff = std::abs(cumSum[1] - target);

    for (size_t i = 2; i < cumSum.size(); ++i) {
      double diff = std::abs(cumSum[i] - target);
      if (diff < minDiff) {
        minDiff = diff;
        bestIdx = i;
      }
    }

    double edge = hGen->GetBinLowEdge(binMin + bestIdx);
    if (edge > bins.back() && edge < xMax) {
      bins.push_back(edge);
    }
  }
  if (bins.back() < xMax) {
    bins.push_back(xMax);
  }

  int nBins = bins.size() - 1;
  double* binArray = &bins[0];

  TH1D* hGenRebin = (TH1D*)hGen->Rebin(nBins, "hGenPt2_rebin", binArray);
  TH1D* hRecoRebin = (TH1D*)hReco->Rebin(nBins, "hRecoPt2_rebin", binArray);

  std::cout << "========================================" << std::endl;
  std::cout << "[Info] Target events per bin: " << step << std::endl;
  std::cout << "[Info] Gen Bin Contents after rebinning:" << std::endl;
  for (int i = 1; i <= nBins; ++i) {
    std::cout << "  Bin " << i << " [" << hGenRebin->GetBinLowEdge(i) << ", "
              << hGenRebin->GetBinLowEdge(i + 1) << "] : "
              << hGenRebin->GetBinContent(i) << " events" << std::endl;
  }
  std::cout << "========================================" << std::endl;

  // --- Save to ROOT ---
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  hGenRebin->Write();
  hRecoRebin->Write();
  fOut->Close();
  std::cout << "[Info] Saved ROOT: " << outputFile << std::endl;

  // --- Draw and Save PNG ---
  TString outBase = outputFile;
  outBase.ReplaceAll(".root", "");

  auto drawAndSave = [&](TH1D* h, const TString& title, const TString& suffix) {
    TCanvas* c = new TCanvas("c_" + suffix, title, 800, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
    c->SetGrid();
    h->SetMinimum(0.0);
    h->SetMaximum(h->GetMaximum() * 1.2);
    h->SetLineColor(suffix == "Gen" ? kRed + 1 : kBlue + 1);
    h->SetLineWidth(2);

    h->Draw("HIST");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.03);
    // tex.DrawLatex(0.55, 0.84, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.55, 0.79, "Method: Flat Distributions");
    tex.DrawLatex(0.55, 0.74, Form("Bins: %d", nBins));

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  drawAndSave(hGenRebin, "Gen p_{T}^{2} (Flat Stats)", "Gen");
  drawAndSave(hRecoRebin, "Reco p_{T}^{2} (Flat Stats)", "Reco");

  std::cout << "[Info] Saved PNGs: " << outBase << "_Gen.png / " << outBase << "_Reco.png" << std::endl;
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
