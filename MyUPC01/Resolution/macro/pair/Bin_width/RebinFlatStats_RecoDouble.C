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

#include "../Make2DHistFromBinnedHists.C"

// Compute flat-statistics bin edges from a histogram over [xMin, xMax]
std::vector<double> ComputeFlatStatsBins_RD(TH1D* h, int nBins, double xMin, double xMax)
{
  int binMin = h->FindBin(xMin);
  int binMax = h->FindBin(xMax);

  std::vector<double> cumSum;
  cumSum.push_back(0.0);
  for (int i = binMin; i <= binMax; ++i)
    cumSum.push_back(cumSum.back() + h->GetBinContent(i));

  double totalIntegral = cumSum.back();
  double step = totalIntegral / nBins;

  std::vector<double> bins;
  bins.push_back(xMin);

  for (int k = 1; k < nBins; ++k) {
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
    double edge = h->GetBinLowEdge(binMin + bestIdx);
    if (edge > bins.back() && edge < xMax)
      bins.push_back(edge);
  }
  if (bins.back() < xMax)
    bins.push_back(xMax);

  return bins;
}

void RebinFlatStats_RecoDouble(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/15697/Coherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/15697/Coherent/flatstats/5bin/Step2_Rebinned_RecoDouble.root",
  int nTargetBins = 5,
  double xMin = 0.0,
  double xMax = 0.06)
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

  // Gen: nTargetBins bins, Reco: 2*nTargetBins bins — both derived from hGen
  std::vector<double> genBins = ComputeFlatStatsBins_RD(hGen, nTargetBins, xMin, xMax);
  std::vector<double> recoBins = ComputeFlatStatsBins_RD(hGen, 2 * nTargetBins, xMin, xMax);

  int nGenBins = (int)genBins.size() - 1;
  int nRecoBins = (int)recoBins.size() - 1;

  TH1D* hGenRebin = (TH1D*)hGen->Rebin(nGenBins, "hGenPt2_rebin", genBins.data());
  TH1D* hRecoRebin = (TH1D*)hReco->Rebin(nRecoBins, "hRecoPt2_rebin", recoBins.data());

  std::cout << "========================================" << std::endl;
  std::cout << "[Info] Gen bins: " << nGenBins << "  Reco bins: " << nRecoBins << std::endl;
  std::cout << "[Info] Gen Bin Contents after rebinning:" << std::endl;
  for (int i = 1; i <= nGenBins; ++i) {
    std::cout << "  Gen Bin " << i << " [" << hGenRebin->GetBinLowEdge(i) << ", "
              << hGenRebin->GetBinLowEdge(i + 1) << "] : "
              << hGenRebin->GetBinContent(i) << " events" << std::endl;
  }
  std::cout << "[Info] Reco Bin Contents after rebinning:" << std::endl;
  for (int i = 1; i <= nRecoBins; ++i) {
    std::cout << "  Reco Bin " << i << " [" << hRecoRebin->GetBinLowEdge(i) << ", "
              << hRecoRebin->GetBinLowEdge(i + 1) << "] : "
              << hRecoRebin->GetBinContent(i) << " events" << std::endl;
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

  auto drawAndSave = [&](TH1D* h, const TString& title, const TString& suffix, int nBinsInfo) {
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
    tex.DrawLatex(0.55, 0.79, "Method: Flat Distributions");
    tex.DrawLatex(0.55, 0.74, Form("Bins: %d", nBinsInfo));

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  drawAndSave(hGenRebin, "Gen p_{T}^{2} (Flat Stats)", "Gen", nGenBins);
  drawAndSave(hRecoRebin, "Reco p_{T}^{2} (Flat Stats x2)", "Reco", nRecoBins);

  std::cout << "[Info] Saved PNGs: " << outBase << "_Gen.png / " << outBase << "_Reco.png" << std::endl;
  fIn->Close();

  // 2D response matrix (non-square: nGenBins x nRecoBins = N x 2N)
  Make2DHistFromBinnedHists(
    inputFile,
    "hResponseMatrixPairPt2",
    outputFile,
    "hGenPt2_rebin",
    "hRecoPt2_rebin",
    outputFile);
}
