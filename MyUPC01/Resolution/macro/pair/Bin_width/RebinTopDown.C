// ============================================================
//  RebinTopDown.C
// ============================================================

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <iostream>
#include <vector>

// -------- Make2DHistFromBinnedHists --------
#include "../Make2DHistFromBinnedHists.C"

// ----------------------------------------------------------
//  TopDown
// ----------------------------------------------------------
std::vector<double> TopDownBinEdges(
  TH2F* hMatrixFine,
  double targetPurity,
  double targetStability,
  double minStats,
  double maxPt2)
{
  std::vector<double> edges;
  edges.push_back(maxPt2);

  int maxBinFine = hMatrixFine->GetXaxis()->FindBin(maxPt2 - 1e-6);
  int upperBin = maxBinFine;

  for (int lowerBin = maxBinFine; lowerBin >= 1; --lowerBin) {
    if (lowerBin % 100 == 0 || lowerBin == 1) {
      int progress = 100 - (int)(100.0 * lowerBin / maxBinFine);
      std::cout << "\r[TopDown] Scanning bins... " << progress << "%" << std::flush;
    }

    double diag = hMatrixFine->Integral(lowerBin, upperBin, lowerBin, upperBin);
    double sumTrue = hMatrixFine->Integral(lowerBin, upperBin, 0, hMatrixFine->GetNbinsY() + 1);
    double sumReco = hMatrixFine->Integral(0, hMatrixFine->GetNbinsX() + 1, lowerBin, upperBin);

    if (sumTrue <= 0 || sumReco <= 0)
      continue;

    double purity = diag / sumReco;
    double stability = diag / sumTrue;

    if (purity >= targetPurity && stability >= targetStability &&
        sumReco >= minStats && sumTrue >= minStats) {
      double lowerEdge = hMatrixFine->GetXaxis()->GetBinLowEdge(lowerBin);

      edges.push_back(lowerEdge);
      upperBin = lowerBin - 1;
    }
  }
  std::cout << "\r[TopDown] Scanning bins... 100% done.       " << std::endl;

  if (upperBin >= 1) {
    if (edges.size() > 1) {
      std::cout << "[TopDown] Final segment [0.0, " << edges.back() << "] failed criteria." << std::endl;
      std::cout << "          -> Merging it with the adjacent upper bin." << std::endl;
      edges.pop_back();
    } else {
      std::cout << "[TopDown] Warning: No valid bins could be formed before reaching 0." << std::endl;
    }
  }

  if (edges.back() != 0.0) {
    edges.push_back(0.0);
  }

  std::reverse(edges.begin(), edges.end());
  return edges;
}

// ----------------------------------------------------------
//  Main
// ----------------------------------------------------------
void RebinTopDown(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0427test/Incoherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0427test/Incoherent/topdown/40/Step2_Rebinned.root",
  double targetPurity = 0.40,
  double targetStability = 0.40,
  double minStats = 500,
  double maxPt2 = 1.5)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open: " << inputFile << std::endl;
    return;
  }
  TH2F* hMatrixFine = (TH2F*)fIn->Get("hResponseMatrixPairPt2");
  if (!hMatrixFine) {
    std::cerr << "[Error] hResponseMatrixPairPt2 not found." << std::endl;
    fIn->Close();
    return;
  }

  std::cout << "\n[Info] Purity threshold    : " << targetPurity << std::endl;
  std::cout << "[Info] Stability threshold : " << targetStability << std::endl;
  std::cout << "[Info] Min statistics      : " << minStats << std::endl;
  std::cout << "[Info] Max pT2             : " << maxPt2 << std::endl;

  std::vector<double> bins = TopDownBinEdges(hMatrixFine, targetPurity, targetStability, minStats, maxPt2);
  const int nBins = (int)bins.size() - 1;

  if (nBins <= 0) {
    std::cerr << "[Error] No bins were determined. "
              << "Try loosening targetPurity / targetStability / minStats." << std::endl;
    fIn->Close();
    return;
  }

  std::cout << "\n[Info] Determined " << nBins << " bins:" << std::endl;
  for (int i = 0; i < (int)bins.size(); ++i)
    std::cout << "  edge[" << i << "] = " << bins[i] << std::endl;

  // Gen / Reco 1D rebin
  TH1D* hGenFine = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hRecoFine = (TH1D*)fIn->Get("hRecoPt2");
  if (!hGenFine || !hRecoFine) {
    std::cerr << "[Error] hGenPt2 or hRecoPt2 not found." << std::endl;
    fIn->Close();
    return;
  }

  TH1D* hGenRebin = (TH1D*)hGenFine->Rebin(nBins, "hGenPt2_rebin", bins.data());
  TH1D* hRecoRebin = (TH1D*)hRecoFine->Rebin(nBins, "hRecoPt2_rebin", bins.data());
  hGenRebin->SetDirectory(nullptr);
  hRecoRebin->SetDirectory(nullptr);
  fIn->Close();

  // Save ROOT file
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  hGenRebin->Write();
  hRecoRebin->Write();
  fOut->Close();
  std::cout << "[Info] Saved ROOT: " << outputFile << std::endl;

  TString outBase = outputFile;
  outBase.ReplaceAll(".root", "");

  auto drawAndSave = [&](TH1D* h, const TString& title, const TString& suffix, int color) {
    TCanvas* c = new TCanvas("c_" + suffix, title, 800, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
    c->SetGrid();
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->Draw("HIST");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.04);
    // tex.DrawLatex(0.35, 0.84, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.35, 0.79, "Method: TopDown (Purity+Stability)");
    tex.DrawLatex(0.35, 0.74, Form("Pur >= %.2f, Stab >= %.2f, minN >= %.0f", targetPurity, targetStability, minStats));
    tex.DrawLatex(0.35, 0.69, Form("Bins: %d", nBins));

    for (int i = 1; i < (int)bins.size() - 1; ++i) {
      TLine* vl = new TLine(bins[i], h->GetMinimum(), bins[i], h->GetMaximum() * 0.9);
      vl->SetLineStyle(2);
      vl->SetLineColor(kGray + 1);
      vl->Draw();
    }

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  drawAndSave(hGenRebin, "Gen p_{T}^{2} (TopDown)", "Gen", kAzure + 2);
  drawAndSave(hRecoRebin, "Reco p_{T}^{2} (TopDown)", "Reco", kOrange + 7);
  std::cout << "[Info] Saved PNGs: " << outBase << "_Gen/Reco.png" << std::endl;

  // Make 2D histogram
  const TString matrixFile = inputFile;
  const TString matrixName = "hResponseMatrixPairPt2";

  Make2DHistFromBinnedHists(
    matrixFile,
    matrixName,
    outputFile,
    "hGenPt2_rebin",
    "hRecoPt2_rebin",
    outputFile);
}
