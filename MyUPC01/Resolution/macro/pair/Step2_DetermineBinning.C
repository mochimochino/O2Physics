// ============================================================
//  Step2_DetermineBinning.C
// ============================================================
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ==========================================================
// Top-down binning
// ==========================================================
std::vector<double> AutoDetermineBinningTopDown(TH2F* hMatrixFine, double targetPurity, double targetStability, double minStats, double maxPt2)
{
  std::vector<double> edges;
  edges.push_back(maxPt2);

  int maxBinFine = hMatrixFine->GetXaxis()->FindBin(maxPt2 - 1e-6);
  int upperBin = maxBinFine;

  for (int lowerBin = maxBinFine; lowerBin >= 1; --lowerBin) {
    if (lowerBin % 100 == 0 || lowerBin == 1) {
      int progress = 100 - (int)(100.0 * lowerBin / maxBinFine);
      std::cout << "\r[1/4] Auto Binning... " << progress << "% completed" << std::flush;
    }

    double diag = hMatrixFine->Integral(lowerBin, upperBin, lowerBin, upperBin);
    double sumTrue = hMatrixFine->Integral(lowerBin, upperBin, 1, hMatrixFine->GetNbinsY());
    double sumReco = hMatrixFine->Integral(1, hMatrixFine->GetNbinsX(), lowerBin, upperBin);

    if (sumTrue <= 0 || sumReco <= 0)
      continue;

    double purity = diag / sumReco;
    double stability = diag / sumTrue;

    if (purity >= targetPurity && stability >= targetStability &&
        sumReco >= minStats && sumTrue >= minStats) {
      double lowerEdge = hMatrixFine->GetXaxis()->GetBinLowEdge(lowerBin);
      if (lowerBin > 1 && lowerEdge > 0.0) {
        edges.push_back(lowerEdge);
        upperBin = lowerBin - 1;
      }
    }
  }
  std::cout << "\r[1/4] Auto Binning... 100% completed!       " << std::endl;

  if (edges.back() != 0.0)
    edges.push_back(0.0);
  std::reverse(edges.begin(), edges.end());
  return edges;
}

void Step2_DetermineBinning()
{
  // ----------------------------------------------------------
  // Settings
  // ----------------------------------------------------------
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/CoherenrtJpsi/";
  const TString inFile = inDir + "Step1_merged.root";
  const TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0413/Coherent/";
  const TString realDataInDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/CoherenrtJpsi/Step1_merged.root";

  gStyle->SetOptStat(0);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Error: Cannot open " << inFile << std::endl;
    return;
  }
  TH2F* hMatrixFine = (TH2F*)fIn->Get("hResponseMatrixPairPt2");
  if (!hMatrixFine) {
    std::cerr << "Error: hResponseMatrixPairPt2 not found." << std::endl;
    return;
  }

  double targetPurity = 0.45;
  double targetStability = 0.45;
  double maxPt2 = 0.065;
  double minStats = 1000;

  std::vector<double> pt2Bins = AutoDetermineBinningTopDown(hMatrixFine, targetPurity, targetStability, minStats, maxPt2);
  const int nBins = pt2Bins.size() - 1;

  // ==========================================================
  // Resolution evaluation
  // ==========================================================
  TGraphErrors* grResolution = new TGraphErrors();
  int nSlices = 20;
  for (int i = 0; i < nSlices; ++i) {
    double pt2Min = i * (maxPt2 / nSlices);
    double pt2Max = (i + 1) * (maxPt2 / nSlices);
    int binMin = hMatrixFine->GetXaxis()->FindBin(pt2Min + 1e-6);
    int binMax = hMatrixFine->GetXaxis()->FindBin(pt2Max - 1e-6);
    TH1D* hSlice = hMatrixFine->ProjectionY(Form("slice_%d", i), binMin, binMax);
    if (hSlice->GetEntries() > 50) {
      hSlice->Fit("gaus", "Q0");
      TF1* fitFunc = hSlice->GetFunction("gaus");
      if (fitFunc) {
        int np = grResolution->GetN();
        grResolution->SetPoint(np, (pt2Min + pt2Max) / 2.0, fitFunc->GetParameter(2));
        grResolution->SetPointError(np, 0.0, fitFunc->GetParError(2));
      }
    }
    delete hSlice;
  }

  // ==========================================================
  // Rebinning Matrix
  // ==========================================================
  TH2D* hMatrixRebinned = new TH2D("hMatrixRebinned", "Response Matrix;Gen p_{T}^{2};Reco p_{T}^{2}", nBins, pt2Bins.data(), nBins, pt2Bins.data());
  for (int ix = 1; ix <= hMatrixFine->GetNbinsX(); ++ix) {
    double truePt2 = hMatrixFine->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= hMatrixFine->GetNbinsY(); ++iy) {
      hMatrixRebinned->Fill(truePt2, hMatrixFine->GetYaxis()->GetBinCenter(iy), hMatrixFine->GetBinContent(ix, iy));
    }
  }

  // ==========================================================
  // Purity and Stability Calculation
  // ==========================================================
  TH1D* hPurity = new TH1D("hPurity", "Purity", nBins, pt2Bins.data());
  TH1D* hStability = new TH1D("hStability", "Stability", nBins, pt2Bins.data());
  for (int i = 1; i <= nBins; ++i) {
    double diag = hMatrixRebinned->GetBinContent(i, i);
    double sumReco = hMatrixRebinned->Integral(1, nBins, i, i);
    double sumTrue = hMatrixRebinned->Integral(i, i, 1, nBins);
    hPurity->SetBinContent(i, (sumReco > 0) ? diag / sumReco : 0);
    hStability->SetBinContent(i, (sumTrue > 0) ? diag / sumTrue : 0);
  }

  // ==========================================================
  // Rebin 1D Histograms
  // ==========================================================
  TH1D* hGenFine = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hRecoFine = (TH1D*)fIn->Get("hRecoPt2");
  TH1D* hGenRebinned = (hGenFine) ? (TH1D*)hGenFine->Rebin(nBins, "hMCGen_Rebinned", pt2Bins.data()) : nullptr;
  TH1D* hRecoRebinned = (hRecoFine) ? (TH1D*)hRecoFine->Rebin(nBins, "hMCReco_Rebinned", pt2Bins.data()) : nullptr;

  TFile* fData = TFile::Open(realDataInDir, "READ");
  TH1D* hDataRebinned = nullptr;
  if (fData && !fData->IsZombie()) {
    TH1D* hDataFine = (TH1D*)fData->Get("hRecoPt2");
    if (hDataFine)
      hDataRebinned = (TH1D*)hDataFine->Rebin(nBins, "hDataReco_Rebinned", pt2Bins.data());
  }

  // ==========================================================
  // Final Output and Validation (ERROR CHECKING)
  // ==========================================================
  std::cout << "\n============================================================" << std::endl;
  std::cout << " Binning Summary & Validation" << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << Form("%-5s | %-15s | %-8s | %-10s | %-8s | %-8s | %-10s | %-12s", "Bin", "pT2 Range", "Events", "StatErr%", "Purity%", "Stab%", "Reco/Gen%", "(Diag+Gen)/Gen%") << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;

  bool hasCriticalError = false;
  std::vector<TString> errorMessages;

  for (int i = 1; i <= nBins; ++i) {
    double events = (hDataRebinned) ? hDataRebinned->GetBinContent(i) : (hRecoRebinned ? hRecoRebinned->GetBinContent(i) : 0);
    double statErr = (events > 0) ? (100.0 / std::sqrt(events)) : 999.0;
    double pur = hPurity->GetBinContent(i) * 100.0;
    double stab = hStability->GetBinContent(i) * 100.0;
    double recoGenRatio = (hGenRebinned && hRecoRebinned && hGenRebinned->GetBinContent(i) > 0) ? (hRecoRebinned->GetBinContent(i) / hGenRebinned->GetBinContent(i)) * 100.0 : 0.0;

    double diag = hMatrixRebinned->GetBinContent(i, i);
    double diagGenRatio = (hGenRebinned && hGenRebinned->GetBinContent(i) > 0) ? ((diag + hGenRebinned->GetBinContent(i)) / hGenRebinned->GetBinContent(i)) * 100.0 : 0.0;

    std::cout << Form("%-5d | %5.4f-%5.4f | %-8.0f | %-10.2f | %-8.1f | %-8.1f | %-10.1f | %-12.1f", i, pt2Bins[i - 1], pt2Bins[i], events, statErr, pur, stab, recoGenRatio, diagGenRatio) << std::endl;

    // --- Validation Checks ---
    if (statErr > 20.0) {
      errorMessages.push_back(Form("[CRITICAL] Bin %d: Statistical error is too high (%.1f%% > 20%%).", i, statErr));
      hasCriticalError = true;
    }
    if (pur < targetPurity * 100.0 || stab < targetStability * 100.0) {
      errorMessages.push_back(Form("[WARNING] Bin %d: Purity or Stability is below target threshold.", i));
    }
  }

  std::cout << "============================================================\n"
            << std::endl;

  if (hasCriticalError || !errorMessages.empty()) {
    std::cerr << "!!! ANALYSIS ALERTS !!!" << std::endl;
    for (const auto& msg : errorMessages) {
      std::cerr << "  " << msg << std::endl;
    }
    if (hasCriticalError) {
      std::cerr << "\n[!] ERROR: Some bins have insufficient statistics. Unfolding may be unstable." << std::endl;
    }
  } else {
    std::cout << "[SUCCESS] All bins passed validation criteria." << std::endl;
  }

  // Save to file
  TFile* fOut = TFile::Open(outDir + "Step2_Response_for_Unfolding.root", "RECREATE");
  hMatrixRebinned->Write("hResponseMatrix");
  hPurity->Write();
  hStability->Write();
  if (hGenRebinned)
    hGenRebinned->Write();
  if (hRecoRebinned)
    hRecoRebinned->Write();
  if (hDataRebinned)
    hDataRebinned->Write();
  fOut->Close();

  std::cout << "\n[Done] Output saved to: " << outDir << "Step2_Response_for_Unfolding.root" << std::endl;
}
