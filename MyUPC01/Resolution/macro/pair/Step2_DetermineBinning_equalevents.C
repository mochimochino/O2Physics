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
// Equal Events Binning (トップダウンを廃止し、Bin数を直接指定する方式)
// ==========================================================
std::vector<double> DetermineEqualEventsBinning(TH1D* hInput, int nTargetBins, double maxPt2)
{
  std::vector<double> edges;
  edges.push_back(0.0);

  int maxBin = hInput->GetXaxis()->FindBin(maxPt2 - 1e-6);
  double totalIntegral = hInput->Integral(1, maxBin);
  double eventsPerBin = totalIntegral / nTargetBins;

  std::cout << "[Info] Total events in range [0, " << maxPt2 << "]: " << totalIntegral << std::endl;
  std::cout << "[Info] Target bins: " << nTargetBins << " (Target events per bin: " << eventsPerBin << ")" << std::endl;

  double currentSum = 0;
  int binsFound = 0;

  for (int i = 1; i <= maxBin; ++i) {
    currentSum += hInput->GetBinContent(i);
    // 最後のBinは残りの全イベントを収容するため条件から外す
    if (currentSum >= eventsPerBin && binsFound < nTargetBins - 1) {
      double edge = hInput->GetXaxis()->GetBinUpEdge(i);
      edges.push_back(edge);
      currentSum = 0;
      binsFound++;
    }
  }

  // 終端を確実にmaxPt2に設定
  if (edges.back() < maxPt2) {
    edges.push_back(maxPt2);
  } else {
    edges.back() = maxPt2;
  }

  return edges;
}

void Step2_DetermineBinning_equalevents()
{
  // ----------------------------------------------------------
  // Settings
  // ----------------------------------------------------------
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/CoherenrtJpsi/";
  const TString inFile = inDir + "Step1_merged.root";
  const TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0413/Coherent/equal/";
  const TString realDataInDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/CoherenrtJpsi/Step1_merged.root";

  // === ここでBinの数を指定します ===
  int nTargetBins = 6;
  double maxPt2 = 0.065;
  bool useDataForBinning = false;

  gStyle->SetOptStat(0);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "Error: Cannot open " << inFile << std::endl;
    return;
  }
  TH2F* hMatrixFine = (TH2F*)fIn->Get("hResponseMatrixPairPt2");
  TH1D* hGenFine = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hRecoFine = (TH1D*)fIn->Get("hRecoPt2");

  if (!hMatrixFine || !hGenFine || !hRecoFine) {
    std::cerr << "Error: Necessary histograms not found in " << inFile << std::endl;
    return;
  }

  TFile* fData = TFile::Open(realDataInDir, "READ");
  TH1D* hDataFine = nullptr;
  if (fData && !fData->IsZombie()) {
    hDataFine = (TH1D*)fData->Get("hRecoPt2");
  }

  // ----------------------------------------------------------
  // Binningの決定 (指定されたBin数で実行)
  // ----------------------------------------------------------
  TH1D* hForBinning = (useDataForBinning && hDataFine) ? hDataFine : hRecoFine;
  std::vector<double> pt2Bins = DetermineEqualEventsBinning(hForBinning, nTargetBins, maxPt2);
  const int nBins = pt2Bins.size() - 1;

  // ==========================================================
  // Rebinning Matrix & 1D Histograms
  // ==========================================================
  TH2D* hMatrixRebinned = new TH2D("hMatrixRebinned", "Response Matrix;Gen p_{T}^{2};Reco p_{T}^{2}", nBins, pt2Bins.data(), nBins, pt2Bins.data());
  for (int ix = 1; ix <= hMatrixFine->GetNbinsX(); ++ix) {
    double truePt2 = hMatrixFine->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= hMatrixFine->GetNbinsY(); ++iy) {
      hMatrixRebinned->Fill(truePt2, hMatrixFine->GetYaxis()->GetBinCenter(iy), hMatrixFine->GetBinContent(ix, iy));
    }
  }

  TH1D* hGenRebinned = (hGenFine) ? (TH1D*)hGenFine->Rebin(nBins, "hMCGen_Rebinned", pt2Bins.data()) : nullptr;
  TH1D* hRecoRebinned = (hRecoFine) ? (TH1D*)hRecoFine->Rebin(nBins, "hMCReco_Rebinned", pt2Bins.data()) : nullptr;
  TH1D* hDataRebinned = (hDataFine) ? (TH1D*)hDataFine->Rebin(nBins, "hDataReco_Rebinned", pt2Bins.data()) : nullptr;

  // ==========================================================
  // 指標の計算 (Purity, Stability, (Rec+Gen)/Gen)
  // ==========================================================
  TH1D* hPurity = new TH1D("hPurity", "Purity [%]", nBins, pt2Bins.data());
  TH1D* hStability = new TH1D("hStability", "Stability [%]", nBins, pt2Bins.data());
  TH1D* hRecPlusGenOverGen = new TH1D("hRecPlusGenOverGen", "(Reco+Gen)/Gen [%]", nBins, pt2Bins.data());

  for (int i = 1; i <= nBins; ++i) {
    double diag = hMatrixRebinned->GetBinContent(i, i);
    double sumReco = hMatrixRebinned->Integral(1, nBins, i, i);
    double sumTrue = hMatrixRebinned->Integral(i, i, 1, nBins);

    // Purity と Stability をパーセント(%)で格納
    hPurity->SetBinContent(i, (sumReco > 0) ? (diag / sumReco) * 100.0 : 0);
    hStability->SetBinContent(i, (sumTrue > 0) ? (diag / sumTrue) * 100.0 : 0);

    // ご要望の (Reco + Gen) / Gen をパーセント(%)で計算して格納
    if (hGenRebinned && hRecoRebinned && hGenRebinned->GetBinContent(i) > 0) {
      double gen = hGenRebinned->GetBinContent(i);
      double reco = hRecoRebinned->GetBinContent(i);
      double ratioPct = ((reco + gen) / gen) * 100.0;
      hRecPlusGenOverGen->SetBinContent(i, ratioPct);
    }
  }

  // ==========================================================
  // ターミナル出力
  // ==========================================================
  std::cout << "\n============================================================" << std::endl;
  std::cout << " Binning Summary (" << nBins << " Bins Specified)" << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << Form("%-5s | %-15s | %-8s | %-10s | %-8s | %-8s | %-15s", "Bin", "pT2 Range", "Events", "StatErr%", "Purity%", "Stab%", "(Rec+Gen)/Gen%") << std::endl;
  std::cout << "---------------------------------------------------------------------------------------" << std::endl;

  for (int i = 1; i <= nBins; ++i) {
    double events = (hDataRebinned) ? hDataRebinned->GetBinContent(i) : (hRecoRebinned ? hRecoRebinned->GetBinContent(i) : 0);
    double statErr = (events > 0) ? (100.0 / std::sqrt(events)) : 999.0;

    double pur = hPurity->GetBinContent(i);
    double stab = hStability->GetBinContent(i);
    double rPlusG_overG = hRecPlusGenOverGen->GetBinContent(i);

    std::cout << Form("%-5d | %5.4f-%5.4f | %-8.0f | %-10.2f | %-8.1f | %-8.1f | %-15.1f",
                      i, pt2Bins[i - 1], pt2Bins[i], events, statErr, pur, stab, rPlusG_overG)
              << std::endl;
  }
  std::cout << "============================================================\n"
            << std::endl;

  // ==========================================================
  // ROOTファイルへの保存
  // ==========================================================
  TFile* fOut = TFile::Open(outDir + "Step2_Response_for_Unfolding.root", "RECREATE");
  hMatrixRebinned->Write("hResponseMatrix");

  // %表記で値が入ったヒストグラムを書き込み
  hPurity->Write();
  hStability->Write();
  hRecPlusGenOverGen->Write();

  if (hGenRebinned)
    hGenRebinned->Write();
  if (hRecoRebinned)
    hRecoRebinned->Write();
  if (hDataRebinned)
    hDataRebinned->Write();
  fOut->Close();

  std::cout << "[Done] Output saved to: " << outDir << "Step2_Response_for_Unfolding.root" << std::endl;
}
