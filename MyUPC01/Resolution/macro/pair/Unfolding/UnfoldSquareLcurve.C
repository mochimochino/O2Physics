// ===========================================================
// UnfoldSquareLcurve.C
// TUnfold with L-curve Optimization using Square Matrix
// ==========================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TUnfoldDensity.h"

#include <iostream>

// ヒストグラム読み込み用ヘルパー関数
template <typename T>
T* LoadHistogram(TFile* file, const TString& histName)
{
  T* h = (T*)file->Get(histName);
  if (!h) {
    std::cerr << "[ERROR] Histogram not found: " << histName << std::endl;
    return nullptr;
  }
  T* hClone = (T*)h->Clone(histName + "_clone");
  hClone->SetDirectory(nullptr);
  return hClone;
}

void UnfoldSquareLcurve()
{
  // ===========================
  // 1. Settings
  // ===========================
  const TString mcFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Sq/flatstats/4bin/Step2_Rebinned.root";
  const TString histNameMatrix = "hResponseMatrix";

  const TString dataFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Step1_merged.root";
  const TString histNameDataFine = "hRecoPt2_Test";
  const TString histNameGenFine = "hGenPt2_Test";

  gStyle->SetOptStat(0);

  // ===========================
  // 2. Load Real Data & MC
  // ===========================
  TFile* fMC = TFile::Open(mcFile, "READ");
  if (!fMC || fMC->IsZombie())
    return;
  TH2D* hMat = LoadHistogram<TH2D>(fMC, histNameMatrix);
  fMC->Close();

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie())
    return;
  TH1D* hDataFine = LoadHistogram<TH1D>(fData, histNameDataFine);
  TH1D* hGenTestFine = LoadHistogram<TH1D>(fData, histNameGenFine);
  fData->Close();

  if (!hMat || !hDataFine || !hGenTestFine)
    return;

  // ===========================
  // 3. Matrix Check & U/O Bin Clearing
  // ===========================
  int nBinsGen = hMat->GetNbinsX();
  int nBinsReco = hMat->GetNbinsY();

  if (nBinsGen != nBinsReco) {
    std::cerr << "[WARNING] The response matrix is not square! Gen Bins: "
              << nBinsGen << ", Reco Bins: " << nBinsReco << std::endl;
  }

  // 【最重要】 hMat の Underflow / Overflow ビンを完全にゼロクリアする
  // TUnfoldは0のビンを自動的に除外するため、Nx=7がNx=5になります
  for (int ix = 0; ix <= nBinsGen + 1; ++ix) {
    for (int iy = 0; iy <= nBinsReco + 1; ++iy) {
      if (ix == 0 || ix == nBinsGen + 1 || iy == 0 || iy == nBinsReco + 1) {
        hMat->SetBinContent(ix, iy, 0.0);
        hMat->SetBinError(ix, iy, 0.0);
      }
    }
  }

  // エッジ取得とテストデータの詰め替え
  double* edgesGen = new double[nBinsGen + 1];
  for (int i = 1; i <= nBinsGen + 1; ++i) {
    edgesGen[i - 1] = hMat->GetXaxis()->GetBinLowEdge(i);
  }

  double* edgesReco = new double[nBinsReco + 1];
  for (int i = 1; i <= nBinsReco + 1; ++i) {
    edgesReco[i - 1] = hMat->GetYaxis()->GetBinLowEdge(i);
  }

  TH1D* hDataTest = new TH1D("hDataTest_Rebinned", "Data Test Rebinned", nBinsReco, edgesReco);
  for (int i = 1; i <= hDataFine->GetNbinsX(); ++i) {
    hDataTest->Fill(hDataFine->GetBinCenter(i), hDataFine->GetBinContent(i));
  }
  // Data側の U/O もゼロクリア
  hDataTest->SetBinContent(0, 0.0);
  hDataTest->SetBinError(0, 0.0);
  hDataTest->SetBinContent(nBinsReco + 1, 0.0);
  hDataTest->SetBinError(nBinsReco + 1, 0.0);

  TH1D* hGenTest = new TH1D("hGenTest_Rebinned", "Gen Test Rebinned", nBinsGen, edgesGen);
  for (int i = 1; i <= hGenTestFine->GetNbinsX(); ++i) {
    hGenTest->Fill(hGenTestFine->GetBinCenter(i), hGenTestFine->GetBinContent(i));
  }

  delete[] edgesReco;
  delete[] edgesGen;

  // ===========================
  // 4. Setup TUnfold
  // ===========================
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;

  // TUnfoldBinning を使わず、シンプルなTH2コンストラクタに戻す
  TUnfoldDensity unfold(hMat, TUnfold::kHistMapOutputHoriz, regMode, constraintMode, densityFlags);

  if (unfold.SetInput(hDataTest) >= 10000) {
    std::cerr << "[WARNING] Unfolding input status error." << std::endl;
  }

  // ===========================
  // 5. ScanTau (Minimizing Global Correlation)
  // ===========================
  Int_t nScan = 30;
  TSpline* scanResult = 0; // log(tau) vs \bar{\rho}
  TSpline *logTauX = 0, *logTauY = 0;
  TGraph* lCurve = 0;

  std::cout << "[INFO] Scanning tau (Minimizing Global Correlation)..." << std::endl;

  // ScanTauを使用し、RhoAvgを最小化しつつ、Lカーブのデータも同時に取得する
  Int_t iBest = unfold.ScanTau(nScan, 1e-5, 1.0, &scanResult, TUnfoldDensity::kEScanTauRhoAvg, 0, 0, &lCurve, &logTauX, &logTauY);

  std::cout << "[RESULT] Best tau = " << unfold.GetTau() << std::endl;
  std::cout << "[RESULT] chi**2_A (Data) = " << unfold.GetChi2A() << std::endl;
  std::cout << "[RESULT] chi**2_L (Reg)  = " << unfold.GetChi2L() << std::endl;
  std::cout << "[RESULT] NDF = " << unfold.GetNdf() << std::endl;

  // ---------------------------
  // 最適点（iBest）のマーカー設定
  // ---------------------------
  Double_t tL[1], xL[1], yL[1];
  logTauX->GetKnot(iBest, tL[0], xL[0]);
  logTauY->GetKnot(iBest, tL[0], yL[0]);
  TGraph* bestLcurve = new TGraph(1, xL, yL); // Lカーブ用マーカー
  bestLcurve->SetMarkerColor(kRed);
  bestLcurve->SetMarkerStyle(29);
  bestLcurve->SetMarkerSize(2.0);

  Double_t tRho[1], rhoAvg[1];
  scanResult->GetKnot(iBest, tRho[0], rhoAvg[0]);
  TGraph* bestRho = new TGraph(1, tRho, rhoAvg); // Rho用マーカー
  bestRho->SetMarkerColor(kRed);
  bestRho->SetMarkerStyle(29);
  bestRho->SetMarkerSize(2.0);

  // ===========================
  // 6. Retrieve Results
  // ===========================
  TH1* hUnfolded = unfold.GetOutput("hUnfolded");

  // ===========================
  // 7. Plotting (2列 x 3行 に拡張)
  // ===========================
  // キャンバスを縦に少し長くする (1200 x 1200)
  TCanvas* c1 = new TCanvas("c1", "TUnfold Results (Square Matrix)", 1200, 1200);
  c1->Divide(2, 3);

  // Pad 1: L-Curve
  c1->cd(1);
  if (lCurve) {
    lCurve->SetTitle("L-Curve;log(#chi^{2}_{L});log(#chi^{2}_{A})");
    lCurve->Draw("AL");
    bestLcurve->Draw("* SAME");
  }

  // Pad 2: tau vs chi2 (L-Curve X)
  c1->cd(2);
  if (logTauX) {
    logTauX->SetTitle("L-Curve X-coordinate vs log(#tau);log(#tau);log(#chi^{2}_{L})");
    logTauX->Draw();
  }

  // Pad 3: log(tau) vs \bar{\rho} (新規追加)
  c1->cd(3);
  if (scanResult) {
    scanResult->SetTitle("Average Global Correlation;log(#tau);#bar{#rho} (Avg Global Correlation)");
    scanResult->SetLineColor(kBlue + 1);
    scanResult->SetLineWidth(2);
    scanResult->Draw();
    bestRho->Draw("* SAME"); // 最小値に赤星を打つ
  }

  // Pad 4: Spectra Comparison
  c1->cd(4);
  gPad->SetLogy();
  hGenTest->SetLineColor(kBlack);
  hGenTest->SetLineWidth(2);
  hGenTest->SetTitle("Spectra Comparison;p_{T}^{2} (GeV^{2}/c^{2});Counts");
  hGenTest->Draw("HIST");

  hDataTest->SetLineColor(kBlue);
  hDataTest->SetMarkerStyle(24);
  hDataTest->SetMarkerColor(kBlue);
  hDataTest->Draw("PE SAME");

  hUnfolded->SetLineColor(kRed);
  hUnfolded->SetMarkerStyle(20);
  hUnfolded->SetMarkerColor(kRed);
  hUnfolded->Draw("PE SAME");

  TLegend* leg = new TLegend(0.5, 0.7, 0.88, 0.88);
  leg->AddEntry(hGenTest, "Truth (Gen)", "l");
  leg->AddEntry(hDataTest, "Measured (Reco)", "pe");
  leg->AddEntry(hUnfolded, "Unfolded", "pe");
  leg->Draw();

  // Pad 5: Ratio (元のPad4)
  c1->cd(5);
  TH1D* hRatio = (TH1D*)hUnfolded->Clone("hRatio");
  hRatio->Divide(hGenTest);
  hRatio->SetTitle("Ratio (Unfolded / Truth);p_{T}^{2} (GeV^{2}/c^{2});Unfolded / Truth");
  hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hRatio->Draw("PE");
  TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw("SAME");

  // Pad 6 は空きスペース（必要に応じてテキストや設定情報を書き込めます）

  c1->SaveAs("TUnfold_AllGraphs_SquareMatrix.png");
}
