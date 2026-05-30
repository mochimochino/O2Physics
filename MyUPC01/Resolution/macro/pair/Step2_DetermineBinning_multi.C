// ============================================================
//  Step2_DetermineBinning_multi.C
//  複数のBin幅決定手法を選択して使用できる統合マクロ。
//  出力ヒストはすべて RooUnfoldBayes.C が読むキー名に統一。
//
//  出力 ROOT ファイル内容
//    hResponseMatrix   : リビン済み 2D レスポンス行列
//    hMCGen_Rebinned   : リビン済み Gen 1D ヒスト
//    hMCReco_Rebinned  : リビン済み Reco 1D ヒスト
//    hDataReco_Rebinned: リビン済み Data Reco 1D ヒスト（任意）
//    hPurity           : Purity [%] per bin
//    hStability        : Stability [%] per bin
// ============================================================

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ============================================================
// Bin幅決定手法の選択
// ============================================================
enum BinMethod {
  kTopDown      = 0, // Purity/Stability 下限を満たす最大幅を上から確定
  kEqualWidth   = 1, // 等幅
  kFlatStats    = 2, // Gen 統計量が各ビンで均一（等統計）
  kConstantDiff = 3, // 等差拡大（初期幅 + k * Δw）
  kEqualEvents  = 4, // Reco 統計量が各ビンで均一
};

// ============================================================
// 各手法の実装
// ============================================================

// --- TopDown ---
std::vector<double> TopDownEdges(TH2F* hMatrix, double targetPurity, double targetStability,
                                 double minStats, double xMax)
{
  std::vector<double> edges;
  edges.push_back(xMax);

  int maxBin = hMatrix->GetXaxis()->FindBin(xMax - 1e-6);
  int upper = maxBin;

  for (int lower = maxBin; lower >= 1; --lower) {
    if (lower % 100 == 0 || lower == 1) {
      int pct = 100 - (int)(100.0 * lower / maxBin);
      std::cout << "\r[TopDown] Scanning... " << pct << "%" << std::flush;
    }
    double diag    = hMatrix->Integral(lower, upper, lower, upper);
    double sumTrue = hMatrix->Integral(lower, upper, 1, hMatrix->GetNbinsY());
    double sumReco = hMatrix->Integral(1, hMatrix->GetNbinsX(), lower, upper);
    if (sumTrue <= 0 || sumReco <= 0) continue;
    if (diag / sumReco >= targetPurity &&
        diag / sumTrue >= targetStability &&
        sumReco >= minStats && sumTrue >= minStats) {
      double lo = hMatrix->GetXaxis()->GetBinLowEdge(lower);
      if (lower > 1 && lo > 0.0) {
        edges.push_back(lo);
        upper = lower - 1;
      }
    }
  }
  std::cout << "\r[TopDown] Scanning... 100%\n";

  if (edges.back() != 0.0) edges.push_back(0.0);
  std::reverse(edges.begin(), edges.end());
  return edges;
}

// --- EqualWidth ---
std::vector<double> EqualWidthEdges(double binWidth, double xMin, double xMax)
{
  std::vector<double> edges;
  int nBins = (int)((xMax - xMin) / binWidth);
  for (int i = 0; i <= nBins; ++i) edges.push_back(xMin + i * binWidth);
  if (edges.back() < xMax) edges.push_back(xMax);
  return edges;
}

// --- FlatStats (Gen) ---
std::vector<double> FlatStatsEdges(TH1D* hGen, int nTargetBins, double xMin, double xMax)
{
  std::vector<double> edges;
  edges.push_back(xMin);

  int binMin = hGen->FindBin(xMin + 1e-9);
  int binMax = hGen->FindBin(xMax - 1e-9);
  double total = hGen->Integral(binMin, binMax);
  double step  = total / nTargetBins;

  double cumSum   = 0;
  double targetSum = step;
  int binsFound    = 0;

  for (int i = binMin; i <= binMax; ++i) {
    cumSum += hGen->GetBinContent(i);
    if (cumSum >= targetSum && binsFound < nTargetBins - 1) {
      edges.push_back(hGen->GetXaxis()->GetBinUpEdge(i));
      targetSum += step;
      binsFound++;
    }
  }
  edges.push_back(xMax);
  return edges;
}

// --- ConstantDiff ---
std::vector<double> ConstantDiffEdges(double w0, double dw, double xMin, double xMax)
{
  std::vector<double> edges;
  edges.push_back(xMin);
  double x = xMin;
  double w = w0;
  while (x < xMax) {
    x += w;
    if (x >= xMax) { edges.push_back(xMax); break; }
    edges.push_back(x);
    w += dw;
  }
  return edges;
}

// --- EqualEvents (Reco) ---
std::vector<double> EqualEventsEdges(TH1D* hReco, int nTargetBins, double xMin, double xMax)
{
  std::vector<double> edges;
  edges.push_back(xMin);

  int binMax = hReco->GetXaxis()->FindBin(xMax - 1e-9);
  double total = hReco->Integral(1, binMax);
  double evPerBin = total / nTargetBins;

  double cumSum = 0;
  int found = 0;
  for (int i = 1; i <= binMax; ++i) {
    cumSum += hReco->GetBinContent(i);
    if (cumSum >= evPerBin && found < nTargetBins - 1) {
      edges.push_back(hReco->GetXaxis()->GetBinUpEdge(i));
      cumSum = 0;
      found++;
    }
  }
  if (edges.back() < xMax) edges.push_back(xMax);
  else edges.back() = xMax;
  return edges;
}

// ============================================================
// Main
// ============================================================
void Step2_DetermineBinning_multi()
{
  // ----------------------------------------------------------
  // ===  ユーザー設定  ===
  // ----------------------------------------------------------

  // --- 入出力 ---
  const TString inFile       = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Incoherent/Step1_merged.root";
  const TString dataFile     = inFile; // Dataが別ファイルなら変更。同じでもOK
  const TString outDir       = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Incoherent/";
  const TString outFile      = outDir + "Step2_Response_for_Unfolding.root";

  // --- Bin幅決定手法の選択 ---
  const BinMethod method = kFlatStats;
  // kTopDown      : Purity/Stability 下限を満たす最大幅を上から確定
  // kEqualWidth   : 等幅ビン
  // kFlatStats    : Gen 統計量が各ビンで均一（等統計）
  // kConstantDiff : 等差拡大（低 pT2 が細かく、高 pT2 が粗い）
  // kEqualEvents  : Reco 統計量が各ビンで均一

  // --- 解析範囲 ---
  const double xMin  = 0.0;
  const double xMax  = 1.2;

  // --- TopDown 用パラメータ ---
  const double targetPurity    = 0.45;
  const double targetStability = 0.45;
  const double minStats        = 500.0;

  // --- EqualWidth 用パラメータ ---
  const double binWidth = 0.08;

  // --- FlatStats / EqualEvents 用パラメータ ---
  const int nTargetBins = 8;

  // --- ConstantDiff 用パラメータ ---
  const double w0 = 0.016;
  const double dw = 0.012;

  // ----------------------------------------------------------

  gStyle->SetOptStat(0);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open: " << inFile << std::endl;
    return;
  }
  TH2F* hMatrixFine = (TH2F*)fIn->Get("hResponseMatrixPairPt2");
  TH1D* hGenFine    = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hRecoFine   = (TH1D*)fIn->Get("hRecoPt2");
  if (!hMatrixFine || !hGenFine || !hRecoFine) {
    std::cerr << "[Error] Required histograms not found." << std::endl;
    return;
  }

  TFile* fData = TFile::Open(dataFile, "READ");
  TH1D* hDataFine = nullptr;
  if (fData && !fData->IsZombie())
    hDataFine = (TH1D*)fData->Get("hRecoPt2");

  // ----------------------------------------------------------
  // ビン境界の決定
  // ----------------------------------------------------------
  std::vector<double> edges;
  TString methodName;

  switch (method) {
    case kTopDown:
      methodName = "TopDown";
      edges = TopDownEdges(hMatrixFine, targetPurity, targetStability, minStats, xMax);
      break;
    case kEqualWidth:
      methodName = "EqualWidth";
      edges = EqualWidthEdges(binWidth, xMin, xMax);
      break;
    case kFlatStats:
      methodName = "FlatStats";
      edges = FlatStatsEdges(hGenFine, nTargetBins, xMin, xMax);
      break;
    case kConstantDiff:
      methodName = "ConstantDiff";
      edges = ConstantDiffEdges(w0, dw, xMin, xMax);
      break;
    case kEqualEvents:
      methodName = "EqualEvents";
      edges = EqualEventsEdges(hRecoFine, nTargetBins, xMin, xMax);
      break;
  }

  const int nBins = (int)edges.size() - 1;
  if (nBins <= 0) {
    std::cerr << "[Error] No bins determined. Adjust parameters." << std::endl;
    return;
  }
  std::cout << "[Info] Method: " << methodName << " -> " << nBins << " bins" << std::endl;

  // ----------------------------------------------------------
  // リビン: 2D レスポンス行列
  // ----------------------------------------------------------
  TH2D* hMatrix = new TH2D("hResponseMatrix",
                            Form("Response Matrix (%s);Gen p_{T}^{2} (GeV^{2}/c^{2});Reco p_{T}^{2} (GeV^{2}/c^{2})", methodName.Data()),
                            nBins, edges.data(), nBins, edges.data());
  for (int ix = 1; ix <= hMatrixFine->GetNbinsX(); ++ix) {
    double gv = hMatrixFine->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= hMatrixFine->GetNbinsY(); ++iy) {
      hMatrix->Fill(gv, hMatrixFine->GetYaxis()->GetBinCenter(iy),
                    hMatrixFine->GetBinContent(ix, iy));
    }
  }

  // ----------------------------------------------------------
  // リビン: 1D ヒスト
  // ----------------------------------------------------------
  TH1D* hGen    = hGenFine  ? (TH1D*)hGenFine ->Rebin(nBins, "hMCGen_Rebinned",  edges.data()) : nullptr;
  TH1D* hReco   = hRecoFine ? (TH1D*)hRecoFine->Rebin(nBins, "hMCReco_Rebinned", edges.data()) : nullptr;
  TH1D* hData   = hDataFine ? (TH1D*)hDataFine->Rebin(nBins, "hDataReco_Rebinned", edges.data()) : nullptr;

  // ----------------------------------------------------------
  // Purity / Stability 計算
  // ----------------------------------------------------------
  TH1D* hPurity    = new TH1D("hPurity",    "Purity [%]",    nBins, edges.data());
  TH1D* hStability = new TH1D("hStability", "Stability [%]", nBins, edges.data());

  for (int i = 1; i <= nBins; ++i) {
    double diag    = hMatrix->GetBinContent(i, i);
    double sumReco = hMatrix->Integral(1, nBins, i, i);
    double sumTrue = hMatrix->Integral(i, i, 1, nBins);
    hPurity   ->SetBinContent(i, sumReco > 0 ? diag / sumReco * 100.0 : 0.0);
    hStability->SetBinContent(i, sumTrue > 0 ? diag / sumTrue * 100.0 : 0.0);
  }

  // ----------------------------------------------------------
  // 検証テーブルのターミナル出力
  // ----------------------------------------------------------
  std::cout << "\n============================================================\n";
  std::cout << Form(" Binning Summary — %s (%d bins)\n", methodName.Data(), nBins);
  std::cout << "============================================================\n";
  std::cout << Form("%-5s | %-15s | %-8s | %-10s | %-8s | %-8s",
                    "Bin", "pT2 Range", "Events", "StatErr%", "Purity%", "Stab%") << "\n";
  std::cout << "------------------------------------------------------------\n";

  bool hasError = false;
  for (int i = 1; i <= nBins; ++i) {
    double ev = hData ? hData->GetBinContent(i) :
               (hReco ? hReco->GetBinContent(i) : 0);
    double se  = ev > 0 ? 100.0 / std::sqrt(ev) : 999.0;
    double pur = hPurity->GetBinContent(i);
    double sta = hStability->GetBinContent(i);
    TString flag = (se > 20.0) ? " [!]" : "";

    std::cout << Form("%-5d | %6.4f-%6.4f | %-8.0f | %-10.2f | %-8.1f | %-8.1f%s",
                      i, edges[i-1], edges[i], ev, se, pur, sta, flag.Data()) << "\n";
    if (se > 20.0) hasError = true;
  }
  std::cout << "============================================================\n";
  if (hasError)
    std::cerr << "[WARNING] Some bins have stat error > 20%. Consider fewer bins.\n";
  else
    std::cout << "[OK] All bins pass validation.\n";
  std::cout << "\n";

  // ----------------------------------------------------------
  // 保存
  // ----------------------------------------------------------
  TFile* fOut = TFile::Open(outFile, "RECREATE");
  hMatrix   ->Write();
  hPurity   ->Write();
  hStability->Write();
  if (hGen)  hGen ->Write();
  if (hReco) hReco->Write();
  if (hData) hData->Write();
  fOut->Close();

  std::cout << "[Done] " << outFile << "\n";
}
