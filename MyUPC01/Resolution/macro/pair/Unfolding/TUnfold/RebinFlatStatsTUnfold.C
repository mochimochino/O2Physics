#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------
// ヘルパー関数: 指定したヒストグラムから等統計(Flat Stats)になるようなビン境界を計算する
// ---------------------------------------------------------
std::vector<double> CalculateFlatStatsBins(TH1D* h, int nTargetBins, double xMin, double xMax)
{
  int binMin = h->FindBin(xMin);
  int binMax = h->FindBin(xMax);

  std::vector<double> cumSum;
  cumSum.push_back(0.0);
  for (int i = binMin; i <= binMax; ++i) {
    cumSum.push_back(cumSum.back() + h->GetBinContent(i));
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

    double edge = h->GetBinLowEdge(binMin + bestIdx);
    if (edge > bins.back() && edge < xMax) {
      bins.push_back(edge);
    }
  }

  if (bins.back() < xMax) {
    bins.push_back(xMax);
  }

  return bins;
}

void RebinFlatStatsTUnfold(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/flatstats/3bin/02/Step2_Rebinned.root",
  int nGenTargetBins = 3,
  double xMin = 0.0,
  double xMax = 0.020)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  int nRecoTargetBins = nGenTargetBins * 2;

  std::cout << "========================================" << std::endl;
  std::cout << "[Info] Target Gen Bins: " << nGenTargetBins << std::endl;
  std::cout << "[Info] Target Reco Bins: " << nRecoTargetBins << " (Gen * 2 for TUnfold)" << std::endl;
  std::cout << "========================================" << std::endl;

  // --- 1. 元データの読み込み ---
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open input file: " << inputFile << std::endl;
    return;
  }

  TH1D* hGenFine = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hRecoFine = (TH1D*)fIn->Get("hRecoPt2");
  TH2D* hMatrixFine = (TH2D*)fIn->Get("hResponseMatrixPairPt2"); // 元の細かい2Dヒストグラム

  if (!hGenFine || !hRecoFine || !hMatrixFine) {
    std::cerr << "[Error] Required histograms (hGenPt2, hRecoPt2, hResponseMatrixPairPt2) not found in input file!" << std::endl;
    fIn->Close();
    return;
  }

  // --- 2. 独立したFlat Statsビン境界の計算 ---
  std::vector<double> binsGen = CalculateFlatStatsBins(hGenFine, nGenTargetBins, xMin, xMax);
  std::vector<double> binsReco = CalculateFlatStatsBins(hRecoFine, nRecoTargetBins, xMin, xMax);

  int nBinsGen = binsGen.size() - 1;
  int nBinsReco = binsReco.size() - 1;

  // 実際のビン数がアルゴリズムの都合でターゲットと異なる場合があるため警告を出す
  if (nBinsGen != nGenTargetBins || nBinsReco != nRecoTargetBins) {
    std::cout << "[Warning] Actual bins slightly differ from target due to edge calculation." << std::endl;
  }

  std::cout << "[Info] Calculated Gen Bins: " << nBinsGen << " | Calculated Reco Bins: " << nBinsReco << std::endl;

  // --- 3. 1DヒストグラムのRebin ---
  TH1D* hGenRebin = (TH1D*)hGenFine->Rebin(nBinsGen, "hGenPt2_rebin", binsGen.data());
  TH1D* hRecoRebin = (TH1D*)hRecoFine->Rebin(nBinsReco, "hRecoPt2_rebin", binsReco.data());

  // --- 4. 非正方(Rectangular)な応答行列(TH2D)の作成 ---
  // X軸: Gen (真の分布)
  // Y軸: Reco (観測分布)
  TH2D* hMatrixRebin = new TH2D("hResponseMatrix", "Response Matrix (Gen:X, Reco:Y)",
                                nBinsGen, binsGen.data(),
                                nBinsReco, binsReco.data());

  // 元の細かいTH2Dの各ビンのコンテンツを、新しい粗いビンのTH2Dへ振り分ける
  for (int i = 1; i <= hMatrixFine->GetNbinsX(); ++i) {
    for (int j = 1; j <= hMatrixFine->GetNbinsY(); ++j) {
      double content = hMatrixFine->GetBinContent(i, j);
      if (content > 0) {
        double xCenter = hMatrixFine->GetXaxis()->GetBinCenter(i); // Gen
        double yCenter = hMatrixFine->GetYaxis()->GetBinCenter(j); // Reco
        hMatrixRebin->Fill(xCenter, yCenter, content);
      }
    }
  }

  // --- 5. 結果の保存 ---
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  hGenRebin->Write();
  hRecoRebin->Write();
  hMatrixRebin->Write();
  fOut->Close();

  std::cout << "[Info] Saved ROOT file: " << outputFile << std::endl;

  // --- 6. 描画とPNG保存 (オプション) ---
  TString outBase = outputFile;
  outBase.ReplaceAll(".root", "");

  // 1D描画ヘルパー
  auto drawAndSave1D = [&](TH1D* h, const TString& title, const TString& suffix) {
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
    tex.DrawLatex(0.55, 0.79, "Method: Flat Stats");
    tex.DrawLatex(0.55, 0.74, Form("Bins: %d", h->GetNbinsX()));

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  // 2D行列描画
  auto drawAndSave2D = [&](TH2D* h2) {
    TCanvas* c2 = new TCanvas("c_Matrix", "Response Matrix", 800, 800);
    c2->SetRightMargin(0.15); // カラーパレット用
    h2->SetStats(0);
    h2->SetTitle("Response Matrix;Gen p_{T}^{2} (GeV^{2}/c^{2});Reco p_{T}^{2} (GeV^{2}/c^{2})");
    gPad->SetLogz(); // 行列は対数スケールが見やすい
    h2->Draw("COLZ");

    c2->SaveAs(outBase + "_Matrix.png");
    delete c2;
  };

  drawAndSave1D(hGenRebin, "Gen p_{T}^{2} (Flat Stats)", "Gen");
  drawAndSave1D(hRecoRebin, "Reco p_{T}^{2} (Flat Stats)", "Reco");
  drawAndSave2D(hMatrixRebin);

  std::cout << "[Info] Saved PNGs to: " << outBase << "_*.png" << std::endl;

  fIn->Close();
}
