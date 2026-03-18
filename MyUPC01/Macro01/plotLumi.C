#include "TAxis.h" // 追加
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>

void plotLumi()
{
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open("AnalysisResults.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open AnalysisResults.root" << std::endl;
    return;
  }

  TH1D* hEvents = (TH1D*)f->Get("my-upc-luminosity/hEvents_All");
  if (!hEvents) {
    std::cerr << "Cannot find hEvents_All in AnalysisResults.root" << std::endl;
    return;
  }

  // 2023年 Pb-Pb衝突用の可視断面積（プレースホルダー）
  // ※必ずALICE公式の内部ノートから該当トリガーの可視断面積を取得して書き換えてください。
  double sigma_visible_mb = 7700.0; // 例: Pb-Pbの全断面積に近い値を仮置き (mb)

  // mb から μb への変換 (1 mb = 1000 μb)
  double sigma_visible_inver_mub = sigma_visible_mb * 1e3;

  TH1D* hLumi = (TH1D*)hEvents->Clone("hLumi");

  // 修正1: タイトルのLaTeX表記を変更 (#mathcal -> #it, #mu b -> #mub)
  hLumi->SetTitle("Integrated Luminosity per Run (Skimmed Data); Run Number; #it{L}_{int} [#mub^{-1}]");
  hLumi->Reset();

  for (int i = 1; i <= hEvents->GetNbinsX(); i++) {
    double nEvents = hEvents->GetBinContent(i);
    if (nEvents > 0) {
      // ルミノシティの計算 (イベント数 / 断面積[μb])
      double lumi = nEvents / sigma_visible_inver_mub;
      double lumiErr = TMath::Sqrt(nEvents) / sigma_visible_inver_mub;
      hLumi->SetBinContent(i, lumi);
      hLumi->SetBinError(i, lumiErr);
    }
  }

  TCanvas* c = new TCanvas("cLumi", "Luminosity per Run", 800, 600);
  c->SetGrid();
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);

  // 修正2: X軸の指数表記（10^x）をオフにして、Run番号をそのまま表示
  hLumi->GetXaxis()->SetNoExponent(kTRUE);

  // 線の設定
  hLumi->SetLineColor(kBlue);
  hLumi->SetLineWidth(2); // 線を少し太くして見やすくする

  int firstBin = hLumi->FindFirstBinAbove(0);
  int lastBin = hLumi->FindLastBinAbove(0);
  if (firstBin > 0 && lastBin >= firstBin) {
    hLumi->GetXaxis()->SetRange(firstBin - 2, lastBin + 2);
  }

  // 修正3: 点(P)と誤差(E)ではなく、ヒストグラムの線(HIST)として描画
  hLumi->Draw("HIST");

  c->SaveAs("LuminosityPerRun.png");

  std::cout << "Successfully plotted Luminosity to LuminosityPerRun.png" << std::endl;
}
