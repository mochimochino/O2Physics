#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>

// ==========================================
// カスタム関数の定義（マクロの外側に記述）
// ==========================================
Double_t BackgroundFitFunction(Double_t* x, Double_t* par)
{
  // 2.9 から 3.3 GeV/c^2 の J/psi シグナル領域をフィット評価から除外
  // 3.55 から 3.8 GeV/c^2 の Psi(2S) シグナル領域も除外
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  // 指数関数 exp(p0 + p1*x)
  return TMath::Exp(par[0] + par[1] * x[0]);
}

void Massfit()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0); // 統計ボックス(Stats/Fitボックス)を消す
  gStyle->SetTitleFontSize(0.04);

  TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/AnalysisResults.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Cannot open AnalysisResults.root" << std::endl;
    return;
  }

  TString dir = "my-upc-mass-01/";
  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dir + "MMuonUnlike"));
  TH1* hLike = dynamic_cast<TH1*>(f->Get(dir + "MMuonLike"));
  if (!hUnlike || !hLike)
    return;

  // 全体の描画範囲を定義
  double drawMin = 2.0;
  double drawMax = 5.0;

  // 1. Raw Histograms
  TCanvas* c1 = new TCanvas("c1", "Raw Histograms", 800, 600);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.12);
  gPad->SetLogy();
  gPad->SetGrid();

  int colUnlike = TColor::GetColor("#2e8b57");
  int colLike = TColor::GetColor("#d2691e");

  hUnlike->SetTitle("Dimuon Invariant Mass ALL; M_{#mu#mu} [GeV/c^{2}]; Counts");
  hUnlike->SetLineColor(colUnlike);
  hUnlike->SetMarkerColor(colUnlike);
  hUnlike->SetLineWidth(2);
  hUnlike->SetMarkerStyle(21);
  hUnlike->GetXaxis()->SetRangeUser(drawMin, drawMax);
  hUnlike->GetYaxis()->SetRangeUser(0.8, hUnlike->GetMaximum() * 5.0);
  hUnlike->Draw("E");

  hLike->SetLineColor(colLike);
  hLike->SetMarkerColor(colLike);
  hLike->SetLineWidth(2);
  hLike->SetMarkerStyle(25);
  hLike->Draw("E SAME");

  TLegend* leg = new TLegend(0.65, 0.80, 0.95, 0.95);
  leg->SetBorderSize(1);
  leg->AddEntry(hUnlike, "Unlike Sign (+-)", "lep");
  leg->AddEntry(hLike, "Like Sign (++ & --)", "lep");
  leg->Draw();
  c1->SaveAs("RawMass_ALL_Massfit.png");

  // ==================================================
  // 2. Signal Extraction + 2段階 Fitting
  // ==================================================
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->SetTitle("Signal (Unlike - Like) ALL; M_{#mu#mu} [GeV/c^{2}]; Counts");
  hSignal->Add(hLike, -1.0);

  TCanvas* c2 = new TCanvas("c2", "Signal Extraction and Mass Fit", 800, 800);

  // 上部パッド（メインのフィットグラフ用）
  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 下軸を消す設定
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);
  pad1->SetTopMargin(0.05);
  pad1->Draw();

  // 下部パッド（(Data-Fit)/Fit用）
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3); // X軸のラベル用
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);
  pad2->Draw();

  // 上のパッドに移動
  pad1->cd();
  gPad->SetGrid();

  double ymin = hSignal->GetMinimum();
  if (ymin > 0)
    ymin = 0;

  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);
  hSignal->GetYaxis()->SetRangeUser(ymin * 1.2, hSignal->GetMaximum() * 1.5);
  // 下部のパッドとX軸を共有するため、上部のX軸ラベル・タイトルは見えなくする
  hSignal->GetXaxis()->SetLabelSize(0);
  hSignal->GetXaxis()->SetTitleSize(0);
  hSignal->Draw("E");

  // --------------------------------------------------
  // Step A: サイドバンドフィット (アクセプタンス低下を避けるため 2.2 〜 7.0 に限定)
  // --------------------------------------------------
  double bkgFitMin = 2.2;
  double bkgFitMax = 7.0;

  TF1* fBkgSideband = new TF1("fBkgSideband", BackgroundFitFunction, bkgFitMin, bkgFitMax, 2);
  fBkgSideband->SetParameters(15.0, -2.0); // BKGの初期値を大きめに設定
  fBkgSideband->SetParNames("Expo_p0", "Expo_p1");
  hSignal->Fit("fBkgSideband", "R0"); // 描画しない(0)

  // --------------------------------------------------
  // Step B: 全体フィット (Crystal Ball(J/psi) + expo + Crystal Ball(Psi(2s)))
  // --------------------------------------------------
  double fitRangeMin = 2.2;
  double fitRangeMax = 5.0;

  // [0]-[4]: J/psi CB, [5]-[6]: Expo Bkg
  // [7]: Psi(2S) Amplitude, [8]: Psi(2S) Mean
  // Psi(2S) の Sigma は J/psi の Sigma * (M_Psi2S / M_Jpsi) と仮定
  // Alpha([3]) と n([4]) は J/psi と共通にする
  TF1* fitFunc = new TF1("fitFunc",
                         "crystalball(0) + expo(5) + [7]*ROOT::Math::crystalball_function(x, [3], [4], [2]*(3.686/3.096), [8])",
                         fitRangeMin, fitRangeMax);

  hSignal->GetXaxis()->SetRangeUser(2.8, 3.4);
  double maxValJpsi = hSignal->GetBinContent(hSignal->GetMaximumBin());
  hSignal->GetXaxis()->SetRangeUser(3.55, 3.8);
  double maxValPsi2s = hSignal->GetBinContent(hSignal->GetMaximumBin());
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax); // 軸の範囲を戻す

  if (maxValJpsi <= 0)
    maxValJpsi = 100;
  if (maxValPsi2s <= 0)
    maxValPsi2s = 10;

  // J/psi パラメータ設定
  fitFunc->SetParameters(maxValJpsi, 3.096, 0.05, 1.0, 5.0, fBkgSideband->GetParameter(0), fBkgSideband->GetParameter(1));
  fitFunc->SetParNames("Jpsi_Const", "Jpsi_Mean", "Jpsi_Sigma", "CB_Alpha", "CB_n", "Expo_p0", "Expo_p1");

  // Psi(2S) パラメータ設定
  fitFunc->SetParameter(7, maxValPsi2s);
  fitFunc->SetParameter(8, 3.686); // Psi(2S) mass
  fitFunc->SetParName(7, "Psi2S_Const");
  fitFunc->SetParName(8, "Psi2S_Mean");

  // J/psi パラメータ制限
  fitFunc->SetParLimits(1, 2.9, 3.25);
  fitFunc->SetParLimits(2, 0.02, 0.15);
  fitFunc->SetParLimits(3, 0.1, 5.0);
  fitFunc->SetParLimits(4, 0.1, 20.0);

  // Psi(2S) パラメータ制限 (収量が負にならないように Amplitude を正に制限)
  fitFunc->SetParLimits(7, 0.0, maxValJpsi);
  fitFunc->SetParLimits(8, 3.55, 3.8);

  // バックグラウンドの傾きを厳しく制限
  double p1_err = fBkgSideband->GetParError(1);
  fitFunc->SetParLimits(6, fBkgSideband->GetParameter(1) - 3 * p1_err, fBkgSideband->GetParameter(1) + 3 * p1_err);

  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineWidth(2);

  // 描画オプション"0"をつけてフィットし、後で全範囲に描画する
  hSignal->Fit("fitFunc", "RM0");

  // ==================================================
  // 3. 全範囲（1.0 - 10.0）への外挿と描画
  // ==================================================
  // ROOTの仕様上、Fit対象の関数はFit範囲でしか描画されないことがあるため、
  // 全範囲描画用に新しく関数を定義します。
  TF1* fitFuncDraw = new TF1("fitFuncDraw",
                             "crystalball(0) + expo(5) + [7]*ROOT::Math::crystalball_function(x, [3], [4], [2]*(3.686/3.096), [8])",
                             drawMin, drawMax);
  for (int i = 0; i < 9; i++) {
    fitFuncDraw->SetParameter(i, fitFunc->GetParameter(i));
  }
  fitFuncDraw->SetLineColor(kBlue);
  fitFuncDraw->SetLineWidth(2);
  fitFuncDraw->Draw("SAME");

  // バックグラウンド成分（expo）を全範囲で描画
  TF1* fitBkgDraw = new TF1("fitBkgDraw", "expo", drawMin, drawMax);
  fitBkgDraw->SetParameters(fitFunc->GetParameter(5), fitFunc->GetParameter(6));
  fitBkgDraw->SetLineColor(kBlue);
  fitBkgDraw->SetLineStyle(2);
  fitBkgDraw->Draw("SAME");

  // J/psi シグナル成分の描画と収量計算
  TF1* fitSigJpsi = new TF1("fitSigJpsi", "crystalball", drawMin, drawMax);
  for (int i = 0; i < 5; i++) {
    fitSigJpsi->SetParameter(i, fitFunc->GetParameter(i));
  }
  fitSigJpsi->SetLineColor(kMagenta);
  fitSigJpsi->SetLineStyle(2);
  fitSigJpsi->SetLineWidth(2);
  fitSigJpsi->Draw("SAME");

  // Psi(2S) シグナル成分の描画(赤色)と収量計算
  TF1* fitSigPsi2s = new TF1("fitSigPsi2s",
                             "[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4])", drawMin, drawMax);
  fitSigPsi2s->SetParameter(0, fitFunc->GetParameter(7));                   // Const
  fitSigPsi2s->SetParameter(1, fitFunc->GetParameter(3));                   // Alpha
  fitSigPsi2s->SetParameter(2, fitFunc->GetParameter(4));                   // n
  fitSigPsi2s->SetParameter(3, fitFunc->GetParameter(2) * (3.686 / 3.096)); // Sigma (Jpsiに質量比をかけたもの)
  fitSigPsi2s->SetParameter(4, fitFunc->GetParameter(8));                   // Mean
  fitSigPsi2s->SetLineColor(kRed);
  fitSigPsi2s->SetLineStyle(2);
  fitSigPsi2s->SetLineWidth(2);
  fitSigPsi2s->Draw("SAME");

  TLegend* leg2 = new TLegend(0.65, 0.65, 0.95, 0.95);
  leg2->SetBorderSize(1);
  leg2->AddEntry(hSignal, "Data (Unlike-Like)", "lep");
  leg2->AddEntry(fitFuncDraw, "Total Fit", "l");
  leg2->AddEntry(fitBkgDraw, "Background Fit", "l");
  leg2->AddEntry(fitSigJpsi, "J/#psi Signal", "l");
  leg2->AddEntry(fitSigPsi2s, "#psi(2S) Signal", "l");
  leg2->Draw();

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.15, 0.85, Form("M_{J/#psi} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
  latex->DrawLatex(0.15, 0.80, Form("#sigma_{J/#psi} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(2), fitFunc->GetParError(2)));
  latex->DrawLatex(0.15, 0.75, Form("M_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(8), fitFunc->GetParError(8)));

  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);
  double nJpsi = fitSigJpsi->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double nPsi2s = fitSigPsi2s->Integral(fitRangeMin, fitRangeMax) / binWidth;
  latex->DrawLatex(0.15, 0.70, Form("N_{J/#psi} #approx %.0f, N_{#psi(2S)} #approx %.0f", nJpsi, nPsi2s));

  // ==================================================
  // 4. (Data - Fit) / Fit の計算と下部パッドへの描画
  // ==================================================
  pad2->cd(); // 下のパッドに移動
  gPad->SetGridy();

  TH1D* hRatio = (TH1D*)hSignal->Clone("hRatio");
  hRatio->SetTitle(""); // タイトルは消す

  // 各ビンについて (Data - Fit) / Fit を計算
  for (int i = 1; i <= hRatio->GetNbinsX(); i++) {
    double x = hRatio->GetBinCenter(i);
    // 描画範囲内のみ計算
    if (x >= drawMin && x <= drawMax) {
      double data = hSignal->GetBinContent(i);
      double dataErr = hSignal->GetBinError(i);
      double fitVal = fitFuncDraw->Eval(x);

      if (fitVal > 0) {
        hRatio->SetBinContent(i, (data - fitVal) / fitVal);
        hRatio->SetBinError(i, dataErr / fitVal); // 誤差の伝播（Fit関数の誤差は無視と仮定）
      } else {
        hRatio->SetBinContent(i, 0);
        hRatio->SetBinError(i, 0);
      }
    } else {
      hRatio->SetBinContent(i, 0);
      hRatio->SetBinError(i, 0);
    }
  }

  // 下部パッドの見栄えの調整（比率に合わせて文字幅などを大きくする）
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->GetYaxis()->SetTitle("(Data-Fit)/Fit");
  hRatio->GetYaxis()->SetRangeUser(-1.5, 1.5); // Y軸の表示範囲は適宜調整
  hRatio->GetYaxis()->SetNdivisions(505);

  // Y軸ラベル
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetTitleOffset(0.4);

  // X軸ラベル
  hRatio->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/c^{2}]");
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetTitleOffset(0.9);

  hRatio->Draw("EP");

  // Y=0 の基準線を描画
  TF1* lineZero = new TF1("lineZero", "0", drawMin, drawMax);
  lineZero->SetLineColor(kRed);
  lineZero->SetLineStyle(2);
  lineZero->Draw("SAME");

  c2->SaveAs("Unlike_bkg_WideRange_Massfit.png");
  f->Close();
}
