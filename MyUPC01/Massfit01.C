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

// --------------------------------------------------
// Dimuon Continuum: Expo3
//   f(m) = p0 * exp(p1*m + p2*m^2)
// par[0]=p0, par[1]=p1, par[2]=p2
// --------------------------------------------------
Double_t Expo3(Double_t* x, Double_t* par)
{
  double m = x[0];
  return par[0] * TMath::Exp(par[1] * m + par[2] * m * m);
}

// --------------------------------------------------
// Background: ExpoPol4
//   f(m) = p0 * exp(p1*m) * (p2 + p3*m + p4*m^2 + p5*m^3 + p6*m^4)
// par[0]=p0, par[1]=p1, par[2]=p2, par[3]=p3, par[4]=p4, par[5]=p5, par[6]=p6
// --------------------------------------------------
Double_t ExpoPol4(Double_t* x, Double_t* par)
{
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m + par[5] * m * m * m + par[6] * m * m * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}

// VMG (Variable width Gaussian) を使ったサイドバンドフィット用関数
// シグナル領域を除外してサイドバンドでExpoPol4をフィット
Double_t ExpoPol4Sideband(Double_t* x, Double_t* par)
{
  // J/psi と psi(2S) のシグナル領域を除外
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m + par[5] * m * m * m + par[6] * m * m * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}

// --------------------------------------------------
// Background: VMG (Variable Width Gaussian)
//   sigma = A + B * (m - mbar) / mbar
//   f(m) = N * exp( -(m - mbar)^2 / (2 * sigma^2) )
// par[0]=N, par[1]=mbar, par[2]=A (sigma_0), par[3]=B
// --------------------------------------------------
Double_t VMG(Double_t* x, Double_t* par)
{
  double m = x[0];
  double mbar = par[1];
  double sigma = par[2] + par[3] * (m - mbar) / mbar;
  if (sigma <= 0)
    return 0;
  return par[0] * TMath::Exp(-(m - mbar) * (m - mbar) / (2.0 * sigma * sigma));
}

// --------------------------------------------------
// J/psi and psi(2S): Double-Sided Crystal Ball (DSCB)
//
// alpha = (m - m0) / sigma
// A1 = (n1/|alpha1|)^n1 * exp(-|alpha1|^2/2),  B1 = n1/|alpha1| - |alpha1|
// A2 = (n2/|alpha2|)^n2 * exp(-|alpha2|^2/2),  D2 = n2/|alpha2| - |alpha2|
//
// f = N * { A1*(B1 - alpha)^(-n1)  if alpha <= -|alpha1|
//         { exp(-alpha^2/2)         if -|alpha1| < alpha < |alpha2|
//         { A2*(D2 + alpha)^(-n2)  if alpha >= |alpha2| }
//
// par[0]=N, par[1]=m0, par[2]=sigma, par[3]=alpha1, par[4]=n1, par[5]=alpha2, par[6]=n2
// --------------------------------------------------
Double_t DSCB(Double_t* x, Double_t* par)
{
  double N      = par[0];
  double m0     = par[1];
  double sigma  = par[2];
  double alpha1 = TMath::Abs(par[3]);
  double n1     = par[4];
  double alpha2 = TMath::Abs(par[5]);
  double n2     = par[6];

  double alpha_var = (x[0] - m0) / sigma;

  double A1 = TMath::Power(n1 / alpha1, n1) * TMath::Exp(-0.5 * alpha1 * alpha1);
  double B1 = n1 / alpha1 - alpha1;
  double A2 = TMath::Power(n2 / alpha2, n2) * TMath::Exp(-0.5 * alpha2 * alpha2);
  double D2 = n2 / alpha2 - alpha2;

  double val = 0;
  if (alpha_var <= -alpha1) {
    val = A1 * TMath::Power(B1 - alpha_var, -n1);
  } else if (alpha_var < alpha2) {
    val = TMath::Exp(-0.5 * alpha_var * alpha_var);
  } else {
    val = A2 * TMath::Power(D2 + alpha_var, -n2);
  }
  return N * val;
}

// --------------------------------------------------
// Total fit function:
//   DSCB(J/psi) + DSCB(psi(2S)) + ExpoPol4 (background)
//
// J/psi  DSCB: par[0..6]  (N, m0, sigma, alpha1, n1, alpha2, n2)
// psi2S  DSCB: par[7..13] (N, m0, sigma, alpha1, n1, alpha2, n2)
//     -> sigma psi2S は J/psi sigma に質量比をかけて固定してもよいが、ここでは独立とする
// ExpoPol4:    par[14..20] (p0, p1, p2, p3, p4, p5, p6)
// --------------------------------------------------
Double_t TotalFit(Double_t* x, Double_t* par)
{
  // J/psi DSCB
  double jpsi = DSCB(x, &par[0]);

  // psi(2S) DSCB
  double psi2s = DSCB(x, &par[7]);

  // ExpoPol4 background
  double bkg = ExpoPol4(x, &par[14]);

  return jpsi + psi2s + bkg;
}

// ==========================================
// メインマクロ
// ==========================================
void Massfit01()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFontSize(0.04);

  TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/AnalysisResults.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Cannot open AnalysisResults.root" << std::endl;
    return;
  }

  TString dir = "my-upc-mass-01/";
  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dir + "MMuonUnlike"));
  TH1* hLike   = dynamic_cast<TH1*>(f->Get(dir + "MMuonLike"));
  if (!hUnlike || !hLike) {
    std::cerr << "Error: Cannot get histograms" << std::endl;
    return;
  }

  // 全体の描画範囲
  double drawMin = 2.0;
  double drawMax = 5.0;

  // =====================================================
  // 1. Raw Histograms
  // =====================================================
  TCanvas* c1 = new TCanvas("c1", "Raw Histograms", 800, 600);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.12);
  gPad->SetLogy();
  gPad->SetGrid();

  int colUnlike = TColor::GetColor("#2e8b57");
  int colLike   = TColor::GetColor("#d2691e");

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
  leg->AddEntry(hLike,   "Like Sign (++ & --)", "lep");
  leg->Draw();
  c1->SaveAs("RawMass_ALL_Massfit01.png");

  // =====================================================
  // 2. Signal (Unlike - Like)
  // =====================================================
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->SetTitle("Signal (Unlike - Like) ALL; M_{#mu#mu} [GeV/c^{2}]; Counts");
  hSignal->Add(hLike, -1.0);

  TCanvas* c2 = new TCanvas("c2", "Signal Extraction and Mass Fit (Massfit01)", 800, 800);

  // 上部パッド
  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);
  pad1->SetTopMargin(0.05);
  pad1->Draw();

  // 下部パッド
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.30);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);
  pad2->Draw();

  pad1->cd();
  gPad->SetGrid();

  double ymin = hSignal->GetMinimum();
  if (ymin > 0) ymin = 0;

  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);
  hSignal->GetYaxis()->SetRangeUser(ymin * 1.2, hSignal->GetMaximum() * 1.5);
  hSignal->GetXaxis()->SetLabelSize(0);
  hSignal->GetXaxis()->SetTitleSize(0);
  hSignal->Draw("E");

  // =====================================================
  // Step A: ExpoPol4 サイドバンドフィット (シグナル領域除外)
  //   フィット範囲: 2.2 - 4.5 GeV
  // =====================================================
  double bkgFitMin = 2.2;
  double bkgFitMax = 4.5;

  TF1* fBkgSideband = new TF1("fBkgSideband", ExpoPol4Sideband, bkgFitMin, bkgFitMax, 7);
  // 初期値設定
  fBkgSideband->SetParameters(1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
  fBkgSideband->SetParNames("BkgP0", "BkgP1", "BkgP2", "BkgP3", "BkgP4", "BkgP5", "BkgP6");
  hSignal->Fit("fBkgSideband", "R0"); // 描画しない

  // =====================================================
  // Step B: 全体フィット
  //   J/psi DSCB [0..6] + psi(2S) DSCB [7..13] + ExpoPol4 bkg [14..20]
  //   フィット範囲: 2.5 - 4.5 GeV
  // =====================================================
  double fitRangeMin = 2.5;
  double fitRangeMax = 4.5;

  // J/psi付近のピーク値を概算
  hSignal->GetXaxis()->SetRangeUser(2.9, 3.3);
  double maxValJpsi = hSignal->GetBinContent(hSignal->GetMaximumBin());
  hSignal->GetXaxis()->SetRangeUser(3.55, 3.8);
  double maxValPsi2s = hSignal->GetBinContent(hSignal->GetMaximumBin());
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);

  if (maxValJpsi <= 0)  maxValJpsi  = 100.0;
  if (maxValPsi2s <= 0) maxValPsi2s = 10.0;

  TF1* fitFunc = new TF1("fitFunc", TotalFit, fitRangeMin, fitRangeMax, 21);

  // --- J/psi DSCB パラメータ (par[0..6]) ---
  // par[0]=N, par[1]=m0, par[2]=sigma, par[3]=alpha1, par[4]=n1, par[5]=alpha2, par[6]=n2
  fitFunc->SetParameter(0,  maxValJpsi); // N
  fitFunc->SetParameter(1,  3.096);      // m0 (J/psi mass)
  fitFunc->SetParameter(2,  0.05);       // sigma
  fitFunc->SetParameter(3,  1.5);        // alpha1
  fitFunc->SetParameter(4,  5.0);        // n1
  fitFunc->SetParameter(5,  1.5);        // alpha2
  fitFunc->SetParameter(6,  5.0);        // n2
  fitFunc->SetParName(0, "Jpsi_N");
  fitFunc->SetParName(1, "Jpsi_m0");
  fitFunc->SetParName(2, "Jpsi_sigma");
  fitFunc->SetParName(3, "Jpsi_alpha1");
  fitFunc->SetParName(4, "Jpsi_n1");
  fitFunc->SetParName(5, "Jpsi_alpha2");
  fitFunc->SetParName(6, "Jpsi_n2");

  // パラメータ制限
  fitFunc->SetParLimits(0, 0.0, maxValJpsi * 5.0);
  fitFunc->SetParLimits(1, 2.9, 3.25);
  fitFunc->SetParLimits(2, 0.02, 0.15);
  fitFunc->SetParLimits(3, 0.1, 10.0);
  fitFunc->SetParLimits(4, 0.1, 50.0);
  fitFunc->SetParLimits(5, 0.1, 10.0);
  fitFunc->SetParLimits(6, 0.1, 50.0);

  // --- psi(2S) DSCB パラメータ (par[7..13]) ---
  // par[7]=N, par[8]=m0, par[9]=sigma, par[10]=alpha1, par[11]=n1, par[12]=alpha2, par[13]=n2
  // sigma の初期値は J/psi の比で設定
  fitFunc->SetParameter(7,  maxValPsi2s);
  fitFunc->SetParameter(8,  3.686);            // psi(2S) mass
  fitFunc->SetParameter(9,  0.05 * (3.686 / 3.096)); // scaled sigma
  fitFunc->SetParameter(10, 1.5);
  fitFunc->SetParameter(11, 5.0);
  fitFunc->SetParameter(12, 1.5);
  fitFunc->SetParameter(13, 5.0);
  fitFunc->SetParName(7,  "Psi2S_N");
  fitFunc->SetParName(8,  "Psi2S_m0");
  fitFunc->SetParName(9,  "Psi2S_sigma");
  fitFunc->SetParName(10, "Psi2S_alpha1");
  fitFunc->SetParName(11, "Psi2S_n1");
  fitFunc->SetParName(12, "Psi2S_alpha2");
  fitFunc->SetParName(13, "Psi2S_n2");

  fitFunc->SetParLimits(7,  0.0, maxValJpsi);
  fitFunc->SetParLimits(8,  3.55, 3.80);
  fitFunc->SetParLimits(9,  0.02, 0.20);
  fitFunc->SetParLimits(10, 0.1, 10.0);
  fitFunc->SetParLimits(11, 0.1, 50.0);
  fitFunc->SetParLimits(12, 0.1, 10.0);
  fitFunc->SetParLimits(13, 0.1, 50.0);

  // --- ExpoPol4 バックグラウンド パラメータ (par[14..20]) ---
  // par[14..20] = p0..p6 of ExpoPol4
  // サイドバンドフィットの結果を初期値として使用
  for (int i = 0; i < 7; i++) {
    fitFunc->SetParameter(14 + i, fBkgSideband->GetParameter(i));
    fitFunc->SetParName(14 + i, Form("Bkg_p%d", i));
  }

  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineWidth(2);

  hSignal->Fit("fitFunc", "RM0");

  // =====================================================
  // 3. 描画（全範囲外挿）
  // =====================================================
  // 全体フィット関数（描画用）
  TF1* fitFuncDraw = new TF1("fitFuncDraw", TotalFit, drawMin, drawMax, 21);
  for (int i = 0; i < 21; i++) {
    fitFuncDraw->SetParameter(i, fitFunc->GetParameter(i));
  }
  fitFuncDraw->SetLineColor(kBlue);
  fitFuncDraw->SetLineWidth(2);
  fitFuncDraw->Draw("SAME");

  // ExpoPol4 バックグラウンド成分（描画用）
  TF1* fitBkgDraw = new TF1("fitBkgDraw", ExpoPol4, drawMin, drawMax, 7);
  for (int i = 0; i < 7; i++) {
    fitBkgDraw->SetParameter(i, fitFunc->GetParameter(14 + i));
  }
  fitBkgDraw->SetLineColor(TColor::GetColor("#e67e22")); // オレンジ
  fitBkgDraw->SetLineStyle(2);
  fitBkgDraw->SetLineWidth(2);
  fitBkgDraw->Draw("SAME");

  // J/psi DSCB シグナル成分（描画用）
  TF1* fitSigJpsi = new TF1("fitSigJpsi", DSCB, drawMin, drawMax, 7);
  for (int i = 0; i < 7; i++) {
    fitSigJpsi->SetParameter(i, fitFunc->GetParameter(i));
  }
  fitSigJpsi->SetLineColor(kMagenta);
  fitSigJpsi->SetLineStyle(2);
  fitSigJpsi->SetLineWidth(2);
  fitSigJpsi->Draw("SAME");

  // psi(2S) DSCB シグナル成分（描画用）
  TF1* fitSigPsi2s = new TF1("fitSigPsi2s", DSCB, drawMin, drawMax, 7);
  for (int i = 0; i < 7; i++) {
    fitSigPsi2s->SetParameter(i, fitFunc->GetParameter(7 + i));
  }
  fitSigPsi2s->SetLineColor(kRed);
  fitSigPsi2s->SetLineStyle(2);
  fitSigPsi2s->SetLineWidth(2);
  fitSigPsi2s->Draw("SAME");

  // Legend
  TLegend* leg2 = new TLegend(0.60, 0.60, 0.95, 0.95);
  leg2->SetBorderSize(1);
  leg2->AddEntry(hSignal,     "Data (Unlike-Like)",    "lep");
  leg2->AddEntry(fitFuncDraw, "Total Fit (DSCB+DSCB+ExpoPol4)", "l");
  leg2->AddEntry(fitBkgDraw,  "Background (ExpoPol4)", "l");
  leg2->AddEntry(fitSigJpsi,  "J/#psi Signal (DSCB)",  "l");
  leg2->AddEntry(fitSigPsi2s, "#psi(2S) Signal (DSCB)", "l");
  leg2->Draw();

  // Latex 情報表示
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.14, 0.88, Form("M_{J/#psi} = %.3f #pm %.3f GeV/c^{2}",
                                     fitFunc->GetParameter(1), fitFunc->GetParError(1)));
  latex->DrawLatex(0.14, 0.83, Form("#sigma_{J/#psi} = %.3f #pm %.3f GeV/c^{2}",
                                     fitFunc->GetParameter(2), fitFunc->GetParError(2)));
  latex->DrawLatex(0.14, 0.78, Form("M_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}",
                                     fitFunc->GetParameter(8), fitFunc->GetParError(8)));
  latex->DrawLatex(0.14, 0.73, Form("#sigma_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}",
                                     fitFunc->GetParameter(9), fitFunc->GetParError(9)));

  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);
  double nJpsi    = fitSigJpsi->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double nPsi2s   = fitSigPsi2s->Integral(fitRangeMin, fitRangeMax) / binWidth;
  latex->DrawLatex(0.14, 0.68, Form("N_{J/#psi} #approx %.0f,  N_{#psi(2S)} #approx %.0f", nJpsi, nPsi2s));

  // =====================================================
  // 4. Residual: (Data - Fit) / Fit
  // =====================================================
  pad2->cd();
  gPad->SetGridy();

  TH1D* hRatio = (TH1D*)hSignal->Clone("hRatio");
  hRatio->SetTitle("");

  for (int i = 1; i <= hRatio->GetNbinsX(); i++) {
    double xc      = hRatio->GetBinCenter(i);
    if (xc >= drawMin && xc <= drawMax) {
      double data    = hSignal->GetBinContent(i);
      double dataErr = hSignal->GetBinError(i);
      double fitVal  = fitFuncDraw->Eval(xc);
      if (fitVal > 0) {
        hRatio->SetBinContent(i, (data - fitVal) / fitVal);
        hRatio->SetBinError(i, dataErr / fitVal);
      } else {
        hRatio->SetBinContent(i, 0);
        hRatio->SetBinError(i, 0);
      }
    } else {
      hRatio->SetBinContent(i, 0);
      hRatio->SetBinError(i, 0);
    }
  }

  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->GetYaxis()->SetTitle("(Data-Fit)/Fit");
  hRatio->GetYaxis()->SetRangeUser(-1.5, 1.5);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetTitleOffset(0.4);
  hRatio->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/c^{2}]");
  hRatio->GetXaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->Draw("EP");

  TF1* lineZero = new TF1("lineZero", "0", drawMin, drawMax);
  lineZero->SetLineColor(kRed);
  lineZero->SetLineStyle(2);
  lineZero->Draw("SAME");

  c2->SaveAs("Unlike_bkg_WideRange_Massfit01.png");
  f->Close();
}
