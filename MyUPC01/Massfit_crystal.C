#include "Math/MinimizerOptions.h"
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
// Function definitions

// ==========================================

// ------------------------------------------
// Dimuon Continuum: Expo3
// f(m) = p0 * exp(p1*m + p2*m^2)

// ------------------------------------------

Double_t Expo3(Double_t* x, Double_t* par)
{
  double m = x[0];
  return par[0] * TMath::Exp(par[1] * m + par[2] * m * m);
}

Double_t Expo3Sideband(Double_t* x, Double_t* par)
{
  // exclude J/\psi and psi(2S) signal regions
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] * TMath::Exp(par[1] * x[0] + par[2] * x[0] * x[0]);
}

// ------------------------------------------
// Background: ExpoPol4
// f(m) = p0 * exp(p1*m) * (p2 + p3*m + p4*m^2 + p5*m^3 + p6*m^4)

// ------------------------------------------

Double_t ExpoPol4(Double_t* x, Double_t* par)
{
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m + par[5] * m * m * m + par[6] * m * m * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}
// ExpoPol4 side band fit

Double_t ExpoPol4Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m + par[5] * m * m * m + par[6] * m * m * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}

// ------------------------------------------
// Background: VMG

// ------------------------------------------

Double_t VMG(Double_t* x, Double_t* par)
{
  double m = x[0];
  double mbar = par[1];
  double sigma = par[2] + par[3] * (m - mbar) / mbar;
  if (sigma <= 0)
    return 0;
  return par[0] * TMath::Exp(-(m - mbar) * (m - mbar) / (2.0 * sigma * sigma));
}

// ------------------------------------------
// Double Sided Crystal Ball

// -------------------------------------------

Double_t DSCB(Double_t* x, Double_t* par)
{
  double N = par[0];
  double m0 = par[1];
  double sigma = par[2];
  double alpha1 = TMath::Abs(par[3]);
  double n1 = par[4];
  double alpha2 = TMath::Abs(par[5]);
  double n2 = par[6];
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

// ------------------------------------------
// Nominal Fit (Expo3)

// ------------------------------------------

Double_t TotalFit_Main(Double_t* x, Double_t* par)
{
  double Jpsi = DSCB(x, &par[0]);
  double psi2s = DSCB(x, &par[7]);
  double continuum = Expo3(x, &par[14]);
  return Jpsi + psi2s + continuum;
}

// ------------------------------------------
// Systematic Check 1 (ExpoPol4)

// ------------------------------------------

Double_t TotalFit_Sys_ExpoPol4(Double_t* x, Double_t* par)
{
  double Jpsi = DSCB(x, &par[0]);
  double psi2s = DSCB(x, &par[7]);
  double bkg = ExpoPol4(x, &par[14]);
  return Jpsi + psi2s + bkg;
}

// ------------------------------------------
// Systematic Check 2 (VMG)

// ------------------------------------------

Double_t TotalFit_Sys_VWG(Double_t* x, Double_t* par)
{
  double Jpsi = DSCB(x, &par[0]);
  double psi2s = DSCB(x, &par[7]);
  double bkg = VMG(x, &par[14]);
  return Jpsi + psi2s + bkg;
}

// ==========================================
// Main macro

// ==========================================

void Massfit_crystal()
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFontSize(0.04);
  TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/AnalysisResults.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open Analysis file" << std::endl;
    return;
  }
  TString dir = "my-upc-mass-01/";
  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dir + "MMuonUnlike"));
  TH1* hLike = dynamic_cast<TH1*>(f->Get(dir + "MMuonLike"));
  if (!hUnlike || !hLike) {
    std::cerr << "ERROR: Cannot get histograms" << std::endl;
    return;
  }
  double drawMin = 1.8;
  double drawMax = 7.0;
  // ========================================
  // Raw Histograms for checking
  // ========================================
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
  hLike->Draw("E Same");
  TLegend* leg = new TLegend(0.65, 0.80, 0.95, 0.95);
  leg->SetBorderSize(1);
  leg->AddEntry(hUnlike, "Unlike Sign (+-)", "lep");
  leg->AddEntry(hLike, "Like Sign (++ & --)", "lep");
  leg->Draw();
  c1->SaveAs("RawMass_ALL_forChecke.png");
  // =======================================
  // Signal (Unlike -Like)
  // =======================================
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->SetTitle("Signal (Unlike - Like) ALL; M_{#mu#mu} [GeV/c^{2}]; Counts");
  hSignal->Add(hLike, -1.0);
  TCanvas* c2 = new TCanvas("c2", "Signal Extraction and Mass Fit (Massfit DSCB)", 800, 800);
  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);
  pad1->SetTopMargin(0.05);
  pad1->Draw();
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.30);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);
  pad2->Draw();
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
  hSignal->GetXaxis()->SetLabelSize(0);
  hSignal->GetXaxis()->SetTitleSize(0);
  hSignal->Draw("E");
  double fitRangeMin = 1.8;
  double fitRangeMax = 7.0;
  // =======================================
  // ExpoPol4 Sideband Fit
  // =======================================
  TF1* fBkgSideband = new TF1("fBkgSideband", ExpoPol4Sideband, fitRangeMin, fitRangeMax, 7);
  double initBkgScale = hSignal->GetBinContent(hSignal->FindBin(fitRangeMin));
  if (initBkgScale <= 0)
    initBkgScale = 100.0;
  fBkgSideband->SetParameters(initBkgScale, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
  fBkgSideband->SetParNames("BkgP0", "BkgP1", "BkgP2", "BkgP3", "BkgP4", "BkgP5", "BkgP6");
  fBkgSideband->FixParameter(2, 1.0);
  hSignal->Fit("fBkgSideband", "R0");
  // =======================================
  // Expo3 Sideband Fit
  // =======================================
  TF1* fExpo3Sideband = new TF1("fExpo3Sideband", Expo3Sideband, fitRangeMin, fitRangeMax, 3);
  double initExpo3Scale = hSignal->GetBinContent(hSignal->FindBin(fitRangeMin));
  if (initExpo3Scale <= 0)
    initExpo3Scale = 100.0;
  fExpo3Sideband->SetParameters(initExpo3Scale, -1.0, 0.01);
  hSignal->Fit("fExpo3Sideband", "R0");
  // =======================================
  // J/\psi DSCB + psi(2S) DSCB + Expo3 bkg
  // =======================================
  hSignal->GetXaxis()->SetRangeUser(2.9, 3.3);
  double maxValJpsi = hSignal->GetBinContent(hSignal->GetMaximumBin());
  hSignal->GetXaxis()->SetRangeUser(3.55, 3.8);
  double maxValPsi2s = hSignal->GetBinContent(hSignal->GetMaximumBin());
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);
  if (maxValJpsi <= 0)
    maxValJpsi = 100.0;
  if (maxValPsi2s <= 0)
    maxValPsi2s = 10.0;
  TF1* fitFunc = new TF1("fitFunc", TotalFit_Main, fitRangeMin, fitRangeMax, 17);
  // ----------------------------------------
  // J/psi DSCB
  // ----------------------------------------
  fitFunc->SetParameter(0, maxValJpsi);
  fitFunc->SetParameter(1, 3.096);
  fitFunc->SetParameter(2, 0.05);
  fitFunc->SetParLimits(0, 0.0, maxValJpsi * 5.0);
  fitFunc->SetParLimits(1, 2.9, 3.25);
  fitFunc->SetParLimits(2, 0.02, 0.15);
  // 【最重要修正】テールパラメータを固定（Fix）して縮退を防ぐ
  fitFunc->FixParameter(3, 1.5); // alpha1
  fitFunc->FixParameter(4, 3.0); // n1
  fitFunc->FixParameter(5, 1.5); // alpha2
  fitFunc->FixParameter(6, 3.0); // n2
  fitFunc->SetParName(0, "Jpsi_N");
  fitFunc->SetParName(1, "Jpsi_m0");
  fitFunc->SetParName(2, "Jpsi_sigma");
  fitFunc->SetParName(3, "Jpsi_alpha1");
  fitFunc->SetParName(4, "Jpsi_n1");
  fitFunc->SetParName(5, "Jpsi_alpha2");
  fitFunc->SetParName(6, "Jpsi_n2");
  // -------------------------------------------
  // psi(2S) DSCB
  // -------------------------------------------
  fitFunc->SetParameter(7, maxValPsi2s);
  fitFunc->SetParameter(8, 3.686);
  fitFunc->SetParameter(9, 0.05 * (3.686 / 3.096));
  fitFunc->SetParLimits(7, 0.0, maxValJpsi);
  fitFunc->SetParLimits(8, 3.65, 3.72);
  fitFunc->SetParLimits(9, 0.03, 0.08);
  // テールパラメータをJ/psiと同じ値に固定
  fitFunc->FixParameter(10, 1.5); // alpha1
  fitFunc->FixParameter(11, 3.0); // n1
  fitFunc->FixParameter(12, 1.5); // alpha2
  fitFunc->FixParameter(13, 3.0); // n2
  fitFunc->SetParName(7, "Psi2S_N");
  fitFunc->SetParName(8, "Psi2S_m0");
  fitFunc->SetParName(9, "Psi2S_sigma");
  fitFunc->SetParName(10, "Psi2S_alpha1");
  fitFunc->SetParName(11, "Psi2S_n1");
  fitFunc->SetParName(12, "Psi2S_alpha2");
  fitFunc->SetParName(13, "Psi2S_n2");
  // -------------------------------------------
  // Dimuon Continuum (Expo3)
  // 事前フィット(Sideband)の結果を初期値として渡す
  // -------------------------------------------
  fitFunc->SetParameter(14, fExpo3Sideband->GetParameter(0));
  fitFunc->SetParameter(15, fExpo3Sideband->GetParameter(1));
  fitFunc->SetParameter(16, fExpo3Sideband->GetParameter(2));
  fitFunc->SetParName(14, "Cont_p0");
  fitFunc->SetParName(15, "Cont_p1");
  fitFunc->SetParName(16, "Cont_p2");
  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineWidth(2);
  // Perform the nominal fit
  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Nominal Fit (TotalFit_Main) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc", "RM0");
  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);
  TF1* fitSigJpsi_Main = new TF1("fitSigJpsi_Main", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigJpsi_Main->SetParameter(i, fitFunc->GetParameter(i));
  double yieldJpsi_Main = fitSigJpsi_Main->Integral(fitRangeMin, fitRangeMax) / binWidth;
  TF1* fitSigPsi2s_Main = new TF1("fitSigPsi2s_Main", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigPsi2s_Main->SetParameter(i, fitFunc->GetParameter(7 + i));
  double yieldPsi2s_Main = fitSigPsi2s_Main->Integral(fitRangeMin, fitRangeMax) / binWidth;
  std::cout << "Yield J/psi (Nominal): " << yieldJpsi_Main << std::endl;
  std::cout << "Yield psi(2S) (Nominal): " << yieldPsi2s_Main << std::endl;
  // -------------------------------------------
  // Systematic Check 1: ExpoPol4 Background
  // -------------------------------------------
  TF1* fitFunc_Sys1 = new TF1("fitFunc_Sys1", TotalFit_Sys_ExpoPol4, fitRangeMin, fitRangeMax, 21);
  for (int i = 0; i < 14; i++) {
    fitFunc_Sys1->SetParameter(i, fitFunc->GetParameter(i));
    double pmin, pmax;
    fitFunc->GetParLimits(i, pmin, pmax);
    if (pmin < pmax)
      fitFunc_Sys1->SetParLimits(i, pmin, pmax);
  }
  // Sys1でもテールを固定
  fitFunc_Sys1->FixParameter(3, 1.5);
  fitFunc_Sys1->FixParameter(4, 3.0);
  fitFunc_Sys1->FixParameter(5, 1.5);
  fitFunc_Sys1->FixParameter(6, 3.0);
  fitFunc_Sys1->FixParameter(10, 1.5);
  fitFunc_Sys1->FixParameter(11, 3.0);
  fitFunc_Sys1->FixParameter(12, 1.5);
  fitFunc_Sys1->FixParameter(13, 3.0);
  for (int i = 0; i < 7; i++) {
    fitFunc_Sys1->SetParameter(14 + i, fBkgSideband->GetParameter(i));
  }
  fitFunc_Sys1->FixParameter(14 + 2, 1.0); // p2固定を引き継ぎ
  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 1 (ExpoPol4) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys1", "RM0");
  TF1* fitSigJpsi_Sys1 = new TF1("fitSigJpsi_Sys1", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigJpsi_Sys1->SetParameter(i, fitFunc_Sys1->GetParameter(i));
  double yieldJpsi_Sys1 = fitSigJpsi_Sys1->Integral(fitRangeMin, fitRangeMax) / binWidth;
  TF1* fitSigPsi2s_Sys1 = new TF1("fitSigPsi2s_Sys1", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigPsi2s_Sys1->SetParameter(i, fitFunc_Sys1->GetParameter(7 + i));
  double yieldPsi2s_Sys1 = fitSigPsi2s_Sys1->Integral(fitRangeMin, fitRangeMax) / binWidth;
  // -------------------------------------------
  // Systematic Check 2: VMG Background
  // -------------------------------------------
  TF1* fitFunc_Sys2 = new TF1("fitFunc_Sys2", TotalFit_Sys_VWG, fitRangeMin, fitRangeMax, 18);
  for (int i = 0; i < 14; i++) {
    fitFunc_Sys2->SetParameter(i, fitFunc->GetParameter(i));
    double pmin, pmax;
    fitFunc->GetParLimits(i, pmin, pmax);
    if (pmin < pmax)
      fitFunc_Sys2->SetParLimits(i, pmin, pmax);
  }
  // Sys2でもテールを固定
  fitFunc_Sys2->FixParameter(3, 1.5);
  fitFunc_Sys2->FixParameter(4, 3.0);
  fitFunc_Sys2->FixParameter(5, 1.5);
  fitFunc_Sys2->FixParameter(6, 3.0);
  fitFunc_Sys2->FixParameter(10, 1.5);
  fitFunc_Sys2->FixParameter(11, 3.0);
  fitFunc_Sys2->FixParameter(12, 1.5);
  fitFunc_Sys2->FixParameter(13, 3.0);
  fitFunc_Sys2->SetParameter(14, 100.0);
  fitFunc_Sys2->SetParameter(15, 3.1);
  fitFunc_Sys2->SetParameter(16, 1.0);
  fitFunc_Sys2->SetParameter(17, 0.0);
  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 2 (VMG) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys2", "RM0");
  TF1* fitSigJpsi_Sys2 = new TF1("fitSigJpsi_Sys2", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigJpsi_Sys2->SetParameter(i, fitFunc_Sys2->GetParameter(i));
  double yieldJpsi_Sys2 = fitSigJpsi_Sys2->Integral(fitRangeMin, fitRangeMax) / binWidth;
  TF1* fitSigPsi2s_Sys2 = new TF1("fitSigPsi2s_Sys2", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigPsi2s_Sys2->SetParameter(i, fitFunc_Sys2->GetParameter(7 + i));
  double yieldPsi2s_Sys2 = fitSigPsi2s_Sys2->Integral(fitRangeMin, fitRangeMax) / binWidth;
  std::cout << "\n------------- YIELD SUMMARY --------------\n";
  std::cout << "J/psi Yield: Nominal = " << yieldJpsi_Main << ", ExpoPol4 = " << yieldJpsi_Sys1 << ", VMG = " << yieldJpsi_Sys2 << "\n";
  std::cout << "psi(2S) Yield: Nominal = " << yieldPsi2s_Main << ", ExpoPol4 = " << yieldPsi2s_Sys1 << ", VMG = " << yieldPsi2s_Sys2 << "\n";
  std::cout << "------------------------------------------\n\n";
  // ===========================================
  // Draw (Using Nominal Fit)
  // ===========================================
  TF1* fitFuncDraw = new TF1("fitFUncDraw", TotalFit_Main, drawMin, drawMax, 17);
  for (int i = 0; i < 17; i++)
    fitFuncDraw->SetParameter(i, fitFunc->GetParameter(i));
  fitFuncDraw->SetLineColor(kBlue);
  fitFuncDraw->SetLineWidth(2);
  fitFuncDraw->Draw("SAME");
  TF1* fitContDraw = new TF1("fitContDraw", Expo3, drawMin, drawMax, 3);
  for (int i = 0; i < 3; i++)
    fitContDraw->SetParameter(i, fitFunc->GetParameter(14 + i));
  fitContDraw->SetLineColor(kCyan);
  fitContDraw->SetLineStyle(3);
  fitContDraw->Draw("SAME");
  TF1* fitSigJpsi = new TF1("fitSigJpsi", DSCB, drawMin, drawMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigJpsi->SetParameter(i, fitFunc->GetParameter(i));
  fitSigJpsi->SetLineColor(kMagenta);
  fitSigJpsi->SetLineStyle(2);
  fitSigJpsi->SetLineWidth(2);
  fitSigJpsi->Draw("SAME");
  TF1* fitSigPsi2s = new TF1("fitSigPsi2s", DSCB, drawMin, drawMax, 7);
  for (int i = 0; i < 7; i++)
    fitSigPsi2s->SetParameter(i, fitFunc->GetParameter(7 + i));
  fitSigPsi2s->SetLineColor(kRed);
  fitSigPsi2s->SetLineStyle(2);
  fitSigPsi2s->SetLineWidth(2);
  fitSigPsi2s->Draw("SAME");
  TLegend* leg2 = new TLegend(0.60, 0.60, 0.95, 0.95);
  leg2->SetBorderSize(1);
  leg2->AddEntry(hSignal, "Data (Unlike-Like)", "lep");
  leg2->AddEntry(fitFuncDraw, "Total Fit (DSCB+DSCB+Expo3)", "l");
  leg2->AddEntry(fitSigJpsi, "J/#psi Signal (DSCB)", "l");
  leg2->AddEntry(fitSigPsi2s, "#psi(2S) Signal (DSCB)", "l");
  leg2->AddEntry(fitContDraw, "#gamma#gamma #rightarrow #mu^{+}#mu^{-} Dimuon Continuum (Expo3)", "l");
  leg2->Draw();
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.14, 0.88, Form("M_{J/#psi} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
  latex->DrawLatex(0.14, 0.83, Form("#sigma_{J/#psi} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(2), fitFunc->GetParError(2)));
  latex->DrawLatex(0.14, 0.78, Form("M_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(8), fitFunc->GetParError(8)));
  latex->DrawLatex(0.14, 0.73, Form("#sigma_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}", fitFunc->GetParameter(9), fitFunc->GetParError(9)));
  latex->DrawLatex(0.14, 0.68, Form("N_{J/#psi} #approx %.0f, N_{#psi(2S)} #approx %.0f", yieldJpsi_Main, yieldPsi2s_Main));
  // ============================================
  // Residual: (Data -Fit) / Fit
  // ============================================
  pad2->cd();
  gPad->SetGridy();
  TH1D* hRatio = (TH1D*)hSignal->Clone("hRatio");
  hRatio->SetTitle("");
  for (int i = 1; i <= hRatio->GetNbinsX(); i++) {
    double xc = hRatio->GetBinCenter(i);
    if (xc >= drawMin && xc <= drawMax) {
      double data = hSignal->GetBinContent(i);
      double dataErr = hSignal->GetBinError(i);
      double fitVal = fitFuncDraw->Eval(xc);
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
  c2->SaveAs("Unlike_bkg_WideRange_crystal.png");
  f->Close();
}
