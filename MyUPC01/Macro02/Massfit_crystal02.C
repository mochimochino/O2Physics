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
// Global Constants & MC Parameters
// ==========================================
const double PDG_M_JPSI = 3.09691;
const double PDG_M_PSI2S = 3.68609;
const double DELTA_M = PDG_M_PSI2S - PDG_M_JPSI;

const double MC_WIDTH_RATIO = 1.0;

// J/psi Tails (MC)
const double MC_JPSI_A1 = 1.0; // Left tail alpha
const double MC_JPSI_N1 = 2.0; // Left tail n
const double MC_JPSI_A2 = 2.5; // Right tail alpha (Negligible impact, fixed)
const double MC_JPSI_N2 = 3.0; // Right tail n (Negligible impact, fixed)

// psi(2S) Tails (MC) -> All fixed
const double MC_PSI2S_A1 = 1.0;
const double MC_PSI2S_N1 = 2.0;
const double MC_PSI2S_A2 = 2.5;
const double MC_PSI2S_N2 = 3.0;

// ==========================================
// Function definitions
// ==========================================

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
// Background Functions
// ------------------------------------------
Double_t Expo3(Double_t* x, Double_t* par)
{
  return par[0] * TMath::Exp(par[1] * x[0] + par[2] * x[0] * x[0]);
}
Double_t Expo3Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] * TMath::Exp(par[1] * x[0] + par[2] * x[0] * x[0]);
}

Double_t ExpoPol4(Double_t* x, Double_t* par)
{
  double poly = par[2] + par[3] * x[0] + par[4] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] * x[0];
  return par[0] * TMath::Exp(par[1] * x[0]) * poly;
}
Double_t ExpoPol4Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  return ExpoPol4(x, par);
}

Double_t VWG(Double_t* x, Double_t* par)
{
  double mbar = par[1];
  double sigma = par[2] + par[3] * (x[0] - mbar) / mbar;
  if (sigma <= 0)
    return 0;
  return par[0] * TMath::Exp(-(x[0] - mbar) * (x[0] - mbar) / (2.0 * sigma * sigma));
}

Double_t ExpoPol2(Double_t* x, Double_t* par)
{
  double poly = par[2] + par[3] * x[0] + par[4] * x[0] * x[0];
  return par[0] * TMath::Exp(par[1] * x[0]) * poly;
}
Double_t ExpoPol2Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint();
    return 0;
  }
  return ExpoPol2(x, par);
}

// ------------------------------------------
// Total Fit Functions (Strict Parameter Tying)
// ------------------------------------------
// Parameter mapping (0-12):
// 0: N_Jpsi, 1: m_Jpsi, 2: sigma_Jpsi, 3: a1_Jpsi, 4: n1_Jpsi, 5: a2_Jpsi, 6: n2_Jpsi
// 7: N_psi2s, 8: a1_psi2s, 9: n1_psi2s, 10: a2_psi2s, 11: n2_psi2s
// 12: MC_WIDTH_RATIO (Fixed, varied for systematics)
// 13+: Background parameters
void SetTiedParameters(Double_t* par, double* parJpsi, double* parPsi2s)
{
  for (int i = 0; i < 7; i++)
    parJpsi[i] = par[i];
  double width_ratio = par[12];
  parPsi2s[0] = par[7];
  parPsi2s[1] = par[1] + DELTA_M;           // Tied Mass
  parPsi2s[2] = par[2] * 1.1 * width_ratio; // Tied Sigma
  parPsi2s[3] = par[8];
  parPsi2s[4] = par[9];
  parPsi2s[5] = par[10];
  parPsi2s[6] = par[11];
}

Double_t TotalFit_Main(Double_t* x, Double_t* par)
{
  double parJpsi[7], parPsi2s[7];
  SetTiedParameters(par, parJpsi, parPsi2s);
  return DSCB(x, parJpsi) + DSCB(x, parPsi2s) + Expo3(x, &par[13]);
}
Double_t TotalFit_Sys_ExpoPol4(Double_t* x, Double_t* par)
{
  double parJpsi[7], parPsi2s[7];
  SetTiedParameters(par, parJpsi, parPsi2s);
  return DSCB(x, parJpsi) + DSCB(x, parPsi2s) + ExpoPol4(x, &par[13]);
}
Double_t TotalFit_Sys_VWG(Double_t* x, Double_t* par)
{
  double parJpsi[7], parPsi2s[7];
  SetTiedParameters(par, parJpsi, parPsi2s);
  return DSCB(x, parJpsi) + DSCB(x, parPsi2s) + VWG(x, &par[13]);
}
Double_t TotalFit_Sys_ExpoPol2(Double_t* x, Double_t* par)
{
  double parJpsi[7], parPsi2s[7];
  SetTiedParameters(par, parJpsi, parPsi2s);
  return DSCB(x, parJpsi) + DSCB(x, parPsi2s) + ExpoPol2(x, &par[13]);
}

// ==========================================
// Main macro
// ==========================================
void Massfit_crystal02()
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleFontSize(0.04);

  TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/0318/AnalysisResults.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open Analysis file\n";
    return;
  }

  TString dir = "my-upc-mass-02/registry/";
  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dir + "hMassUnlike"));
  if (!hUnlike) {
    std::cerr << "ERROR: Cannot get histogram\n";
    return;
  }

  hUnlike->Rebin(25); // 10 MeV -> 25 MeV

  double drawMin = 2.0;
  double drawMax = 4.5;
  double fitRangeMin = 2.5;
  double fitRangeMax = 4.5;

  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->SetTitle("Signal (Unlike) ALL; m_{#mu#mu} [GeV/c^{2}]; Counts");
  double ymin = hSignal->GetMinimum();
  if (ymin > 0)
    ymin = 0;
  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);
  hSignal->GetYaxis()->SetRangeUser(ymin * 1.2, hSignal->GetMaximum() * 1.5);

  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);

  // --- Raw Continuum Outside Fit Region ---
  int binMinFit = hSignal->GetXaxis()->FindBin(fitRangeMin);
  int binMaxFit = hSignal->GetXaxis()->FindBin(fitRangeMax) - 1;
  double rawYieldLow = hSignal->Integral(1, binMinFit - 1);
  double rawYieldHigh = hSignal->Integral(binMaxFit + 1, hSignal->GetNbinsX());
  double rawContinuumOutside = rawYieldLow + rawYieldHigh;
  double errRawContinuumOutside = TMath::Sqrt(rawContinuumOutside);

  // --- Sideband Fits ---
  double initScale = hSignal->GetBinContent(hSignal->FindBin(fitRangeMin));
  if (initScale <= 0)
    initScale = 100.0;

  TF1* fExpo3Sideband = new TF1("fExpo3Sideband", Expo3Sideband, fitRangeMin, fitRangeMax, 3);
  fExpo3Sideband->SetParameters(initScale, -1.0, 0.01);
  hSignal->Fit("fExpo3Sideband", "R0Q");

  TF1* fExpoPol4Sideband = new TF1("fExpoPol4Sideband", ExpoPol4Sideband, fitRangeMin, fitRangeMax, 7);
  fExpoPol4Sideband->SetParameters(initScale, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
  fExpoPol4Sideband->FixParameter(2, 1.0);
  hSignal->Fit("fExpoPol4Sideband", "R0Q");

  TF1* fExpoPol2Sideband = new TF1("fExpoPol2Sideband", ExpoPol2Sideband, fitRangeMin, fitRangeMax, 5);
  fExpoPol2Sideband->SetParameters(initScale, -1.0, 1.0, 0.0, 0.0);
  fExpoPol2Sideband->FixParameter(2, 1.0);
  hSignal->Fit("fExpoPol2Sideband", "R0Q");

  double maxValJpsi = hSignal->GetMaximum();

  // Helper lambda to setup base parameters for any TotalFit
  auto setupBaseParams = [&](TF1* f) {
    f->SetParameter(0, maxValJpsi);
    f->SetParLimits(0, 0.0, maxValJpsi * 5.0); // N_Jpsi
    f->SetParameter(1, PDG_M_JPSI);
    f->SetParLimits(1, 3.05, 3.15); // m_Jpsi (Free)
    f->SetParameter(2, 0.08);
    f->SetParLimits(2, 0.05, 0.12); // sigma_Jpsi (Free)

    f->SetParameter(3, MC_JPSI_A1); // J/psi a1 (Free in Nominal)
    f->SetParameter(4, MC_JPSI_N1); // J/psi n1 (Free in Nominal)

    f->FixParameter(5, MC_JPSI_A2); // J/psi a2 (Fixed)
    f->FixParameter(6, MC_JPSI_N2); // J/psi n2 (Fixed)

    f->SetParameter(7, maxValJpsi * 0.02);
    f->SetParLimits(7, 0.0, maxValJpsi); // N_psi2s

    f->FixParameter(8, MC_PSI2S_A1); // psi(2S) tails (All Fixed)
    f->FixParameter(9, MC_PSI2S_N1);
    f->FixParameter(10, MC_PSI2S_A2);
    f->FixParameter(11, MC_PSI2S_N2);

    f->FixParameter(12, MC_WIDTH_RATIO); // Width Ratio (Fixed in Nominal)
  };

  // =======================================
  // 0. Nominal Fit (Expo3)
  // =======================================
  TF1* fitFunc = new TF1("fitFunc", TotalFit_Main, fitRangeMin, fitRangeMax, 16);
  setupBaseParams(fitFunc);
  fitFunc->SetParameter(13, fExpo3Sideband->GetParameter(0));
  fitFunc->SetParameter(14, fExpo3Sideband->GetParameter(1));
  fitFunc->SetParameter(15, fExpo3Sideband->GetParameter(2));

  std::cout << "\n=== Nominal Fit (Expo3) ===" << std::endl;
  hSignal->Fit("fitFunc", "RM0");
  double yieldJpsi_Main = fitFunc->GetParameter(0) / binWidth;
  double yieldPsi2s_Main = fitFunc->GetParameter(7) / binWidth;
  TF1* fBkg0 = new TF1("fBkg0", Expo3, drawMin, drawMax, 3);
  for (int i = 0; i < 3; i++)
    fBkg0->SetParameter(i, fitFunc->GetParameter(13 + i));
  double totalContinuum_Main = rawContinuumOutside + fBkg0->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Sys1 Fit (ExpoPol4)
  // =======================================
  TF1* fitFunc_Sys1 = new TF1("fitFunc_Sys1", TotalFit_Sys_ExpoPol4, fitRangeMin, fitRangeMax, 20);
  setupBaseParams(fitFunc_Sys1);
  for (int i = 0; i < 7; i++) {
    double sb_val = fExpoPol4Sideband->GetParameter(i);
    fitFunc_Sys1->SetParameter(13 + i, sb_val);
    if (i > 2)
      fitFunc_Sys1->SetParLimits(13 + i, sb_val - std::abs(sb_val) * 5.0 - 0.1, sb_val + std::abs(sb_val) * 5.0 + 0.1);
  }
  fitFunc_Sys1->FixParameter(13 + 2, 1.0);

  std::cout << "\n=== Sys1: ExpoPol4 ===" << std::endl;
  hSignal->Fit("fitFunc_Sys1", "RM0");
  double yieldJpsi_Sys1 = fitFunc_Sys1->GetParameter(0) / binWidth;
  double yieldPsi2s_Sys1 = fitFunc_Sys1->GetParameter(7) / binWidth;
  TF1* fBkg1 = new TF1("fBkg1", ExpoPol4, drawMin, drawMax, 7);
  for (int i = 0; i < 7; i++)
    fBkg1->SetParameter(i, fitFunc_Sys1->GetParameter(13 + i));
  double totalContinuum_Sys1 = rawContinuumOutside + fBkg1->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Sys2 Fit (VWG)
  // =======================================
  TF1* fitFunc_Sys2 = new TF1("fitFunc_Sys2", TotalFit_Sys_VWG, fitRangeMin, fitRangeMax, 17);
  setupBaseParams(fitFunc_Sys2);
  fitFunc_Sys2->SetParameter(13, 100.0);
  fitFunc_Sys2->SetParameter(14, 3.1);
  fitFunc_Sys2->SetParameter(15, 1.0);
  fitFunc_Sys2->SetParameter(16, 0.0);

  std::cout << "\n=== Sys2: VWG ===" << std::endl;
  hSignal->Fit("fitFunc_Sys2", "RM0");
  double yieldJpsi_Sys2 = fitFunc_Sys2->GetParameter(0) / binWidth;
  double yieldPsi2s_Sys2 = fitFunc_Sys2->GetParameter(7) / binWidth;
  TF1* fBkg2 = new TF1("fBkg2", VWG, drawMin, drawMax, 4);
  for (int i = 0; i < 4; i++)
    fBkg2->SetParameter(i, fitFunc_Sys2->GetParameter(13 + i));
  double totalContinuum_Sys2 = rawContinuumOutside + fBkg2->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Sys3 Fit (ExpoPol2)
  // =======================================
  TF1* fitFunc_Sys3 = new TF1("fitFunc_Sys3", TotalFit_Sys_ExpoPol2, fitRangeMin, fitRangeMax, 18);
  setupBaseParams(fitFunc_Sys3);
  for (int i = 0; i < 5; i++)
    fitFunc_Sys3->SetParameter(13 + i, fExpoPol2Sideband->GetParameter(i));
  fitFunc_Sys3->FixParameter(13 + 2, 1.0);

  std::cout << "\n=== Sys3: ExpoPol2 ===" << std::endl;
  hSignal->Fit("fitFunc_Sys3", "RM0");
  double yieldJpsi_Sys3 = fitFunc_Sys3->GetParameter(0) / binWidth;
  double yieldPsi2s_Sys3 = fitFunc_Sys3->GetParameter(7) / binWidth;
  TF1* fBkg3 = new TF1("fBkg3", ExpoPol2, drawMin, drawMax, 5);
  for (int i = 0; i < 5; i++)
    fBkg3->SetParameter(i, fitFunc_Sys3->GetParameter(13 + i));
  double totalContinuum_Sys3 = rawContinuumOutside + fBkg3->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Sys4 Fit (J/psi Left Tails Fixed to MC)
  // =======================================
  TF1* fitFunc_Sys4 = new TF1("fitFunc_Sys4", TotalFit_Main, fitRangeMin, fitRangeMax, 16);
  setupBaseParams(fitFunc_Sys4);
  fitFunc_Sys4->FixParameter(3, MC_JPSI_A1); // FIXED for Systematic check
  fitFunc_Sys4->FixParameter(4, MC_JPSI_N1); // FIXED for Systematic check
  for (int i = 0; i < 3; i++)
    fitFunc_Sys4->SetParameter(13 + i, fExpo3Sideband->GetParameter(i));

  std::cout << "\n=== Sys4: J/psi Left Tails Fixed ===" << std::endl;
  hSignal->Fit("fitFunc_Sys4", "RM0");
  double yieldJpsi_Sys4 = fitFunc_Sys4->GetParameter(0) / binWidth;
  double yieldPsi2s_Sys4 = fitFunc_Sys4->GetParameter(7) / binWidth;
  TF1* fBkg4 = new TF1("fBkg4", Expo3, drawMin, drawMax, 3);
  for (int i = 0; i < 3; i++)
    fBkg4->SetParameter(i, fitFunc_Sys4->GetParameter(13 + i));
  double totalContinuum_Sys4 = rawContinuumOutside + fBkg4->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Sys5 Fit (Width Ratio + 0.1)
  // =======================================
  TF1* fitFunc_Sys5 = new TF1("fitFunc_Sys5", TotalFit_Main, fitRangeMin, fitRangeMax, 16);
  setupBaseParams(fitFunc_Sys5);
  fitFunc_Sys5->FixParameter(12, MC_WIDTH_RATIO + 0.1); // Systematic Variation
  for (int i = 0; i < 3; i++)
    fitFunc_Sys5->SetParameter(13 + i, fExpo3Sideband->GetParameter(i));

  std::cout << "\n=== Sys5: Width Ratio + 0.1 ===" << std::endl;
  hSignal->Fit("fitFunc_Sys5", "RM0");
  double yieldJpsi_Sys5 = fitFunc_Sys5->GetParameter(0) / binWidth;
  double yieldPsi2s_Sys5 = fitFunc_Sys5->GetParameter(7) / binWidth;
  TF1* fBkg5 = new TF1("fBkg5", Expo3, drawMin, drawMax, 3);
  for (int i = 0; i < 3; i++)
    fBkg5->SetParameter(i, fitFunc_Sys5->GetParameter(13 + i));
  double totalContinuum_Sys5 = rawContinuumOutside + fBkg5->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Sys6 Fit (Width Ratio - 0.1)
  // =======================================
  TF1* fitFunc_Sys6 = new TF1("fitFunc_Sys6", TotalFit_Main, fitRangeMin, fitRangeMax, 16);
  setupBaseParams(fitFunc_Sys6);
  fitFunc_Sys6->FixParameter(12, MC_WIDTH_RATIO - 0.1); // Systematic Variation
  for (int i = 0; i < 3; i++)
    fitFunc_Sys6->SetParameter(13 + i, fExpo3Sideband->GetParameter(i));

  std::cout << "\n=== Sys6: Width Ratio - 0.1 ===" << std::endl;
  hSignal->Fit("fitFunc_Sys6", "RM0");
  double yieldJpsi_Sys6 = fitFunc_Sys6->GetParameter(0) / binWidth;
  double yieldPsi2s_Sys6 = fitFunc_Sys6->GetParameter(7) / binWidth;
  TF1* fBkg6 = new TF1("fBkg6", Expo3, drawMin, drawMax, 3);
  for (int i = 0; i < 3; i++)
    fBkg6->SetParameter(i, fitFunc_Sys6->GetParameter(13 + i));
  double totalContinuum_Sys6 = rawContinuumOutside + fBkg6->Integral(fitRangeMin, fitRangeMax) / binWidth;

  // =======================================
  // Yield Summary Output
  // =======================================
  std::cout << "\n------------- YIELD SUMMARY --------------\n";
  std::cout << "--- J/psi Yield ---" << std::endl;
  std::cout << Form(" Nominal (Expo3)     = %.1f \n", yieldJpsi_Main);
  std::cout << Form(" Sys1 (ExpoPol4)     = %.1f \n", yieldJpsi_Sys1);
  std::cout << Form(" Sys2 (VWG)          = %.1f \n", yieldJpsi_Sys2);
  std::cout << Form(" Sys3 (ExpoPol2)     = %.1f \n", yieldJpsi_Sys3);
  std::cout << Form(" Sys4 (J/psi Tails Fix)= %.1f \n", yieldJpsi_Sys4);
  std::cout << Form(" Sys5 (WidthRatio+0.1)= %.1f \n", yieldJpsi_Sys5);
  std::cout << Form(" Sys6 (WidthRatio-0.1)= %.1f \n", yieldJpsi_Sys6);

  std::cout << "\n--- psi(2S) Yield ---" << std::endl;
  std::cout << Form(" Nominal (Expo3)     = %.1f \n", yieldPsi2s_Main);
  std::cout << Form(" Sys1 (ExpoPol4)     = %.1f \n", yieldPsi2s_Sys1);
  std::cout << Form(" Sys2 (VWG)          = %.1f \n", yieldPsi2s_Sys2);
  std::cout << Form(" Sys3 (ExpoPol2)     = %.1f \n", yieldPsi2s_Sys3);
  std::cout << Form(" Sys4 (J/psi Tails Fix)= %.1f \n", yieldPsi2s_Sys4);
  std::cout << Form(" Sys5 (WidthRatio+0.1)= %.1f \n", yieldPsi2s_Sys5);
  std::cout << Form(" Sys6 (WidthRatio-0.1)= %.1f \n", yieldPsi2s_Sys6);

  std::cout << "\n--- Continuum Yield (Raw Outside + Fit Inside) ---" << std::endl;
  std::cout << Form(" Nominal (Expo3)     = %.1f \n", totalContinuum_Main);
  std::cout << Form(" Sys1 (ExpoPol4)     = %.1f \n", totalContinuum_Sys1);
  std::cout << Form(" Sys2 (VWG)          = %.1f \n", totalContinuum_Sys2);
  std::cout << Form(" Sys3 (ExpoPol2)     = %.1f \n", totalContinuum_Sys3);
  std::cout << Form(" Sys4 (J/psi Tails Fix)= %.1f \n", totalContinuum_Sys4);
  std::cout << Form(" Sys5 (WidthRatio+0.1)= %.1f \n", totalContinuum_Sys5);
  std::cout << Form(" Sys6 (WidthRatio-0.1)= %.1f \n", totalContinuum_Sys6);
  std::cout << "------------------------------------------\n\n";

  // ===========================================
  // Draw Macros Helper Function
  // ===========================================
  auto drawCanvas = [&](TString cname, TString title, TF1* fTot, TF1* fBkg, TString bkgName, TString outName) {
    TCanvas* c = new TCanvas(cname, title, 800, 800);
    TPad* p1 = new TPad("p1_" + cname, "p1", 0, 0.3, 1, 1.0);
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(0.12);
    p1->SetRightMargin(0.05);
    p1->SetTopMargin(0.05);
    p1->Draw();
    TPad* p2 = new TPad("p2_" + cname, "p2", 0, 0.0, 1, 0.3);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.30);
    p2->SetLeftMargin(0.12);
    p2->SetRightMargin(0.05);
    p2->Draw();

    p1->cd();
    gPad->SetGrid();
    TH1* hS = (TH1*)hSignal->Clone("hS_" + cname);
    hS->GetXaxis()->SetLabelSize(0);
    hS->GetXaxis()->SetTitleSize(0);
    hS->Draw("E");

    fTot->SetLineColor(kBlue);
    fTot->SetLineWidth(2);
    fTot->Draw("SAME");
    fBkg->SetLineColor(kCyan);
    fBkg->SetLineStyle(3);
    fBkg->Draw("SAME");

    // Recreate isolated Signal functions for drawing
    double parJ[7], parP[7];
    SetTiedParameters((Double_t*)fTot->GetParameters(), parJ, parP);
    TF1* fJpsi = new TF1("fJpsi_" + cname, DSCB, drawMin, drawMax, 7);
    fJpsi->SetParameters(parJ);
    fJpsi->SetLineColor(kMagenta);
    fJpsi->SetLineStyle(2);
    fJpsi->SetLineWidth(2);
    fJpsi->Draw("SAME");

    TF1* fPsi2s = new TF1("fPsi2s_" + cname, DSCB, drawMin, drawMax, 7);
    fPsi2s->SetParameters(parP);
    fPsi2s->SetLineColor(kRed);
    fPsi2s->SetLineStyle(2);
    fPsi2s->SetLineWidth(2);
    fPsi2s->Draw("SAME");

    TLegend* leg = new TLegend(0.65, 0.70, 0.95, 0.95);
    leg->SetBorderSize(1);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(hS, "Data", "lep");
    leg->AddEntry(fTot, "Total Fit", "l");
    leg->AddEntry(fJpsi, "J/#psi Signal", "l");
    leg->AddEntry(fPsi2s, "#psi(2S) Signal", "l");
    leg->AddEntry(fBkg, bkgName, "l");
    leg->Draw();

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextFont(42);
    latex->DrawLatex(0.60, 0.56, Form("#chi^{2}/ndf = %.1f / %d = %.2f", fTot->GetChisquare(), fTot->GetNDF(), fTot->GetChisquare() / fTot->GetNDF()));
    latex->DrawLatex(0.60, 0.50, Form("M_{J/#psi} = %.3f #pm %.3f", fTot->GetParameter(1), fTot->GetParError(1)));
    latex->DrawLatex(0.60, 0.44, Form("#sigma_{J/#psi} = %.3f #pm %.3f", fTot->GetParameter(2), fTot->GetParError(2)));

    p2->cd();
    gPad->SetGridy();
    TH1D* hR = (TH1D*)hS->Clone("hR_" + cname);
    hR->SetTitle("");
    for (int i = 1; i <= hR->GetNbinsX(); i++) {
      double xc = hR->GetBinCenter(i);
      if (xc >= drawMin && xc <= drawMax) {
        double d = hS->GetBinContent(i);
        double fv = fTot->Eval(xc);
        hR->SetBinContent(i, fv > 0 ? (d - fv) / fv : 0);
        hR->SetBinError(i, fv > 0 ? hS->GetBinError(i) / fv : 0);
      } else {
        hR->SetBinContent(i, 0);
        hR->SetBinError(i, 0);
      }
    }
    hR->SetLineColor(kBlack);
    hR->SetMarkerStyle(20);
    hR->GetYaxis()->SetTitle("(Data-Fit)/Fit");
    hR->GetYaxis()->SetRangeUser(-1.5, 1.5);
    hR->GetYaxis()->SetNdivisions(505);
    hR->GetYaxis()->SetLabelSize(0.10);
    hR->GetYaxis()->SetTitleSize(0.12);
    hR->GetYaxis()->SetTitleOffset(0.4);
    hR->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
    hR->GetXaxis()->SetLabelSize(0.10);
    hR->GetXaxis()->SetTitleSize(0.12);
    hR->GetXaxis()->SetTitleOffset(0.9);
    hR->Draw("EP");
    TF1* l0 = new TF1("l0_" + cname, "0", drawMin, drawMax);
    l0->SetLineColor(kRed);
    l0->SetLineStyle(2);
    l0->Draw("SAME");
    c->SaveAs(outName);
  };

  // ===========================================
  // Draw All Fits
  // ===========================================
  drawCanvas("c_nom", "Nominal (Expo3)", fitFunc, fBkg0, "Background (Expo3)", "Fit_Nominal_Expo3.png");
  drawCanvas("c_sys1", "Sys1 (ExpoPol4)", fitFunc_Sys1, fBkg1, "Background (ExpoPol4)", "Fit_Sys1_ExpoPol4.png");
  drawCanvas("c_sys2", "Sys2 (VWG)", fitFunc_Sys2, fBkg2, "Background (VWG)", "Fit_Sys2_VWG.png");
  drawCanvas("c_sys3", "Sys3 (ExpoPol2)", fitFunc_Sys3, fBkg3, "Background (ExpoPol2)", "Fit_Sys3_ExpoPol2.png");
  drawCanvas("c_sys4", "Sys4 (Jpsi Tails Fixed)", fitFunc_Sys4, fBkg4, "Background (Expo3)", "Fit_Sys4_JpsiTails.png");
  drawCanvas("c_sys5", "Sys5 (Width +0.1)", fitFunc_Sys5, fBkg5, "Background (Expo3)", "Fit_Sys5_WidthPlus.png");
  drawCanvas("c_sys6", "Sys6 (Width -0.1)", fitFunc_Sys6, fBkg6, "Background (Expo3)", "Fit_Sys6_WidthMinus.png");

  f->Close();
}
