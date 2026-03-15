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
// Dimuon Continuum: Expo3
// ------------------------------------------
Double_t Expo3(Double_t* x, Double_t* par)
{
  double m = x[0];
  return par[0] * TMath::Exp(par[1] * m + par[2] * m * m);
}

Double_t Expo3Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint(); return 0;
  }
  return par[0] * TMath::Exp(par[1] * x[0] + par[2] * x[0] * x[0]);
}

// ------------------------------------------
// Background: ExpoPol4
// ------------------------------------------
Double_t ExpoPol4(Double_t* x, Double_t* par)
{
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m + par[5] * m * m * m + par[6] * m * m * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}

Double_t ExpoPol4Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint(); return 0;
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
  if (sigma <= 0) return 0;
  return par[0] * TMath::Exp(-(m - mbar) * (m - mbar) / (2.0 * sigma * sigma));
}

// ------------------------------------------
// Background: ExpoPol2
// ------------------------------------------
Double_t ExpoPol2(Double_t* x, Double_t* par)
{
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}

Double_t ExpoPol2Sideband(Double_t* x, Double_t* par)
{
  if ((x[0] > 2.9 && x[0] < 3.3) || (x[0] > 3.55 && x[0] < 3.8)) {
    TF1::RejectPoint(); return 0;
  }
  double m = x[0];
  double poly = par[2] + par[3] * m + par[4] * m * m;
  return par[0] * TMath::Exp(par[1] * m) * poly;
}

// ------------------------------------------
// Total Fit Functions
// ------------------------------------------
Double_t TotalFit_Main(Double_t* x, Double_t* par) {
  return DSCB(x, &par[0]) + DSCB(x, &par[7]) + Expo3(x, &par[14]);
}
Double_t TotalFit_Sys_ExpoPol4(Double_t* x, Double_t* par) {
  return DSCB(x, &par[0]) + DSCB(x, &par[7]) + ExpoPol4(x, &par[14]);
}
Double_t TotalFit_Sys_VWG(Double_t* x, Double_t* par) {
  return DSCB(x, &par[0]) + DSCB(x, &par[7]) + VMG(x, &par[14]);
}
Double_t TotalFit_Sys_ExpoPol2(Double_t* x, Double_t* par) {
  return DSCB(x, &par[0]) + DSCB(x, &par[7]) + ExpoPol2(x, &par[14]);
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

  hUnlike->Rebin(5);
  hLike->Rebin(5);
  
  double drawMin = 1.8;
  double drawMax = 6.0;
  double fitRangeMin = 1.8;
  double fitRangeMax = 6.0;
  
  // ========================================
  // Signal Extraction
  // ========================================
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->SetTitle("Signal (Unlike - Like) ALL; m_{#mu#mu} [GeV/c^{2}]; Counts / (50 MeV/c^{2})");
  hSignal->Add(hLike, -1.0);
  double ymin = hSignal->GetMinimum();
  if (ymin > 0) ymin = 0;
  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);
  hSignal->GetYaxis()->SetRangeUser(ymin * 1.2, hSignal->GetMaximum() * 1.5);

  // =======================================
  // Background Sideband Fits
  // =======================================
  double initScale = hSignal->GetBinContent(hSignal->FindBin(fitRangeMin));
  if (initScale <= 0) initScale = 100.0;

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
  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);

  // =======================================
  // Nominal Fit (Expo3 - Psi2S Fixed)
  // =======================================
  TF1* fitFunc = new TF1("fitFunc", TotalFit_Main, fitRangeMin, fitRangeMax, 17);
  fitFunc->SetParameter(0, maxValJpsi);
  fitFunc->SetParameter(1, 3.075);
  fitFunc->SetParameter(2, 0.08);
  fitFunc->SetParLimits(0, 0.0, maxValJpsi * 5.0);
  fitFunc->SetParLimits(1, 3.05, 3.10);
  fitFunc->SetParLimits(2, 0.05, 0.12);
  fitFunc->FixParameter(3, 1.5);
  fitFunc->FixParameter(4, 3.0);
  fitFunc->FixParameter(5, 1.5);
  fitFunc->FixParameter(6, 3.0);
  
  fitFunc->SetParameter(7, maxValJpsi * 0.02);
  fitFunc->SetParLimits(7, 0.0, maxValJpsi);
  fitFunc->FixParameter(8, 3.664); 
  fitFunc->FixParameter(9, 0.095); 
  fitFunc->FixParameter(10, 1.5);
  fitFunc->FixParameter(11, 3.0);
  fitFunc->FixParameter(12, 1.5);
  fitFunc->FixParameter(13, 3.0);
  
  fitFunc->SetParameter(14, fExpo3Sideband->GetParameter(0));
  fitFunc->SetParameter(15, fExpo3Sideband->GetParameter(1));
  fitFunc->SetParameter(16, fExpo3Sideband->GetParameter(2));

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Nominal Fit (Expo3, Psi2S Fixed) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc", "RM0");
  
  TF1* fitSigJpsi_Main = new TF1("fitSigJpsi_Main", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Main->SetParameter(i, fitFunc->GetParameter(i));
  double yieldJpsi_Main = fitSigJpsi_Main->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Main   = yieldJpsi_Main * (fitFunc->GetParError(0) / fitFunc->GetParameter(0));

  TF1* fitSigPsi2s_Main = new TF1("fitSigPsi2s_Main", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Main->SetParameter(i, fitFunc->GetParameter(7 + i));
  double yieldPsi2s_Main = fitSigPsi2s_Main->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Main   = yieldPsi2s_Main * (fitFunc->GetParError(7) / fitFunc->GetParameter(7));

  // =======================================
  // Sys1 Fit (ExpoPol4)
  // =======================================
  TF1* fitFunc_Sys1 = new TF1("fitFunc_Sys1", TotalFit_Sys_ExpoPol4, fitRangeMin, fitRangeMax, 21);
  for (int i = 0; i < 14; i++) {
    fitFunc_Sys1->SetParameter(i, fitFunc->GetParameter(i));
    if (i == 3 || i == 4 || i == 5 || i == 6 || i == 8 || i == 9 || i == 10 || i == 11 || i == 12 || i == 13) {
      fitFunc_Sys1->FixParameter(i, fitFunc->GetParameter(i));
    } else if (i == 0 || i == 1 || i == 2 || i == 7) {
      double pmin, pmax; fitFunc->GetParLimits(i, pmin, pmax);
      if (pmin < pmax) fitFunc_Sys1->SetParLimits(i, pmin, pmax);
    }
  }
  for (int i = 0; i < 7; i++) fitFunc_Sys1->SetParameter(14 + i, fExpoPol4Sideband->GetParameter(i));
  fitFunc_Sys1->FixParameter(14 + 2, 1.0);

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 1 (ExpoPol4) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys1", "RM0");
  
  TF1* fitSigJpsi_Sys1 = new TF1("fitSigJpsi_Sys1", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Sys1->SetParameter(i, fitFunc_Sys1->GetParameter(i));
  double yieldJpsi_Sys1 = fitSigJpsi_Sys1->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Sys1   = yieldJpsi_Sys1 * (fitFunc_Sys1->GetParError(0) / fitFunc_Sys1->GetParameter(0));

  TF1* fitSigPsi2s_Sys1 = new TF1("fitSigPsi2s_Sys1", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Sys1->SetParameter(i, fitFunc_Sys1->GetParameter(7 + i));
  double yieldPsi2s_Sys1 = fitSigPsi2s_Sys1->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Sys1   = yieldPsi2s_Sys1 * (fitFunc_Sys1->GetParError(7) / fitFunc_Sys1->GetParameter(7));

  // =======================================
  // Sys2 Fit (VMG)
  // =======================================
  TF1* fitFunc_Sys2 = new TF1("fitFunc_Sys2", TotalFit_Sys_VWG, fitRangeMin, fitRangeMax, 18);
  for (int i = 0; i < 14; i++) {
    fitFunc_Sys2->SetParameter(i, fitFunc->GetParameter(i));
    if (i == 3 || i == 4 || i == 5 || i == 6 || i == 8 || i == 9 || i == 10 || i == 11 || i == 12 || i == 13) {
      fitFunc_Sys2->FixParameter(i, fitFunc->GetParameter(i));
    } else if (i == 0 || i == 1 || i == 2 || i == 7) {
      double pmin, pmax; fitFunc->GetParLimits(i, pmin, pmax);
      if (pmin < pmax) fitFunc_Sys2->SetParLimits(i, pmin, pmax);
    }
  }
  fitFunc_Sys2->SetParameter(14, 100.0);
  fitFunc_Sys2->SetParameter(15, 3.1);
  fitFunc_Sys2->SetParameter(16, 1.0);
  fitFunc_Sys2->SetParameter(17, 0.0);

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 2 (VMG) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys2", "RM0");
  
  TF1* fitSigJpsi_Sys2 = new TF1("fitSigJpsi_Sys2", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Sys2->SetParameter(i, fitFunc_Sys2->GetParameter(i));
  double yieldJpsi_Sys2 = fitSigJpsi_Sys2->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Sys2   = yieldJpsi_Sys2 * (fitFunc_Sys2->GetParError(0) / fitFunc_Sys2->GetParameter(0));

  TF1* fitSigPsi2s_Sys2 = new TF1("fitSigPsi2s_Sys2", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Sys2->SetParameter(i, fitFunc_Sys2->GetParameter(7 + i));
  double yieldPsi2s_Sys2 = fitSigPsi2s_Sys2->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Sys2   = yieldPsi2s_Sys2 * (fitFunc_Sys2->GetParError(7) / fitFunc_Sys2->GetParameter(7));

  // =======================================
  // Sys3 Fit (ExpoPol2)
  // =======================================
  TF1* fitFunc_Sys3 = new TF1("fitFunc_Sys3", TotalFit_Sys_ExpoPol2, fitRangeMin, fitRangeMax, 19);
  for (int i = 0; i < 14; i++) {
    fitFunc_Sys3->SetParameter(i, fitFunc->GetParameter(i));
    if (i == 3 || i == 4 || i == 5 || i == 6 || i == 8 || i == 9 || i == 10 || i == 11 || i == 12 || i == 13) {
      fitFunc_Sys3->FixParameter(i, fitFunc->GetParameter(i));
    } else if (i == 0 || i == 1 || i == 2 || i == 7) {
      double pmin, pmax; fitFunc->GetParLimits(i, pmin, pmax);
      if (pmin < pmax) fitFunc_Sys3->SetParLimits(i, pmin, pmax);
    }
  }
  for (int i = 0; i < 5; i++) fitFunc_Sys3->SetParameter(14 + i, fExpoPol2Sideband->GetParameter(i));
  fitFunc_Sys3->FixParameter(14 + 2, 1.0);

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 3 (ExpoPol2) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys3", "RM0");
  
  TF1* fitSigJpsi_Sys3 = new TF1("fitSigJpsi_Sys3", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Sys3->SetParameter(i, fitFunc_Sys3->GetParameter(i));
  double yieldJpsi_Sys3 = fitSigJpsi_Sys3->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Sys3   = yieldJpsi_Sys3 * (fitFunc_Sys3->GetParError(0) / fitFunc_Sys3->GetParameter(0));

  TF1* fitSigPsi2s_Sys3 = new TF1("fitSigPsi2s_Sys3", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Sys3->SetParameter(i, fitFunc_Sys3->GetParameter(7 + i));
  double yieldPsi2s_Sys3 = fitSigPsi2s_Sys3->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Sys3   = yieldPsi2s_Sys3 * (fitFunc_Sys3->GetParError(7) / fitFunc_Sys3->GetParameter(7));

  // =======================================
  // Sys4 Fit (Expo3, Free Psi2S Mass/Sigma)
  // =======================================
  TF1* fitFunc_Sys4 = new TF1("fitFunc_Sys4", TotalFit_Main, fitRangeMin, fitRangeMax, 17);
  for (int i = 0; i < 17; i++) {
    fitFunc_Sys4->SetParameter(i, fitFunc->GetParameter(i));
    if (i == 3 || i == 4 || i == 5 || i == 6 || i == 10 || i == 11 || i == 12 || i == 13) {
      fitFunc_Sys4->FixParameter(i, fitFunc->GetParameter(i)); // テールパラメータのみ固定
    } else if (i == 8) {
      fitFunc_Sys4->SetParLimits(i, 3.60, 3.75); // psi(2S) 質量をフリーに
    } else if (i == 9) {
      fitFunc_Sys4->SetParLimits(i, 0.05, 0.15); // psi(2S) シグマをフリーに
    } else if (i == 0 || i == 1 || i == 2 || i == 7) {
      double pmin, pmax; fitFunc->GetParLimits(i, pmin, pmax);
      if (pmin < pmax) fitFunc_Sys4->SetParLimits(i, pmin, pmax);
    }
  }

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 4 (Expo3, Free Psi2S) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys4", "RM0");
  
  TF1* fitSigJpsi_Sys4 = new TF1("fitSigJpsi_Sys4", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Sys4->SetParameter(i, fitFunc_Sys4->GetParameter(i));
  double yieldJpsi_Sys4 = fitSigJpsi_Sys4->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Sys4   = yieldJpsi_Sys4 * (fitFunc_Sys4->GetParError(0) / fitFunc_Sys4->GetParameter(0));

  TF1* fitSigPsi2s_Sys4 = new TF1("fitSigPsi2s_Sys4", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Sys4->SetParameter(i, fitFunc_Sys4->GetParameter(7 + i));
  double yieldPsi2s_Sys4 = fitSigPsi2s_Sys4->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Sys4   = yieldPsi2s_Sys4 * (fitFunc_Sys4->GetParError(7) / fitFunc_Sys4->GetParameter(7));

  // =======================================
  // Sys5 Fit (Expo3, Tail Parameters Low)
  // alpha = 1.2, n = 2.5
  // =======================================
  TF1* fitFunc_Sys5 = new TF1("fitFunc_Sys5", TotalFit_Main, fitRangeMin, fitRangeMax, 17);
  for (int i = 0; i < 17; i++) {
    fitFunc_Sys5->SetParameter(i, fitFunc->GetParameter(i));
    if (i == 8 || i == 9) {
      fitFunc_Sys5->FixParameter(i, fitFunc->GetParameter(i)); // Psi2Sの質量とシグマは固定
    } else if (i == 0 || i == 1 || i == 2 || i == 7) {
      double pmin, pmax; fitFunc->GetParLimits(i, pmin, pmax);
      if (pmin < pmax) fitFunc_Sys5->SetParLimits(i, pmin, pmax);
    }
  }
  fitFunc_Sys5->FixParameter(3, 1.2); fitFunc_Sys5->FixParameter(4, 2.5);
  fitFunc_Sys5->FixParameter(5, 1.2); fitFunc_Sys5->FixParameter(6, 2.5);
  fitFunc_Sys5->FixParameter(10, 1.2); fitFunc_Sys5->FixParameter(11, 2.5);
  fitFunc_Sys5->FixParameter(12, 1.2); fitFunc_Sys5->FixParameter(13, 2.5);

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 5 (Expo3, Tail Low) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys5", "RM0");
  
  TF1* fitSigJpsi_Sys5 = new TF1("fitSigJpsi_Sys5", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Sys5->SetParameter(i, fitFunc_Sys5->GetParameter(i));
  double yieldJpsi_Sys5 = fitSigJpsi_Sys5->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Sys5   = yieldJpsi_Sys5 * (fitFunc_Sys5->GetParError(0) / fitFunc_Sys5->GetParameter(0));

  TF1* fitSigPsi2s_Sys5 = new TF1("fitSigPsi2s_Sys5", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Sys5->SetParameter(i, fitFunc_Sys5->GetParameter(7 + i));
  double yieldPsi2s_Sys5 = fitSigPsi2s_Sys5->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Sys5   = yieldPsi2s_Sys5 * (fitFunc_Sys5->GetParError(7) / fitFunc_Sys5->GetParameter(7));

  // =======================================
  // Sys6 Fit (Expo3, Tail Parameters High)
  // alpha = 1.8, n = 3.5
  // =======================================
  TF1* fitFunc_Sys6 = new TF1("fitFunc_Sys6", TotalFit_Main, fitRangeMin, fitRangeMax, 17);
  for (int i = 0; i < 17; i++) {
    fitFunc_Sys6->SetParameter(i, fitFunc->GetParameter(i));
    if (i == 8 || i == 9) {
      fitFunc_Sys6->FixParameter(i, fitFunc->GetParameter(i)); // Psi2Sの質量とシグマは固定
    } else if (i == 0 || i == 1 || i == 2 || i == 7) {
      double pmin, pmax; fitFunc->GetParLimits(i, pmin, pmax);
      if (pmin < pmax) fitFunc_Sys6->SetParLimits(i, pmin, pmax);
    }
  }
  fitFunc_Sys6->FixParameter(3, 1.8); fitFunc_Sys6->FixParameter(4, 3.5);
  fitFunc_Sys6->FixParameter(5, 1.8); fitFunc_Sys6->FixParameter(6, 3.5);
  fitFunc_Sys6->FixParameter(10, 1.8); fitFunc_Sys6->FixParameter(11, 3.5);
  fitFunc_Sys6->FixParameter(12, 1.8); fitFunc_Sys6->FixParameter(13, 3.5);

  std::cout << "\n==========================================" << std::endl;
  std::cout << "=== Systematic Fit 6 (Expo3, Tail High) ===" << std::endl;
  std::cout << "==========================================" << std::endl;
  hSignal->Fit("fitFunc_Sys6", "RM0");
  
  TF1* fitSigJpsi_Sys6 = new TF1("fitSigJpsi_Sys6", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigJpsi_Sys6->SetParameter(i, fitFunc_Sys6->GetParameter(i));
  double yieldJpsi_Sys6 = fitSigJpsi_Sys6->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errJpsi_Sys6   = yieldJpsi_Sys6 * (fitFunc_Sys6->GetParError(0) / fitFunc_Sys6->GetParameter(0));

  TF1* fitSigPsi2s_Sys6 = new TF1("fitSigPsi2s_Sys6", DSCB, fitRangeMin, fitRangeMax, 7);
  for (int i = 0; i < 7; i++) fitSigPsi2s_Sys6->SetParameter(i, fitFunc_Sys6->GetParameter(7 + i));
  double yieldPsi2s_Sys6 = fitSigPsi2s_Sys6->Integral(fitRangeMin, fitRangeMax) / binWidth;
  double errPsi2s_Sys6   = yieldPsi2s_Sys6 * (fitFunc_Sys6->GetParError(7) / fitFunc_Sys6->GetParameter(7));

  // =======================================
  // Yield Summary Output (WITH ERRORS)
  // =======================================
  std::cout << "\n------------- YIELD SUMMARY --------------\n";
  std::cout << "--- J/psi Yield ---" << std::endl;
  std::cout << Form(" Nominal (Expo3)     = %.1f +/- %.1f\n", yieldJpsi_Main, errJpsi_Main);
  std::cout << Form(" Sys1 (ExpoPol4)     = %.1f +/- %.1f  <-- *NOT POSDEF*\n", yieldJpsi_Sys1, errJpsi_Sys1);
  std::cout << Form(" Sys2 (VMG)          = %.1f +/- %.1f\n", yieldJpsi_Sys2, errJpsi_Sys2);
  std::cout << Form(" Sys3 (ExpoPol2)     = %.1f +/- %.1f\n", yieldJpsi_Sys3, errJpsi_Sys3);
  std::cout << Form(" Sys4 (Free Psi2S)   = %.1f +/- %.1f\n", yieldJpsi_Sys4, errJpsi_Sys4);
  std::cout << Form(" Sys5 (Tail Low)     = %.1f +/- %.1f\n", yieldJpsi_Sys5, errJpsi_Sys5);
  std::cout << Form(" Sys6 (Tail High)    = %.1f +/- %.1f\n", yieldJpsi_Sys6, errJpsi_Sys6);
  
  std::cout << "\n--- psi(2S) Yield ---" << std::endl;
  std::cout << Form(" Nominal (Expo3)     = %.1f +/- %.1f\n", yieldPsi2s_Main, errPsi2s_Main);
  std::cout << Form(" Sys1 (ExpoPol4)     = %.1f +/- %.1f  <-- *NOT POSDEF*\n", yieldPsi2s_Sys1, errPsi2s_Sys1);
  std::cout << Form(" Sys2 (VMG)          = %.1f +/- %.1f\n", yieldPsi2s_Sys2, errPsi2s_Sys2);
  std::cout << Form(" Sys3 (ExpoPol2)     = %.1f +/- %.1f\n", yieldPsi2s_Sys3, errPsi2s_Sys3);
  std::cout << Form(" Sys4 (Free Psi2S)   = %.1f +/- %.1f\n", yieldPsi2s_Sys4, errPsi2s_Sys4);
  std::cout << Form(" Sys5 (Tail Low)     = %.1f +/- %.1f\n", yieldPsi2s_Sys5, errPsi2s_Sys5);
  std::cout << Form(" Sys6 (Tail High)    = %.1f +/- %.1f\n", yieldPsi2s_Sys6, errPsi2s_Sys6);
  std::cout << "------------------------------------------\n\n";

  // ===========================================
  // Draw Macros Helper Function (Chi2/ndf直接取得版)
  // ===========================================
  // 引数に origChi2 と origNdf を追加し、ダミーのFit処理を削除しました。
  auto drawCanvas = [&](TString cname, TString title, TF1* fTot, TF1* fJpsi, TF1* fPsi2s, TF1* fBkg, double yJ, double eJ, double yP, double eP, double origChi2, int origNdf, TString bkgName, TString outName) {
    TCanvas* c = new TCanvas(cname, title, 800, 800);
    TPad* p1 = new TPad("p1_"+cname, "p1", 0, 0.3, 1, 1.0);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.12); p1->SetRightMargin(0.05); p1->SetTopMargin(0.05); p1->Draw();
    TPad* p2 = new TPad("p2_"+cname, "p2", 0, 0.0, 1, 0.3);
    p2->SetTopMargin(0.02); p2->SetBottomMargin(0.30); p2->SetLeftMargin(0.12); p2->SetRightMargin(0.05); p2->Draw();

    p1->cd(); gPad->SetGrid();
    TH1* hS = (TH1*)hSignal->Clone("hS_"+cname);
    hS->GetXaxis()->SetLabelSize(0); hS->GetXaxis()->SetTitleSize(0);
    hS->Draw("E");

    fTot->SetLineColor(kBlue); fTot->SetLineWidth(2); fTot->Draw("SAME");
    fBkg->SetLineColor(kCyan); fBkg->SetLineStyle(3); fBkg->Draw("SAME");
    fJpsi->SetLineColor(kMagenta); fJpsi->SetLineStyle(2); fJpsi->SetLineWidth(2); fJpsi->Draw("SAME");
    fPsi2s->SetLineColor(kRed); fPsi2s->SetLineStyle(2); fPsi2s->SetLineWidth(2); fPsi2s->Draw("SAME");

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

    double chi2ndf = (origNdf > 0) ? origChi2 / origNdf : 0.0;

    TLatex* latex = new TLatex(); latex->SetNDC(); latex->SetTextSize(0.04); latex->SetTextColor(kBlack);
    latex->SetTextFont(42);
    latex->DrawLatex(0.60, 0.56, Form("p^{#mu#mu}_{T} < 0.25 GeV/c"));
    latex->DrawLatex(0.60, 0.51, Form("-4.00 < y < 4.00"));
    
    latex->DrawLatex(0.60, 0.46, Form("#chi^{2}/ndf = %.1f / %d = %.2f", origChi2, origNdf, chi2ndf));

    latex->DrawLatex(0.60, 0.41, Form("M_{J/#psi} = %.3f #pm %.3f GeV/c^{2}", fTot->GetParameter(1), fTot->GetParError(1)));
    latex->DrawLatex(0.60, 0.36, Form("#sigma_{J/#psi} = %.3f #pm %.3f GeV/c^{2}", fTot->GetParameter(2), fTot->GetParError(2)));
    
    if (fTot->GetParError(8) > 0) {
      latex->DrawLatex(0.60, 0.31, Form("M_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}", fTot->GetParameter(8), fTot->GetParError(8)));
    } else {
      latex->DrawLatex(0.60, 0.31, Form("M_{#psi(2S)} = %.3f GeV/c^{2} (Fixed)", fTot->GetParameter(8)));
    }
    if (fTot->GetParError(9) > 0) {
      latex->DrawLatex(0.60, 0.26, Form("#sigma_{#psi(2S)} = %.3f #pm %.3f GeV/c^{2}", fTot->GetParameter(9), fTot->GetParError(9)));
    } else {
      latex->DrawLatex(0.60, 0.26, Form("#sigma_{#psi(2S)} = %.3f GeV/c^{2} (Fixed)", fTot->GetParameter(9)));
    }
    
    latex->DrawLatex(0.60, 0.21, Form("N_{J/#psi} = %.0f #pm %.0f", yJ, eJ));
    latex->DrawLatex(0.60, 0.16, Form("N_{#psi(2S)} = %.0f #pm %.0f", yP, eP));
    latex->DrawLatex(0.60, 0.11, Form("N_{#psi(2S)}/N_{J/#psi} = %.3f #pm %.3f", yP/yJ, yP/yJ * TMath::Sqrt((eP/yP)*(eP/yP) + (eJ/yJ)*(eJ/yJ))));

    p2->cd(); gPad->SetGridy();
    TH1D* hR = (TH1D*)hS->Clone("hR_"+cname); hR->SetTitle("");
    for (int i = 1; i <= hR->GetNbinsX(); i++) {
      double xc = hR->GetBinCenter(i);
      if (xc >= drawMin && xc <= drawMax) {
        double d = hS->GetBinContent(i);
        double fv = fTot->Eval(xc);
        if (fv > 0) { hR->SetBinContent(i, (d - fv)/fv); hR->SetBinError(i, hS->GetBinError(i)/fv); }
        else { hR->SetBinContent(i, 0); hR->SetBinError(i, 0); }
      } else { hR->SetBinContent(i, 0); hR->SetBinError(i, 0); }
    }
    hR->SetLineColor(kBlack); hR->SetMarkerColor(kBlack); hR->SetMarkerStyle(20);
    hR->GetYaxis()->SetTitle("(Data-Fit)/Fit"); hR->GetYaxis()->SetRangeUser(-1.5, 1.5);
    hR->GetYaxis()->SetNdivisions(505); hR->GetYaxis()->SetLabelSize(0.10); hR->GetYaxis()->SetTitleSize(0.12); hR->GetYaxis()->SetTitleOffset(0.4);
    hR->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]"); hR->GetXaxis()->SetLabelSize(0.10); hR->GetXaxis()->SetTitleSize(0.12); hR->GetXaxis()->SetTitleOffset(0.9);
    hR->Draw("EP");
    TF1* l0 = new TF1("l0_"+cname, "0", drawMin, drawMax); l0->SetLineColor(kRed); l0->SetLineStyle(2); l0->Draw("SAME");
    c->SaveAs(outName);
  };

  // ===========================================
  // Draw All Fits (直接パラメータと誤差、Chi2/ndfを渡す)
  // ===========================================
  
  // Nominal
  TF1* fTot0 = new TF1("fTot0", TotalFit_Main, drawMin, drawMax, 17); 
  for(int i=0;i<17;i++) { fTot0->SetParameter(i, fitFunc->GetParameter(i)); fTot0->SetParError(i, fitFunc->GetParError(i)); }
  TF1* fBkg0 = new TF1("fBkg0", Expo3, drawMin, drawMax, 3); for(int i=0;i<3;i++) fBkg0->SetParameter(i, fitFunc->GetParameter(14+i));
  drawCanvas("c2", "Nominal (Expo3)", fTot0, fitSigJpsi_Main, fitSigPsi2s_Main, fBkg0, yieldJpsi_Main, errJpsi_Main, yieldPsi2s_Main, errPsi2s_Main, fitFunc->GetChisquare(), fitFunc->GetNDF(), "Background (Expo3)", "Unlike_bkg_WideRange_crystal.png");

  // Sys1 (ExpoPol4)
  TF1* fTot1 = new TF1("fTot1", TotalFit_Sys_ExpoPol4, drawMin, drawMax, 21); 
  for(int i=0;i<21;i++) { fTot1->SetParameter(i, fitFunc_Sys1->GetParameter(i)); fTot1->SetParError(i, fitFunc_Sys1->GetParError(i)); }
  TF1* fBkg1 = new TF1("fBkg1", ExpoPol4, drawMin, drawMax, 7); for(int i=0;i<7;i++) fBkg1->SetParameter(i, fitFunc_Sys1->GetParameter(14+i));
  drawCanvas("c_sys1", "Sys1 (ExpoPol4)", fTot1, fitSigJpsi_Sys1, fitSigPsi2s_Sys1, fBkg1, yieldJpsi_Sys1, errJpsi_Sys1, yieldPsi2s_Sys1, errPsi2s_Sys1, fitFunc_Sys1->GetChisquare(), fitFunc_Sys1->GetNDF(), "Background (ExpoPol4)", "Unlike_bkg_WideRange_crystal_ExpoPol4.png");

  // Sys2 (VMG)
  TF1* fTot2 = new TF1("fTot2", TotalFit_Sys_VWG, drawMin, drawMax, 18); 
  for(int i=0;i<18;i++) { fTot2->SetParameter(i, fitFunc_Sys2->GetParameter(i)); fTot2->SetParError(i, fitFunc_Sys2->GetParError(i)); }
  TF1* fBkg2 = new TF1("fBkg2", VMG, drawMin, drawMax, 4); for(int i=0;i<4;i++) fBkg2->SetParameter(i, fitFunc_Sys2->GetParameter(14+i));
  drawCanvas("c_sys2", "Sys2 (VMG)", fTot2, fitSigJpsi_Sys2, fitSigPsi2s_Sys2, fBkg2, yieldJpsi_Sys2, errJpsi_Sys2, yieldPsi2s_Sys2, errPsi2s_Sys2, fitFunc_Sys2->GetChisquare(), fitFunc_Sys2->GetNDF(), "Background (VMG)", "Unlike_bkg_WideRange_crystal_VMG.png");

  // Sys3 (ExpoPol2)
  TF1* fTot3 = new TF1("fTot3", TotalFit_Sys_ExpoPol2, drawMin, drawMax, 19); 
  for(int i=0;i<19;i++) { fTot3->SetParameter(i, fitFunc_Sys3->GetParameter(i)); fTot3->SetParError(i, fitFunc_Sys3->GetParError(i)); }
  TF1* fBkg3 = new TF1("fBkg3", ExpoPol2, drawMin, drawMax, 5); for(int i=0;i<5;i++) fBkg3->SetParameter(i, fitFunc_Sys3->GetParameter(14+i));
  drawCanvas("c_sys3", "Sys3 (ExpoPol2)", fTot3, fitSigJpsi_Sys3, fitSigPsi2s_Sys3, fBkg3, yieldJpsi_Sys3, errJpsi_Sys3, yieldPsi2s_Sys3, errPsi2s_Sys3, fitFunc_Sys3->GetChisquare(), fitFunc_Sys3->GetNDF(), "Background (ExpoPol2)", "Unlike_bkg_WideRange_crystal_ExpoPol2.png");

  // Sys4 (Free Psi2S)
  TF1* fTot4 = new TF1("fTot4", TotalFit_Main, drawMin, drawMax, 17); 
  for(int i=0;i<17;i++) { fTot4->SetParameter(i, fitFunc_Sys4->GetParameter(i)); fTot4->SetParError(i, fitFunc_Sys4->GetParError(i)); }
  TF1* fBkg4 = new TF1("fBkg4", Expo3, drawMin, drawMax, 3); for(int i=0;i<3;i++) fBkg4->SetParameter(i, fitFunc_Sys4->GetParameter(14+i));
  drawCanvas("c_sys4", "Sys4 (Free Psi2S)", fTot4, fitSigJpsi_Sys4, fitSigPsi2s_Sys4, fBkg4, yieldJpsi_Sys4, errJpsi_Sys4, yieldPsi2s_Sys4, errPsi2s_Sys4, fitFunc_Sys4->GetChisquare(), fitFunc_Sys4->GetNDF(), "Background (Expo3)", "Unlike_bkg_WideRange_crystal_FreePsi2S.png");

  // Sys5 (Tail Low)
  TF1* fTot5 = new TF1("fTot5", TotalFit_Main, drawMin, drawMax, 17); 
  for(int i=0;i<17;i++) { fTot5->SetParameter(i, fitFunc_Sys5->GetParameter(i)); fTot5->SetParError(i, fitFunc_Sys5->GetParError(i)); }
  TF1* fBkg5 = new TF1("fBkg5", Expo3, drawMin, drawMax, 3); for(int i=0;i<3;i++) fBkg5->SetParameter(i, fitFunc_Sys5->GetParameter(14+i));
  drawCanvas("c_sys5", "Sys5 (Tail Low)", fTot5, fitSigJpsi_Sys5, fitSigPsi2s_Sys5, fBkg5, yieldJpsi_Sys5, errJpsi_Sys5, yieldPsi2s_Sys5, errPsi2s_Sys5, fitFunc_Sys5->GetChisquare(), fitFunc_Sys5->GetNDF(), "Background (Expo3)", "Unlike_bkg_WideRange_crystal_TailLow.png");

  // Sys6 (Tail High)
  TF1* fTot6 = new TF1("fTot6", TotalFit_Main, drawMin, drawMax, 17); 
  for(int i=0;i<17;i++) { fTot6->SetParameter(i, fitFunc_Sys6->GetParameter(i)); fTot6->SetParError(i, fitFunc_Sys6->GetParError(i)); }
  TF1* fBkg6 = new TF1("fBkg6", Expo3, drawMin, drawMax, 3); for(int i=0;i<3;i++) fBkg6->SetParameter(i, fitFunc_Sys6->GetParameter(14+i));
  drawCanvas("c_sys6", "Sys6 (Tail High)", fTot6, fitSigJpsi_Sys6, fitSigPsi2s_Sys6, fBkg6, yieldJpsi_Sys6, errJpsi_Sys6, yieldPsi2s_Sys6, errPsi2s_Sys6, fitFunc_Sys6->GetChisquare(), fitFunc_Sys6->GetNDF(), "Background (Expo3)", "Unlike_bkg_WideRange_crystal_TailHigh.png");

  f->Close();
}