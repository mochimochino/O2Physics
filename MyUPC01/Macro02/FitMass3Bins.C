#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>

// ======================================================================
// Parameter Settings
// ======================================================================
struct FitConfig {
  // Mass from PDG
  double m_Jpsi_PDG = 3.0969;  // GeV/c^2
  double m_Psi2S_PDG = 3.6861; // GeV/c^2

  // Tail parameters from MC
  double mc_Jpsi_alpha2 = 2.0; // Right tail alpha
  double mc_Jpsi_n2 = 3.0;     // Right tail n

  double mc_Psi2S_alpha1 = 1.5; // Left tail alpha
  double mc_Psi2S_n1 = 2.0;     // Left tail n
  double mc_Psi2S_alpha2 = 2.0; // Right tail alpha
  double mc_Psi2S_n2 = 3.0;     // Right tail n

  // Resolution Parameters
  double mc_sigma_ratio = 1.05; // (sigma_Psi2S / sigma_Jpsi)_MC
  double sigma_factor = 1.1;

  // Fit range
  double fitMin = 2.5;
  double fitMax = 4.5;

  // Background function type (1: Expo3, 2: ExpoPol4, 3: VWG)
  // Note: ExpoPol4 and VWG are not implemented yet.
  int bgType = 1;
};

FitConfig cfg;

// ======================================================================
// Fit Functions
// ======================================================================

// Double-sided Crystal Ball Function
double DSCB(double* x, double* p)
{
  double m = x[0];
  double N = p[0];
  double m0 = p[1];
  double sigma = p[2];
  double alpha1 = p[3]; // Left tail alpha
  double n1 = p[4];     // Left tail n
  double alpha2 = p[5]; // Right tail alpha
  double n2 = p[6];     // Right tail n

  double alpha = (m - m0) / sigma;

  if (alpha < -alpha1) {
    double A = TMath::Power(n1 / TMath::Abs(alpha1), n1) * TMath::Exp(-0.5 * alpha1 * alpha1);
    double B = n1 / TMath::Abs(alpha1) - TMath::Abs(alpha1);
    return N * A * TMath::Power(B - alpha, -n1);
  } else if (alpha > alpha2) {
    double C = TMath::Power(n2 / TMath::Abs(alpha2), n2) * TMath::Exp(-0.5 * alpha2 * alpha2);
    double D = n2 / TMath::Abs(alpha2) - TMath::Abs(alpha2);
    return N * C * TMath::Power(D + alpha, -n2);
  } else {
    return N * TMath::Exp(-0.5 * alpha * alpha);
  }
}

// Expo3 Functions
double Expo3(double* x, double* p)
{
  double m = x[0];
  return TMath::Exp(p[0] + p[1] * m + p[2] * m * m);
}

// Variable Width Gaussian (VWG) Functions
double VWG(double* x, double* p)
{
  double m = x[0];
  double N = p[0];
  double m_bar = p[1];
  double A = p[2];
  double B = p[3];

  double sigma = A + B * (m - m_bar) / m_bar;
  if (sigma <= 0)
    return 0;

  return N * TMath::Exp(-0.5 * TMath::Power((m - m_bar) / sigma, 2));
}

// -----------------------
// ExpoPol4 Functions
// To be added
// -----------------------

// J/psi + Psi(2S) Signal Functions
double SignalJpsiPsi2s(double* x, double* p)
{
  // J/psi parameters
  double p_Jpsi[7];
  for (int i = 0; i < 7; i++)
    p_Jpsi[i] = p[i];

  // Psi(2S) parameters
  double p_Psi2S[7];
  p_Psi2S[0] = p[7];                                              // Psi(2S) Yield N
  p_Psi2S[1] = p_Jpsi[1] + (cfg.m_Psi2S_PDG - cfg.m_Jpsi_PDG);    // m0
  p_Psi2S[2] = p_Jpsi[2] * cfg.sigma_factor * cfg.mc_sigma_ratio; // sigma
  p_Psi2S[3] = cfg.mc_Psi2S_alpha1;
  p_Psi2S[4] = cfg.mc_Psi2S_n1;
  p_Psi2S[5] = cfg.mc_Psi2S_alpha2;
  p_Psi2S[6] = cfg.mc_Psi2S_n2;

  return DSCB(x, p_Jpsi) + DSCB(x, p_Psi2S);
}

// Total Fit FUnction: J/psi + Psi(2S) + Background
double TotalFit(double* x, double* p)
{
  double signal = SignalJpsiPsi2s(x, p);
  double bg = 0;

  if (cfg.bgType == 1) { // Expo3 (p[8], p[9], p[10])
    bg = Expo3(x, &p[8]);
  }
  // ExpoPol4
  // if (cfg.bgType == 2) { // ExpoPol4 (p[8], p[9], p[10], p[11], p[12])
  // bg = ExpoPol4(x, &p[8]);
  // }
  // VWG
  // if (cfg.bgType == 3) { // VWG (p[8], p[9], p[10], p[11])
  // bg = VWG(x, &p[8]);
  // }
  return signal + bg;
}

// ======================================================================
// Main macro
// ======================================================================
void FitMass3Bins(const char* filename = "AnalysisResults.root")
{
  gStyle->SetOptStat(0);
  int fontCode = 42;
  gStyle->SetLabelFont(fontCode, "XYZ");
  gStyle->SetTitleFont(fontCode, "XYZ");
  gStyle->SetTextFont(fontCode);
  gStyle->SetLegendFont(fontCode);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Error: " << filename << " can't open" << std::endl;
    return;
  }
  TH2D* h2 = (TH2D*)f->Get("my-upc-mass-02/registry/hMassVsRapidityUnlike");

  TCanvas* c1 = new TCanvas("c1", "Fit Invariant Mass", 800, 800);
  TPad* pad1 = new TPad("pad1", "pad1", 0.00, 0.50, 0.50, 1.00);
  TPad* pad2 = new TPad("pad2", "pad2", 0.50, 0.50, 1.00, 1.00);
  TPad* pad3 = new TPad("pad3", "pad3", 0.25, 0.00, 0.75, 0.50);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  TPad* pads[3] = {pad1, pad2, pad3};

  double y_bins[4] = {-4.00, -3.50, -3.00, -2.50};

  for (int i = 0; i < 3; ++i) {
    pads[i]->cd();
    gPad->SetMargin(0.18, 0.04, 0.15, 0.07);

    double ymin = y_bins[i];
    double ymax = y_bins[i + 1];
    int rebin = 25;

    int binMin = h2->GetXaxis()->FindBin(ymin + 1e-4);
    int binMax = h2->GetXaxis()->FindBin(ymax - 1e-4);
    TH1D* hUnlike = h2->ProjectionY(Form("hMass_bin%d", i), binMin, binMax);
    hUnlike->Rebin(rebin);
    hUnlike->SetTitle("");

    hUnlike->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
    hUnlike->GetYaxis()->SetTitle(Form("Counts / (%d MeV/#it{c}^{2})", rebin));
    hUnlike->GetXaxis()->SetRangeUser(2.0, 5.0);
    hUnlike->GetYaxis()->SetRangeUser(0, hUnlike->GetMaximum() * 1.5);
    hUnlike->SetMarkerStyle(kFullCircle);
    hUnlike->SetMarkerSize(0.5);
    hUnlike->GetYaxis()->SetTitleOffset(1.8);
    hUnlike->Draw("PE");

    // ==========================================================
    // Fitting Setting and Execution
    // ==========================================================
    // Number of parameters: J/psi(7) + Psi2S_N(1) + BG_Expo3(3) = 11
    TF1* fTotal = new TF1(Form("fTotal_bin%d", i), TotalFit, cfg.fitMin, cfg.fitMax, 11);
    fTotal->SetNpx(1000);
    fTotal->SetLineColor(kBlue);
    fTotal->SetLineWidth(1);

    // Initial values
    fTotal->SetParameter(0, hUnlike->GetMaximum()); // N_Jpsi Free parameter
    fTotal->SetParameter(1, cfg.m_Jpsi_PDG);        // m0_Jpsi Free parameter
    fTotal->SetParameter(2, 0.07);                  // sigma_Jpsi Free parameter

    /*
    fTotal->SetParLimits(0, 0.0, hUnlike->GetMaximum() * 2.0);
    fTotal->SetParLimits(1, cfg.m_Jpsi_PDG - 0.1, cfg.m_Jpsi_PDG + 0.1);
    fTotal->SetParLimits(2, 0.01, 0.20);
    */

    // J/psi Left tails (Free parameters)
    fTotal->SetParameter(3, 1.0); // alpha1 Free parameter
    fTotal->SetParameter(4, 5.0); // n1 Free parameter

    /*
    fTotal->SetParLimits(3, 0.0, 5.0);
    fTotal->SetParLimits(4, 0.0, 50.0);
    */

    // J/psi Right tails (Fix from MC)
    fTotal->FixParameter(5, cfg.mc_Jpsi_alpha2); // alpha2 Fix parameter
    fTotal->FixParameter(6, cfg.mc_Jpsi_n2);     // n2 Fix parameter

    /*
    fTotal->SetParLimits(5, 0.0, 9999.0);
    fTotal->SetParLimits(6, 0.0, 100.0);
    */

    // Psi(2S) Normalization (
    fTotal->SetParameter(7, hUnlike->GetMaximum() * 0.05);
    // fTotal->SetParLimits(7, 0.0, hUnlike->GetMaximum() * 1.0);

    // BG (Expo3)
    fTotal->SetParameter(8, 5.0);
    fTotal->SetParameter(9, -1.0);
    fTotal->SetParameter(10, 0.1);

    // Fitting (L: Log-likelihood, R: Range, S: Save)
    hUnlike->Fit(fTotal, "L R S");

    // ==========================================================
    // J/psi, Psi(2S), Background Drawing
    // ==========================================================

    // Background Expo3
    TF1* fBg = new TF1(Form("fBg_bin%d", i), Expo3, cfg.fitMin, cfg.fitMax, 3);
    fBg->SetParameters(fTotal->GetParameter(8), fTotal->GetParameter(9), fTotal->GetParameter(10));
    fBg->SetLineColor(kGreen);
    fBg->SetLineStyle(1);
    fBg->SetLineWidth(1);
    fBg->Draw("SAME");

    // J/psi DSCB
    TF1* fJpsi = new TF1(Form("fJpsi_bin%d", i), DSCB, cfg.fitMin, cfg.fitMax, 7);
    for (int p = 0; p < 7; ++p) {
      fJpsi->SetParameter(p, fTotal->GetParameter(p));
    }
    fJpsi->SetLineColor(kMagenta);
    fJpsi->SetLineStyle(1);
    fJpsi->SetLineWidth(1);
    fJpsi->Draw("SAME");

    // Psi(2S) DSCB
    TF1* fPsi2S = new TF1(Form("fPsi2S_bin%d", i), DSCB, cfg.fitMin, cfg.fitMax, 7);
    fPsi2S->SetParameter(0, fTotal->GetParameter(7));                                         // N
    fPsi2S->SetParameter(1, fTotal->GetParameter(1) + (cfg.m_Psi2S_PDG - cfg.m_Jpsi_PDG));    // m0
    fPsi2S->SetParameter(2, fTotal->GetParameter(2) * cfg.sigma_factor * cfg.mc_sigma_ratio); // sigma
    fPsi2S->SetParameter(3, cfg.mc_Psi2S_alpha1);
    fPsi2S->SetParameter(4, cfg.mc_Psi2S_n1);
    fPsi2S->SetParameter(5, cfg.mc_Psi2S_alpha2);
    fPsi2S->SetParameter(6, cfg.mc_Psi2S_n2);
    fPsi2S->SetLineColor(kRed);
    fPsi2S->SetLineStyle(1);
    fPsi2S->SetLineWidth(1);
    fPsi2S->Draw("SAME");

    fTotal->Draw("SAME");

    // ==========================================================
    // Draw additional information
    // ==========================================================
    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(fontCode);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.86, Form("%.2f < y_{#mu#mu} < %.2f", ymin, ymax));

    double jpsi_yield = fTotal->GetParameter(0);
    latex.DrawLatex(0.60, 0.86, Form("J/#psi N = %.0f", jpsi_yield));

    double chi2 = fTotal->GetChisquare();
    int ndf = fTotal->GetNDF();
    double chi2_dof = (ndf > 0) ? chi2 / ndf : 0.0;

    latex.DrawLatex(0.60, 0.80, Form("#chi^{2}/dof = %.2f", chi2_dof));
  }

  c1->SaveAs("Mass_Fit_3Bins.png");
}
