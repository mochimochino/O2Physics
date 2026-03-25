#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
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

  // J/psi Tail parameters from MC
  double mc_Jpsi_alpha1[3] = {1.091, 0.969, 0.944};  // Left tail alpha
  double mc_Jpsi_n1[3] = {21.436, 30.811, 20.212};   // Left tail n
  double mc_Jpsi_alpha2[3] = {8.551, 11.939, 3.683}; // Right tail alpha
  double mc_Jpsi_n2[3] = {41.997, 5.289, 5.609};     // Right tail n

  // Psi(2S) Tail parameters from MC
  double mc_Psi2S_alpha1[3] = {1.143, 1.108, 1.025};
  double mc_Psi2S_n1[3] = {26.046, 13.743, 15.865};
  double mc_Psi2S_alpha2[3] = {7.781, 7.479, 1904.071};
  double mc_Psi2S_n2[3] = {27.004, 9.630, 12.503};

  // Resolution Parameters
  double mc_sigma_Jpsi[3] = {0.064, 0.065, 0.066};
  double mc_sigma_Psi2S[3] = {0.07, 0.071, 0.072};
  double mc_sigma_ratio[3] = {mc_sigma_Jpsi[0] / mc_sigma_Psi2S[0], mc_sigma_Jpsi[1] / mc_sigma_Psi2S[1], mc_sigma_Jpsi[2] / mc_sigma_Psi2S[2]};
  double sigma_factor[3] = {1.1, 1.1, 1.1};

  // Fit range
  double fitMin = 2.5;
  double fitMax = 4.5;

  // Background function type (1: Expo3, 2: ExpoPol4, 3: VWG)
  int bgType = 1;

  int current_bin = 0;
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
  return p[0] * TMath::Exp(p[1] * m + p[2] * m * m);
}

// ExpoPol4 Functions
double ExpoPol4(double* x, double* p)
{
  double m = x[0];
  return p[0] * TMath::Exp((p[1] * m) * (p[2] + p[3] * m + p[4] * m * m + p[5] * m * m * m + p[6] * m * m * m * m));
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

// J/psi + Psi(2S) Signal Functions
double SignalJpsiPsi2s(double* x, double* p)
{
  // J/psi parameters
  double p_Jpsi[7];
  for (int j = 0; j < 7; j++)
    p_Jpsi[j] = p[j];

  int bin = cfg.current_bin;

  // Psi(2S) parameters
  double p_Psi2S[7];
  p_Psi2S[0] = p[7];
  p_Psi2S[1] = p_Jpsi[1] + (cfg.m_Psi2S_PDG - cfg.m_Jpsi_PDG);
  p_Psi2S[2] = p_Jpsi[2] * cfg.sigma_factor[bin] * cfg.mc_sigma_ratio[bin];
  p_Psi2S[3] = cfg.mc_Psi2S_alpha1[bin];
  p_Psi2S[4] = cfg.mc_Psi2S_n1[bin];
  p_Psi2S[5] = cfg.mc_Psi2S_alpha2[bin];
  p_Psi2S[6] = cfg.mc_Psi2S_n2[bin];

  return DSCB(x, p_Jpsi) + DSCB(x, p_Psi2S);
}

// Total Fit FUnction: J/psi + Psi(2S) + Background
double TotalFit(double* x, double* p)
{
  double signal = SignalJpsiPsi2s(x, p);
  double bg = 0;

  if (cfg.bgType == 1) {
    bg = Expo3(x, &p[8]);
  }
  return signal + bg;
}

// ======================================================================
// Main macro
// ======================================================================
void FitMass3Bins(const char* filename = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/0319/AnalysisResults.root")
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

  double total_jpsi_yield = 0.0;
  double total_jpsi_err_sq = 0.0;

  double total_psi2s_yield = 0.0;
  double total_psi2s_err_sq = 0.0;

  for (int i = 0; i < 3; ++i) {
    cfg.current_bin = i;

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

    // ==========================================================
    // Outside
    // ==========================================================
    TH1D* hOutside = (TH1D*)hUnlike->Clone(Form("hOutside_bin%d", i));
    for (int b = 1; b <= hOutside->GetNbinsX(); ++b) {
      double binCenter = hOutside->GetBinCenter(b);
      if (binCenter >= cfg.fitMin && binCenter <= cfg.fitMax) {
        hOutside->SetBinContent(b, 0);
        hOutside->SetBinError(b, 0);
      }
    }
    hOutside->SetFillColorAlpha(kGray + 1, 0.8);
    hOutside->SetFillStyle(3004);
    hOutside->SetLineColor(kGray + 1);
    hOutside->SetMarkerSize(0);

    // ==========================================================
    // Fitting Setting and Execution
    // ==========================================================
    TF1* fTotal = new TF1(Form("fTotal_bin%d", i), TotalFit, cfg.fitMin, cfg.fitMax, 11);
    fTotal->SetNpx(1000);
    fTotal->SetLineColor(kBlue);
    fTotal->SetLineWidth(1);

    // J/psi parameters
    fTotal->SetParameter(0, hUnlike->GetMaximum());
    fTotal->SetParameter(1, cfg.m_Jpsi_PDG);
    fTotal->SetParLimits(1, 2.9, 3.2);
    fTotal->SetParameter(2, cfg.mc_sigma_Jpsi[i]);
    // fTotal->SetParLimits(2, 0.04, 0.15);

    fTotal->SetParameter(3, cfg.mc_Jpsi_alpha1[i]);
    // fTotal->SetParLimits(3, 0.1, 5.0);
    fTotal->SetParameter(4, cfg.mc_Jpsi_n1[i]);
    fTotal->SetParLimits(4, 1.0, 30.0);
    fTotal->FixParameter(5, cfg.mc_Jpsi_alpha2[i]);
    fTotal->FixParameter(6, cfg.mc_Jpsi_n2[i]);

    // Psi(2S) parameters
    fTotal->SetParameter(7, hUnlike->GetMaximum() * 0.1); // N

    // fTotal->SetParameter(8, 10.7128);
    // fTotal->SetParameter(9, -2.15838);
    // fTotal->SetParameter(10, 0.136308);

    fTotal->SetParameter(8, 15000.0);
    fTotal->SetParameter(9, -2.15838);
    fTotal->SetParameter(10, 0.136308);

    hUnlike->Fit(fTotal, "L R S B 0");

    // ==========================================================
    // Draw
    // ==========================================================
    hUnlike->Draw("PE");
    hOutside->Draw("HIST SAME");
    hUnlike->Draw("PE SAME");

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
    fPsi2S->SetParameter(0, fTotal->GetParameter(7));
    fPsi2S->SetParameter(1, fTotal->GetParameter(1) + (cfg.m_Psi2S_PDG - cfg.m_Jpsi_PDG));
    fPsi2S->SetParameter(2, fTotal->GetParameter(2) * cfg.sigma_factor[i] * cfg.mc_sigma_ratio[i]);
    fPsi2S->SetParameter(3, cfg.mc_Psi2S_alpha1[i]);
    fPsi2S->SetParameter(4, cfg.mc_Psi2S_n1[i]);
    fPsi2S->SetParameter(5, cfg.mc_Psi2S_alpha2[i]);
    fPsi2S->SetParameter(6, cfg.mc_Psi2S_n2[i]);
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

    // double jpsi_yield = fTotal->GetParameter(0);
    // latex.DrawLatex(0.60, 0.88, "p_{T}^{#mu#mu} < 0.25 GeV/c");
    // latex.DrawLatex(0.60, 0.82, Form("N_{J/#psi} = %.0f", jpsi_yield));

    double binWidth = hUnlike->GetBinWidth(1);

    double jpsi_yield = fJpsi->Integral(cfg.fitMin, cfg.fitMax) / binWidth;
    double yield_err = jpsi_yield * (fTotal->GetParError(0) / fTotal->GetParameter(0));
    latex.DrawLatex(0.60, 0.82, Form("N_{J/#psi} = %.0f #pm %.0f", jpsi_yield, yield_err));
    std::cout << Form("Bin %d (%.2f < y < %.2f): N_J/psi = %.1f +/- %.1f",
                      i, ymin, ymax, jpsi_yield, yield_err)
              << std::endl;

    total_jpsi_yield += jpsi_yield;
    total_jpsi_err_sq += (yield_err * yield_err);

    double psi2s_yield = fPsi2S->Integral(cfg.fitMin, cfg.fitMax) / binWidth;
    double psi2s_err = psi2s_yield * (fTotal->GetParError(7) / fTotal->GetParameter(7));
    latex.DrawLatex(0.60, 0.76, Form("N_{#psi(2S)} = %.0f #pm %.0f", psi2s_yield, psi2s_err));
    std::cout << Form("Bin %d (%.2f < y < %.2f): N_psi(2S) = %.1f +/- %.1f",
                      i, ymin, ymax, psi2s_yield, psi2s_err)
              << std::endl;

    total_psi2s_yield += psi2s_yield;
    total_psi2s_err_sq += (psi2s_err * psi2s_err);

    double chi2 = fTotal->GetChisquare();
    int ndf = fTotal->GetNDF();
    double chi2_dof = (ndf > 0) ? chi2 / ndf : 0.0;
    latex.DrawLatex(0.60, 0.70, Form("#chi^{2}/dof = %.2f", chi2_dof));

    TLegend* legend = new TLegend(0.60, 0.30, 0.85, 0.63);
    legend->SetBorderSize(0);
    legend->SetTextFont(fontCode);
    legend->SetTextSize(0.035);
    legend->AddEntry(hUnlike, "Data", "P");
    legend->AddEntry(hOutside, "Continuum (Raw)", "F");
    legend->AddEntry(fJpsi, "J/#psi", "L");
    legend->AddEntry(fPsi2S, "#psi(2S)", "L");
    legend->AddEntry(fBg, "Background (Expo3)", "L");
    legend->AddEntry(fTotal, "Total Fit", "L");
    legend->Draw();
  }

  c1->SaveAs("Mass_Fit_3Bins.png");
  double total_jpsi_err = TMath::Sqrt(total_jpsi_err_sq);
  double total_psi2s_err = TMath::Sqrt(total_psi2s_err_sq);
  std::cout << "=========================================" << std::endl;
  std::cout << "Total J/psi Yield (Sum of 3 bins)" << std::endl;
  std::cout << Form("N_J/psi_total = %.1f +/- %.1f", total_jpsi_yield, total_jpsi_err) << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << "Total Psi(2S) Yield (Sum of 3 bins)" << std::endl;
  std::cout << Form("N_psi(2S)_total = %.1f +/- %.1f", total_psi2s_yield, total_psi2s_err) << std::endl;
  std::cout << "=========================================" << std::endl;
}
