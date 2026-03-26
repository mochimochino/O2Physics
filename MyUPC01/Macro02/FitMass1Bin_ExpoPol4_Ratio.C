#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
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
  double mc_Jpsi_alpha1[1] = {1.000}; // Left tail alpha
  double mc_Jpsi_n1[1] = {24.000};    // Left tail n
  double mc_Jpsi_alpha2[1] = {8.000}; // Right tail alpha
  double mc_Jpsi_n2[1] = {17.000};    // Right tail n

  // Psi(2S) Tail parameters from MC
  double mc_Psi2S_alpha1[1] = {1.100};
  double mc_Psi2S_n1[1] = {18.000};
  double mc_Psi2S_alpha2[1] = {7.500};
  double mc_Psi2S_n2[1] = {16.000};

  // Resolution Parameters
  double mc_sigma_Jpsi[1] = {0.065};
  double mc_sigma_Psi2S[1] = {0.071};

  double mc_sigma_ratio[1] = {mc_sigma_Jpsi[0] / mc_sigma_Psi2S[0]};
  double sigma_factor[1] = {1.1};

  // Fit range
  double fitMin = 2.5;
  double fitMax = 4.5;
  int bgType = 1; // 1: Expo3
  int current_bin = 0;
};

FitConfig cfg;

// ======================================================================
// Fit Functions
// ======================================================================
double DSCB(double* x, double* p)
{
  double m = x[0], N = p[0], m0 = p[1], sigma = p[2];
  double alpha1 = p[3], n1 = p[4], alpha2 = p[5], n2 = p[6];
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

double Expo3(double* x, double* p)
{
  return p[0] * TMath::Exp(p[1] * x[0] + p[2] * x[0] * x[0]);
}

// ExpoPol4 Functions
double ExpoPol4(double* x, double* p)
{
  double m = x[0];
  return p[0] * TMath::Exp((p[1] * m) * (p[2] + p[3] * m + p[4] * m * m + p[5] * m * m * m + p[6] * m * m * m * m));
}

double SignalJpsiPsi2s(double* x, double* p)
{
  double p_Jpsi[7];
  for (int j = 0; j < 7; j++)
    p_Jpsi[j] = p[j];
  int bin = cfg.current_bin;
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

double TotalFit(double* x, double* p)
{
  return SignalJpsiPsi2s(x, p) + ExpoPol4(x, &p[8]);
}

// ======================================================================
// Main macro
// ======================================================================
void FitMass1Bin_ExpoPol4_Ratio(const char* filename = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/0325/AnalysisResults.root")
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

  TCanvas* c1 = new TCanvas("c1", "Fit Invariant Mass with Ratio", 900, 900);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetTopMargin(0.08);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.05);
  pad1->Draw();

  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.05);
  pad2->Draw();

  // -----------------------------------------------------------------
  // Pad1
  // -----------------------------------------------------------------
  pad1->cd();

  double ymin = -4.00;
  double ymax = -2.50;
  int rebin = 25;

  int binMin = h2->GetXaxis()->FindBin(ymin + 1e-4);
  int binMax = h2->GetXaxis()->FindBin(ymax - 1e-4);
  TH1D* hUnlike = h2->ProjectionY("hMass", binMin, binMax);
  hUnlike->Rebin(rebin);
  hUnlike->SetTitle("");

  hUnlike->GetXaxis()->SetLabelSize(0);
  hUnlike->GetXaxis()->SetTitleSize(0);
  hUnlike->GetYaxis()->SetTitle(Form("Counts / (%d MeV/#it{c}^{2})", rebin));
  hUnlike->GetXaxis()->SetRangeUser(2.0, 5.0);
  hUnlike->GetYaxis()->SetRangeUser(0, hUnlike->GetMaximum() * 1.5);

  hUnlike->SetMarkerStyle(kFullCircle);
  hUnlike->SetMarkerSize(0.8);
  hUnlike->SetLineColor(kBlack);
  hUnlike->GetYaxis()->SetTitleSize(0.05);
  hUnlike->GetYaxis()->SetLabelSize(0.045);
  hUnlike->GetYaxis()->SetTitleOffset(1.3);

  // Continuum Background
  TH1D* hOutside = (TH1D*)hUnlike->Clone("hOutside");
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

  // Fit
  TF1* fTotal = new TF1("fTotal", TotalFit, cfg.fitMin, cfg.fitMax, 15);
  fTotal->SetNpx(1000);
  fTotal->SetLineColor(kBlue);
  fTotal->SetLineWidth(2);

  fTotal->SetParameter(0, hUnlike->GetMaximum());
  fTotal->SetParameter(1, cfg.m_Jpsi_PDG);
  fTotal->SetParLimits(1, 2.9, 3.2);
  fTotal->SetParameter(2, cfg.mc_sigma_Jpsi[0]);

  fTotal->SetParameter(3, cfg.mc_Jpsi_alpha1[0]);
  fTotal->SetParameter(4, cfg.mc_Jpsi_n1[0]);
  fTotal->SetParLimits(4, 1.0, 30.0);
  fTotal->FixParameter(5, cfg.mc_Jpsi_alpha2[0]);
  fTotal->FixParameter(6, cfg.mc_Jpsi_n2[0]);

  fTotal->SetParameter(7, hUnlike->GetMaximum() * 0.1);
  fTotal->SetParLimits(7, 0.0, hUnlike->GetMaximum());

  fTotal->SetParameter(8, 30000.0);
  fTotal->SetParLimits(8, 0.0, hUnlike->GetMaximum() * 10.0);
  fTotal->SetParameter(9, -1.0);
  fTotal->SetParameter(10, 1.0);
  fTotal->SetParameter(11, 0.0);
  fTotal->SetParameter(12, 0.0);
  fTotal->SetParameter(13, 0.0);
  fTotal->SetParameter(14, 0.0);

  TFitResultPtr r = hUnlike->Fit(fTotal, "L R S B 0");

  hUnlike->Draw("PE");
  hOutside->Draw("HIST SAME");
  hUnlike->Draw("PE SAME");

  TF1* fBg = new TF1("fBg", ExpoPol4, cfg.fitMin, cfg.fitMax, 7);
  for (int i = 0; i < 7; i++) {
    fBg->SetParameter(i, fTotal->GetParameter(i + 8));
  }
  fBg->SetLineColor(kGreen);
  fBg->SetLineStyle(2);
  fBg->SetLineWidth(2);
  fBg->Draw("SAME");

  TF1* fJpsi = new TF1("fJpsi", DSCB, cfg.fitMin, cfg.fitMax, 7);
  for (int p = 0; p < 7; ++p)
    fJpsi->SetParameter(p, fTotal->GetParameter(p));
  fJpsi->SetLineColor(kMagenta);
  fJpsi->SetLineWidth(2);
  fJpsi->Draw("SAME");

  TF1* fPsi2S = new TF1("fPsi2S", DSCB, cfg.fitMin, cfg.fitMax, 7);
  fPsi2S->SetParameter(0, fTotal->GetParameter(7));
  fPsi2S->SetParameter(1, fTotal->GetParameter(1) + (cfg.m_Psi2S_PDG - cfg.m_Jpsi_PDG));
  fPsi2S->SetParameter(2, fTotal->GetParameter(2) * cfg.sigma_factor[0] * cfg.mc_sigma_ratio[0]);
  fPsi2S->SetParameter(3, cfg.mc_Psi2S_alpha1[0]);
  fPsi2S->SetParameter(4, cfg.mc_Psi2S_n1[0]);
  fPsi2S->SetParameter(5, cfg.mc_Psi2S_alpha2[0]);
  fPsi2S->SetParameter(6, cfg.mc_Psi2S_n2[0]);
  fPsi2S->SetLineColor(kRed);
  fPsi2S->SetLineWidth(2);
  fPsi2S->Draw("SAME");

  fTotal->Draw("SAME");

  TLatex latex;
  latex.SetNDC(true);
  latex.SetTextFont(fontCode);
  latex.SetTextSize(0.045);

  latex.DrawLatex(0.20, 0.82, Form("%.2f < y_{#mu#mu} < %.2f", ymin, ymax));
  latex.DrawLatex(0.65, 0.82, "p_{T}^{#mu#mu} < 0.25 GeV/c");

  double binWidth = hUnlike->GetBinWidth(1);
  double jpsi_yield = fJpsi->Integral(cfg.fitMin, cfg.fitMax) / binWidth;
  double jpsi_err = jpsi_yield * (fTotal->GetParError(0) / fTotal->GetParameter(0));
  latex.DrawLatex(0.65, 0.76, Form("N_{J/#psi} = %.0f #pm %.0f", jpsi_yield, jpsi_err));

  double psi2s_yield = fPsi2S->Integral(cfg.fitMin, cfg.fitMax) / binWidth;
  double psi2s_err = psi2s_yield * (fTotal->GetParError(7) / fTotal->GetParameter(7));
  latex.DrawLatex(0.65, 0.70, Form("N_{#psi(2S)} = %.0f #pm %.0f", psi2s_yield, psi2s_err));

  double chi2 = fTotal->GetChisquare();
  int ndf = fTotal->GetNDF();
  double chi2_dof = (ndf > 0) ? chi2 / ndf : 0.0;
  latex.DrawLatex(0.65, 0.64, Form("#chi^{2}/dof = %.2f", chi2_dof));

  TLegend* legend = new TLegend(0.65, 0.34, 0.89, 0.59);
  legend->SetBorderSize(0);
  legend->SetTextFont(fontCode);
  legend->SetTextSize(0.03);
  legend->AddEntry(hUnlike, "Data", "P");
  legend->AddEntry(hOutside, "Continuum", "F");
  legend->AddEntry(fJpsi, "J/#psi", "L");
  legend->AddEntry(fPsi2S, "#psi(2S)", "L");
  legend->AddEntry(fBg, "Background (ExpoPol4)", "L");
  legend->AddEntry(fTotal, "Total Fit", "L");
  legend->Draw();

  // -----------------------------------------------------------------
  // Pad2
  // -----------------------------------------------------------------
  pad2->cd();
  pad2->SetGridy();

  TH1D* hRatio = (TH1D*)hUnlike->Clone("hRatio");
  hRatio->SetTitle("");

  // (Data - Fit) / Fit
  for (int i = 1; i <= hRatio->GetNbinsX(); i++) {
    double x_center = hRatio->GetBinCenter(i);
    if (x_center >= cfg.fitMin && x_center <= cfg.fitMax) {
      double data = hUnlike->GetBinContent(i);
      double err = hUnlike->GetBinError(i);
      double fit_val = fTotal->Eval(x_center);

      if (fit_val > 0) {
        hRatio->SetBinContent(i, (data - fit_val) / fit_val);
        hRatio->SetBinError(i, err / fit_val);
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
  hRatio->SetMarkerStyle(20);
  hRatio->GetYaxis()->SetTitle("[Data - Fit]/Fit");
  hRatio->GetYaxis()->SetRangeUser(-0.25, 0.25);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetTitleSize(0.10);
  hRatio->GetYaxis()->SetTitleOffset(0.5);

  hRatio->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
  hRatio->GetXaxis()->SetLabelSize(0.10);
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetTitleOffset(1.0);

  hRatio->Draw("EP");

  TLine* line0 = new TLine(2.0, 0.0, 5.0, 0.0);
  line0->SetLineColor(kRed);
  line0->SetLineStyle(2);
  line0->Draw("SAME");

  c1->SaveAs("Mass_Fit_1Bin_ExpoPol4_Ratio.png");

  pad1->cd();
  pad1->SetLogy();
  hUnlike->GetYaxis()->SetRangeUser(10.0, hUnlike->GetMaximum() * 10.0);
  pad1->Update();
  c1->SaveAs("Mass_Fit_1Bin_ExpoPol4_Ratio_Log.png");

  // =========================================================
  // Error
  // =========================================================
  double R_N = psi2s_yield / jpsi_yield;
  double rel_err_psi2s = psi2s_err / psi2s_yield;
  double rel_err_jpsi = jpsi_err / jpsi_yield;

  // Wihout covariance
  double err_independent = R_N * TMath::Sqrt(rel_err_psi2s * rel_err_psi2s + rel_err_jpsi * rel_err_jpsi);

  // With covariance
  double cov_0_7 = r->CovMatrix(0, 7);
  double param_0 = fTotal->GetParameter(0);
  double param_7 = fTotal->GetParameter(7);
  double rel_cov = cov_0_7 / (param_0 * param_7);

  double err_covariance = R_N * TMath::Sqrt(rel_err_psi2s * rel_err_psi2s + rel_err_jpsi * rel_err_jpsi - 2 * rel_cov);

  // ==========================================================
  // Print Total Yields to Terminal
  // ==========================================================
  std::cout << "=========================================" << std::endl;
  std::cout << "Fit Results (-4.00 < y < -2.50)" << std::endl;
  std::cout << Form("N_J/psi   = %.1f +/- %.1f", jpsi_yield, jpsi_err) << std::endl;
  std::cout << Form("N_psi(2S) = %.1f +/- %.1f", psi2s_yield, psi2s_err) << std::endl;
  std::cout << Form("Chi2/NDF  = %.2f", chi2_dof) << std::endl;
  std::cout << "Max bin content: " << hUnlike->GetMaximum() << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "[ Ratio Error Comparison]" << std::endl;
  std::cout << Form("Ratio (R)       = %.5f", R_N) << std::endl;
  std::cout << Form("1. User Formula = %.5f (Independent)", err_independent) << std::endl;
  std::cout << Form("2. Cov Formula  = %.5f (With Covariance)", err_covariance) << std::endl;
  std::cout << Form("-> Covariance(0,7) value = %e", cov_0_7) << std::endl;
  std::cout << "=========================================" << std::endl;
}
