#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>

// ======================================================================
// Expo3 Background Function
// ======================================================================
/*
double Expo3(double* x, double* p)
{
  double m = x[0];
  return TMath::Exp(p[0] + p[1] * m + p[2] * m * m);
}

double Expo3_Sideband(double* x, double* p)
{
  double m = x[0];
  // J/psi (2.7 - 3.4) - Psi(2S) (3.5 - 3.9)
  if ((m > 2.7 && m < 3.4) || (m > 3.5 && m < 3.9)) {
    TF1::RejectPoint();
    return 0;
  }
  return TMath::Exp(p[0] + p[1] * m + p[2] * m * m);
}
*/

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

double VWG_Sideband(double* x, double* p)
{
  double m = x[0];
  // J/psi (2.7 - 3.4) - Psi(2S) (3.5 - 3.9)
  if ((m > 2.7 && m < 3.4) || (m > 3.5 && m < 3.9)) {
    TF1::RejectPoint();
    return 0;
  }
  return VWG(x, p);
}

// ======================================================================
// Main macro
// ======================================================================
void FitMassBinCounting_VWG(const char* filename = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/0325/AnalysisResults.root")
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

  double ymin = -4.00;
  double ymax = -2.50;
  int binMinY = h2->GetXaxis()->FindBin(ymin + 1e-4);
  int binMaxY = h2->GetXaxis()->FindBin(ymax - 1e-4);

  TH1D* hUnlike = h2->ProjectionY("hMass_Inclusive", binMinY, binMaxY);

  int rebin = 25; // 25 MeV/c^2 per bin
  hUnlike->Rebin(rebin);
  hUnlike->SetTitle("");

  double fitMin = 2.2;
  double fitMax = 4.5;
  double srMin = 2.80; // J/psi Signal Region Min
  double srMax = 3.30; // J/psi Signal Region Max
  double binWidth = hUnlike->GetBinWidth(1);

  // ==========================================================
  // Background Fit
  // ==========================================================
  TF1* fBgFit = new TF1("fBgFit", VWG_Sideband, fitMin, fitMax, 4);
  fBgFit->SetParameters(10.0, -2.0, 0.1);
  hUnlike->Fit(fBgFit, "L R S 0");

  TF1* fBgDraw = new TF1("fBgDraw", VWG, fitMin, fitMax, 4);
  fBgDraw->SetParameters(fBgFit->GetParameters());
  fBgDraw->SetLineColor(kGreen + 2);
  fBgDraw->SetLineStyle(2);
  fBgDraw->SetLineWidth(2);

  // ==========================================================
  // Yield
  // ==========================================================
  int binSRMin = hUnlike->FindBin(srMin + 1e-4);
  int binSRMax = hUnlike->FindBin(srMax - 1e-4);

  // Systematic Error
  double errTotal;
  double nTotal = hUnlike->IntegralAndError(binSRMin, binSRMax, errTotal);

  // Background J/psi (2.7 - 3.4) - Psi(2S) (3.5 - 3.9)
  double nBg = fBgDraw->Integral(srMin, srMax) / binWidth;
  double yieldJpsi = nTotal - nBg;
  double errJpsi = TMath::Sqrt(nTotal + nBg);

  // ==========================================================
  // Draw
  // ==========================================================
  TCanvas* c1 = new TCanvas("c1", "Bin Counting", 900, 750);
  // gPad->SetMargin(0.12, 0.05, 0.12, 0.08);
  gPad->SetMargin(0.14, 0.05, 0.14, 0.08);

  hUnlike->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
  hUnlike->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binWidth * 1000));
  hUnlike->GetXaxis()->SetRangeUser(2.0, 5.0);
  hUnlike->GetYaxis()->SetRangeUser(0, hUnlike->GetMaximum() * 1.5);
  hUnlike->SetMarkerStyle(kFullCircle);
  hUnlike->SetMarkerSize(0.8);
  hUnlike->SetLineColor(kBlack);

  hUnlike->Draw("PE");
  fBgDraw->Draw("SAME");

  TH1D* hSignalRegion = (TH1D*)hUnlike->Clone("hSignalRegion");
  for (int b = 1; b <= hUnlike->GetNbinsX(); ++b) {
    if (b < binSRMin || b > binSRMax) {
      hSignalRegion->SetBinContent(b, 0);
      hSignalRegion->SetBinError(b, 0);
    }
  }
  // hSignalRegion->SetFillColorAlpha(kBlue, 0.3);
  hSignalRegion->SetFillColor(kBlue);
  hSignalRegion->SetFillStyle(3004);
  hSignalRegion->Draw("HIST SAME");
  hUnlike->Draw("PE SAME");

  // ==========================================================
  // Text and Legend
  // ==========================================================
  TLatex latex;
  latex.SetNDC(true);
  latex.SetTextFont(fontCode);
  latex.SetTextSize(0.03);

  latex.DrawLatex(0.65, 0.89, Form("%.2f < y_{#mu#mu} < %.2f", ymin, ymax));
  latex.DrawLatex(0.65, 0.84, "p_{T}^{#mu#mu} < 0.25 GeV/c");

  double chi2 = fBgFit->GetChisquare();
  int ndf = fBgFit->GetNDF();
  double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;

  // latex.DrawLatex(0.20, 0.89, "Bin Counting Method");
  latex.DrawLatex(0.20, 0.89, Form("Signal Region: %.2f - %.2f GeV/#it{c}^{2}", srMin, srMax));
  latex.DrawLatex(0.20, 0.84, Form("N_{Total} = %.0f", nTotal));
  latex.DrawLatex(0.20, 0.79, Form("N_{BG} = %.0f", nBg));
  latex.DrawLatex(0.20, 0.74, Form("#chi^{2}/ndf = %.2f / %d = %.2f", chi2, ndf, chi2_ndf));

  latex.SetTextColor(kRed);
  latex.DrawLatex(0.20, 0.69, Form("N_{J/#psi} = %.0f #pm %.0f", yieldJpsi, errJpsi));

  TLegend* leg = new TLegend(0.65, 0.35, 0.90, 0.55);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(fontCode);
  leg->AddEntry(hUnlike, "Data", "P");
  leg->AddEntry(hSignalRegion, "Signal Region", "F");
  leg->AddEntry(fBgDraw, "Background Fit (VWG)", "L");
  leg->Draw();

  c1->SaveAs("Mass_BinCounting_VWG.png");

  std::cout << "=========================================" << std::endl;
  std::cout << "Bin Counting Results (17% Statistics)" << std::endl;
  std::cout << "Max bin content: " << hUnlike->GetMaximum() << std::endl;
  std::cout << "N_Total : " << nTotal << std::endl;
  std::cout << "N_BG    : " << nBg << std::endl;
  std::cout << "N_J/psi : " << yieldJpsi << " +/- " << errJpsi << std::endl;
  std::cout << "=========================================" << std::endl;
}
