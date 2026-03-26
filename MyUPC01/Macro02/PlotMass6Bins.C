#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>

void PlotMass6Bins(const char* filename = "AnalysisResults.root")
{
  gStyle->SetOptStat(0);

  int fontCode = 42;
  gStyle->SetLabelFont(fontCode, "XYZ");
  gStyle->SetTitleFont(fontCode, "XYZ");
  gStyle->SetTextFont(fontCode);
  gStyle->SetLegendFont(fontCode);

  gStyle->SetTitleFont(12, "t");

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Error: " << filename << " can't open" << std::endl;
    return;
  }
  TH2D* h2 = (TH2D*)f->Get("my-upc-mass-02/registry/hMassVsRapidityUnlike");
  if (!h2) {
    std::cerr << "Error: 2D histogram 'hMassVsRapidityUnlike' not found." << std::endl;
    return;
  }

  TCanvas* c1 = new TCanvas("c1", "Invariant Mass in Rapidity Bins", 800, 1200);
  c1->Divide(2, 3, 1e-5, 1e-5);

  double y_bins[7] = {-4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50};

  for (int i = 0; i < 6; ++i) {
    c1->cd(i + 1);

    gPad->SetMargin(0.15, 0.07, 0.15, 0.07);

    double ymin = y_bins[i];
    double ymax = y_bins[i + 1];
    double drawMin = 2.0;
    double drawMax = 5.0;
    int rebin = 25;

    int binMin = h2->GetXaxis()->FindBin(ymin + 1e-4);
    int binMax = h2->GetXaxis()->FindBin(ymax - 1e-4);

    TString histName = Form("hMass_bin%d", i);
    TH1D* hUnlike = h2->ProjectionY(histName, binMin, binMax);

    hUnlike->Rebin(rebin);
    hUnlike->SetTitle(Form("%.2f < y < %.2f", ymin, ymax));
    hUnlike->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
    hUnlike->GetYaxis()->SetTitle(Form("Counts / (%d MeV/#it{c}^{2})", rebin));
    hUnlike->GetXaxis()->SetRangeUser(drawMin, drawMax);
    hUnlike->GetYaxis()->SetRangeUser(0, hUnlike->GetMaximum() * 1.5);

    hUnlike->SetLineColor(kBlack);
    hUnlike->SetMarkerColor(kBlack);
    hUnlike->SetMarkerStyle(kFullCircle);
    hUnlike->SetMarkerSize(0.5);
    hUnlike->SetLineWidth(1);

    hUnlike->GetXaxis()->SetTitleSize(0.045);
    hUnlike->GetYaxis()->SetTitleSize(0.045);
    hUnlike->GetXaxis()->SetLabelSize(0.045);
    hUnlike->GetYaxis()->SetLabelSize(0.045);

    hUnlike->GetYaxis()->SetTitleOffset(1.3);

    hUnlike->Draw("PE");

    double counts = hUnlike->Integral();

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(fontCode);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.80, Form("%.2f < y_{#mu#mu} < %.2f", ymin, ymax));
    latex.DrawLatex(0.60, 0.89, Form("Unlike = %.0f", counts));
  }

  c1->SaveAs("Mass_Rapidity_6Bins.png");
}
