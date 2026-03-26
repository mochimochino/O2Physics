#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>

void PlotMass3Bins(const char* filename = "AnalysisResults.root")
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

  TCanvas* c1 = new TCanvas("c1", "Invariant Mass in Rapidity Bins", 800, 800);

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

    hUnlike->GetYaxis()->SetTitleOffset(1.8);

    hUnlike->Draw("PE");

    double counts = hUnlike->Integral();

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextFont(fontCode);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.86, Form("%.2f < y_{#mu#mu} < %.2f", ymin, ymax));
    latex.DrawLatex(0.60, 0.86, Form("Unlike = %.0f", counts));
  }

  c1->SaveAs("Mass_Rapidity_3Bins.png");
}
