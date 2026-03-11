#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>

void AnalysisPtReal()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);

  // Open the file
  TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/AnalysisResults.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Cannot open AnalysisResults.root" << std::endl;
    return;
  }

  TString dir = "my-upc-mass-jpsi/";

  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dir + "PtMuonUnlike"));
  TH1* hLike = dynamic_cast<TH1*>(f->Get(dir + "PtMuonLike"));

  if (!hUnlike || !hLike) {
    std::cerr << "Error: Histograms PtMuonUnlike or PtMuonLike not found in " << dir << std::endl;
    f->ls();
    return;
  }

  // Colors from AnalysisMassReal.C
  int colUnlike = TColor::GetColor("#2e8b57"); // SeaGreen
  int colLike = TColor::GetColor("#d2691e");   // Chocolate

  // ==========================================
  // 1. Raw Histograms (Unlike and Like Overlaid)
  // ==========================================
  TCanvas* c1 = new TCanvas("c1", "Raw Dimuon pT", 800, 600);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.12);

  // gPad->SetLogy(); // y-axis to log scale
  gPad->SetGrid(); // show grid

  // Unlike Sign settings
  hUnlike->SetTitle("Dimuon p_{T} ALL; p_{T} [GeV/c]; Counts");
  hUnlike->SetLineColor(colUnlike);
  hUnlike->SetMarkerColor(colUnlike);
  hUnlike->SetLineWidth(2);
  hUnlike->SetMarkerStyle(21); // Full Square
  hUnlike->GetYaxis()->SetTitleOffset(1.2);
  hUnlike->GetXaxis()->SetRangeUser(0.0, 1.5);
  hUnlike->GetYaxis()->SetRangeUser(0.8, hUnlike->GetMaximum() * 5.0);
  hUnlike->Draw("E");

  // Like Sign settings
  hLike->SetLineColor(colLike);
  hLike->SetMarkerColor(colLike);
  hLike->SetLineWidth(2);
  hLike->SetMarkerStyle(25); // Open Square
  hLike->Draw("E SAME");

  // Legend
  TLegend* leg = new TLegend(0.65, 0.80, 0.95, 0.95);
  leg->SetBorderSize(1);
  leg->AddEntry(hUnlike, "Unlike Sign (+-)", "lep");
  leg->AddEntry(hLike, "Like Sign (++, --)", "lep");
  leg->Draw();

  c1->SaveAs("RawPt_Dimuon_ALL.png");

  // ==========================================
  // 2. Signal Extraction (Unlike - Like)
  // ==========================================
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  hSignal->SetTitle("Signal Dimuon p_{T} ALL; p_{T} [GeV/c]; Counts");
  hSignal->Add(hLike, -1.0);

  TCanvas* c2 = new TCanvas("c2", "Signal Dimuon pT", 800, 600);
  c2->SetLeftMargin(0.12);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.05);
  c2->SetBottomMargin(0.12);

  gPad->SetGrid();
  // gPad->SetLogy();

  double ymin = hSignal->GetMinimum();
  if (ymin < 0.1)
    ymin = 0.1;

  hSignal->SetLineColor(colUnlike);
  hSignal->SetMarkerColor(colUnlike);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetYaxis()->SetRangeUser(ymin, hSignal->GetMaximum() * 2.0);
  hSignal->GetXaxis()->SetRangeUser(0.0, 1.5);
  hSignal->Draw("E");

  TLegend* leg2 = new TLegend(0.65, 0.85, 0.95, 0.95);
  leg2->SetBorderSize(1);
  leg2->AddEntry(hSignal, "Unlike - Like", "lep");
  leg2->Draw();

  c2->SaveAs("SignalPt_Dimuon_ALL.png");

  std::cout << "Done. Saved plots to RawPt_Dimuon_ALL.png and SignalPt_Dimuon_ALL.png" << std::endl;

  f->Close();
}
