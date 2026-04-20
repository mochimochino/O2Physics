#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TString.h>

#include <algorithm>
#include <iostream>

void drawpT_compare(const char* filename = "AnalysisResults.root")
{
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return;
  }

  TH1F* h_pt_reco = (TH1F*)f->Get("hPtReco");
  if (!h_pt_reco)
    h_pt_reco = (TH1F*)f->Get("my-upc-muon-resolution/registry/hPtReco");

  TH1F* h_pt_true = (TH1F*)f->Get("hPtTrue");
  if (!h_pt_true)
    h_pt_true = (TH1F*)f->Get("my-upc-muon-resolution/registry/hPtTrue");

  if (!h_pt_reco || !h_pt_true) {
    std::cerr << "Error: Could not find hPtReco or hPtTrue in " << filename << std::endl;
    f->Close();
    return;
  }

  TCanvas* c1 = new TCanvas("c1", "pT Comparison Canvas", 800, 600);
  c1->SetLeftMargin(0.13);

  h_pt_true->SetStats(0);
  h_pt_reco->SetStats(0);

  h_pt_true->SetLineColor(kRed);
  h_pt_true->SetLineWidth(2);

  h_pt_reco->SetLineColor(kBlue);
  h_pt_reco->SetLineWidth(2);

  double maxTrue = h_pt_true->GetMaximum();
  double maxReco = h_pt_reco->GetMaximum();
  h_pt_true->SetMaximum(std::max(maxTrue, maxReco) * 1.2);
  h_pt_true->GetXaxis()->SetRangeUser(0, 3.0);

  h_pt_true->SetTitle("Single p_{T} Distribution;p_{T} (GeV/c);Entries");

  h_pt_true->Draw("HIST");
  h_pt_reco->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.15, 0.75, 0.48, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TString trueLabel = Form("p_{T} True (Entries: %lld)", (long long)h_pt_true->GetEntries());
  TString recoLabel = Form("p_{T} Reco (Entries: %lld)", (long long)h_pt_reco->GetEntries());

  leg->AddEntry(h_pt_true, trueLabel, "l");
  leg->AddEntry(h_pt_reco, recoLabel, "l");
  leg->Draw();

  c1->SaveAs("pT_comparison.png");

  std::cout << "Successfully generated pT_comparison.png" << std::endl;
}
