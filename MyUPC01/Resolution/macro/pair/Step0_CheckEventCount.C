// ============================================================
//  Step0_CheckEventCount.C
// ============================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

#include <algorithm>
#include <iostream>
#include <vector>

void Step0_CheckEventCount()
{
  // ----------------------------------------------------------
  // Settings
  // ----------------------------------------------------------
  const TString dataDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0408/Resolution/CoherentCut/";
  const TString taskReg = "my-upc-muon-pair-resolution/registry/";
  const TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0414/Coherent/";

  const int rebinFactor = 5;

  struct Sample {
    TString file;
    TString label;
  };
  std::vector<Sample> samples = {
    // {dataDir + "jpsi-incoh.root", "J/#psi incoh."},
   {dataDir + "jpsi-coh.root", "J/#psi coh."},
    // {dataDir + "psi2s-incoh.root", "#psi(2S) incoh."},
    // {dataDir + "psi2s-coh.root", "#psi(2S) coh."},
    // {dataDir + "psi2s-incoh-fd.root", "#psi(2S) incoh. fd"},
    // {dataDir + "psi2s-coh-fd.root", "#psi(2S) coh. fd"},
    // {dataDir + "mumu-low.root", "#mu#mu low"},
    // {dataDir + "mumu-mid.root", "#mu#mu mid"},
    // {dataDir + "mumu-high.root", "#mu#mu high"},
  };

  // ----------------------------------------------------------
  // Global style
  // ----------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  const int fontCode = 42;
  for (const char* ax : {"X", "Y", "Z"}) {
    gStyle->SetLabelFont(fontCode, ax);
    gStyle->SetTitleFont(fontCode, ax);
    gStyle->SetTitleSize(0.048, ax);
    gStyle->SetLabelSize(0.042, ax);
  }

  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                  kOrange + 1, kCyan + 2, kViolet + 1, kTeal + 2, kPink + 1};

  // ==========================================================
  // Helper: merge a 1D histogram across all samples
  // ==========================================================
  auto mergeHist = [&](const TString& histName) -> TH1F* {
    TH1F* hTotal = nullptr;
    for (size_t is = 0; is < samples.size(); ++is) {
      TFile* f = TFile::Open(samples[is].file, "READ");
      if (!f || f->IsZombie()) {
        if (f)
          f->Close();
        continue;
      }
      TH1F* h = (TH1F*)f->Get(taskReg + histName);
      if (!h) {
        f->Close();
        continue;
      }
      h->SetDirectory(0);
      if (!hTotal) {
        hTotal = (TH1F*)h->Clone("h_" + histName + "_total");
        hTotal->SetDirectory(0);
      } else {
        hTotal->Add(h);
      }
      delete h;
      f->Close();
    }
    return hTotal;
  };

  // ==========================================================
  // Helper: draw MC (blue) vs Reco (red) overlaid on one canvas
  // ==========================================================
  auto drawComparison = [&](const TString& histMC,
                            const TString& histReco,
                            const TString& title,
                            const TString& xTitle,
                            const TString& canName,
                            const TString& outFile) {
    TH1F* hMC = mergeHist(histMC);
    TH1F* hReco = mergeHist(histReco);
    if (!hMC || !hReco) {
      std::cerr << "[Error] Could not load histograms for " << title << std::endl;
      return;
    }

    hMC->Rebin(rebinFactor);
    hReco->Rebin(rebinFactor);

    TCanvas* c = new TCanvas(canName, title, 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.06);
    c->SetRightMargin(0.05);
    c->SetGrid();
    c->SetLogy();

    // Style MC
    hMC->SetLineColor(kBlue + 1);
    hMC->SetLineWidth(2);
    hMC->SetMarkerColor(kBlue + 1);
    hMC->SetMarkerStyle(24);
    hMC->SetMarkerSize(0.7);

    // Style Reco
    hReco->SetLineColor(kRed + 1);
    hReco->SetLineWidth(2);
    hReco->SetMarkerColor(kRed + 1);
    hReco->SetMarkerStyle(20);
    hReco->SetMarkerSize(0.7);

    // Axis
    double ymax = std::max(hMC->GetMaximum(), hReco->GetMaximum());
    hMC->GetYaxis()->SetRangeUser(1e0, ymax * 100.0);
    hMC->GetXaxis()->SetRangeUser(0.0, 0.40);
    hMC->GetXaxis()->SetTitle(xTitle);
    hMC->GetYaxis()->SetTitle("Counts");

    hMC->Draw("HIST");
    hReco->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.84, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(fontCode);
    leg->SetTextSize(0.036);
    leg->AddEntry(hMC, Form("MC truth  (%.0f)", hMC->GetEntries()), "l");
    leg->AddEntry(hReco, Form("Reco      (%.0f)", hReco->GetEntries()), "l");
    leg->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.036);
    tex.DrawLatex(0.16, 0.90, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.16, 0.85, "with Coherent J/#psi Cut");
    tex.DrawLatex(0.16, 0.80, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");

    c->SaveAs(outDir + outFile);
    std::cout << "[Info] Saved: " << outDir + outFile << std::endl;

    delete leg;
    delete c;
    delete hMC;
    delete hReco;
  };

  // ----------------------------------------------------------
  // Helper: per-sample overlay for one histogram
  // ----------------------------------------------------------
  auto drawPerSample = [&](const TString& histName,
                           const TString& title,
                           const TString& xTitle,
                           const TString& canName,
                           const TString& outFile) {
    TCanvas* c = new TCanvas(canName, title, 1000, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.06);
    c->SetRightMargin(0.05);
    c->SetGrid();
    c->SetLogy();

    TLegend* leg = new TLegend(0.64, 0.55, 0.94, 0.94);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(fontCode);
    leg->SetTextSize(0.030);

    std::vector<TH1F*> hists;
    std::vector<size_t> sampleIndices;
    double overallMax = 0.0;

    for (size_t is = 0; is < samples.size(); ++is) {
      TFile* f = TFile::Open(samples[is].file, "READ");
      if (!f || f->IsZombie()) {
        if (f)
          f->Close();
        continue;
      }
      TH1F* h = (TH1F*)f->Get(taskReg + histName);
      if (!h) {
        f->Close();
        continue;
      }
      h->SetDirectory(0);
      f->Close();

      h->Rebin(rebinFactor);

      h->SetLineColor(colors[is % 9]);
      h->SetLineWidth(2);
      h->GetXaxis()->SetTitle(xTitle);
      h->GetYaxis()->SetTitle("Counts");

      if (h->GetMaximum() > overallMax) {
        overallMax = h->GetMaximum();
      }

      hists.push_back(h);
      sampleIndices.push_back(is);
    }

    bool firstDraw = true;
    for (size_t i = 0; i < hists.size(); ++i) {
      TH1F* h = hists[i];
      size_t is = sampleIndices[i];

      if (firstDraw) {
        h->GetYaxis()->SetRangeUser(1e0, overallMax * 100.0);
        h->GetXaxis()->SetRangeUser(0.0, 0.40);
        h->Draw("HIST");
        firstDraw = false;
      } else {
        h->Draw("HIST SAME");
      }
      leg->AddEntry(h, Form("%s  (%.0f)", samples[is].label.Data(), h->GetEntries()), "l");
    }

    leg->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.036);
    tex.DrawLatex(0.16, 0.90, Form("#bf{%s (per sample)}", title.Data()));
    tex.DrawLatex(0.16, 0.85, "without Coherent J/#psi Cut");
    tex.DrawLatex(0.16, 0.80, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");

    c->SaveAs(outDir + outFile);
    std::cout << "[Info] Saved: " << outDir + outFile << std::endl;

    delete leg;
    delete c;
    for (auto h : hists) {
      delete h;
    }
  };

  // ==========================================================
  // Figure 1a: Pair pT  MC vs Reco  (PostCut, merged)
  // ==========================================================
  drawComparison("hPairPtMC_PostCut", "hPairPtReco_PostCut",
                 "Pair p_{T} (MC vs Reco)",
                 "p_{T,#mu#mu} (GeV/c)",
                 "c_pt_comparison",
                 "Step0_PairPt_MCvsReco.png");

  // ==========================================================
  // Figure 1b: Pair pT^2  MC vs Reco  (PostCut, merged)
  // ==========================================================
  drawComparison("hPairPt2MC_PostCut", "hPairPt2Reco_PostCut",
                 "Pair p_{T}^{2} (MC vs Reco, |t| #approx p_{T}^{2})",
                 "p_{T,#mu#mu}^{2} (GeV^{2}/c^{2})",
                 "c_pt2_comparison",
                 "Step0_PairPt2_MCvsReco.png");

  // ==========================================================
  // Figure 2a: Pair pT  MC per sample (PostCut)
  // ==========================================================
  drawPerSample("hPairPtMC_PostCut",
                "Pair p_{T}^{MC}",
                "p_{T,#mu#mu} (GeV/c)",
                "c_pt_mc_each",
                "Step0_PairPtMC_each.png");

  // ==========================================================
  // Figure 2b: Pair pT^2  MC per sample (PostCut)
  // ==========================================================
  drawPerSample("hPairPt2MC_PostCut",
                "Pair p_{T}^{2,MC} (|t| #approx p_{T}^{2})",
                "p_{T,#mu#mu}^{2} (GeV^{2}/c^{2})",
                "c_pt2_mc_each",
                "Step0_PairPt2MC_each.png");

  std::cout << "[Done] Step 0 finished." << std::endl;
}
