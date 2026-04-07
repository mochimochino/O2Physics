// ============================================================
//  Step1_MergeAndView2D.C
// ============================================================
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>
#include <vector>

void Step1_MergeAndView2D()
{
  // ----------------------------------------------------------
  // Settings
  // ----------------------------------------------------------
  const TString dataDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0403/without/";
  const TString taskReg = "my-upc-muon-pair-resolution/registry/";
  const TString outDir = dataDir;

  std::vector<TString> fileNames = {
    dataDir + "jpsi-incoh.root",
    dataDir + "jpsi-coh.root",
    dataDir + "psi2s-incoh.root",
    dataDir + "psi2s-coh.root",
    dataDir + "psi2s-incoh-fd.root",
    dataDir + "psi2s-coh-fd.root",
    dataDir + "mumu-low.root",
    dataDir + "mumu-mid.root",
    dataDir + "mumu-high.root",
  };

  // ----------------------------------------------------------
  // Global style
  // ----------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(256);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  const int fontCode = 42;
  for (const char* ax : {"X", "Y", "Z"}) {
    gStyle->SetLabelFont(fontCode, ax);
    gStyle->SetTitleFont(fontCode, ax);
    gStyle->SetTitleSize(0.045, ax);
    gStyle->SetLabelSize(0.040, ax);
  }

  // ----------------------------------------------------------
  // Helper: merge a 2D histogram from all files
  // ----------------------------------------------------------
  auto merge2D = [&](const TString& histName) -> TH2F* {
    TH2F* h2 = nullptr;
    for (const auto& fname : fileNames) {
      TFile* f = TFile::Open(fname, "READ");
      if (!f || f->IsZombie()) {
        std::cerr << "[Warning] Cannot open: " << fname << std::endl;
        if (f)
          f->Close();
        continue;
      }
      TH2F* htmp = (TH2F*)f->Get(taskReg + histName);
      if (!htmp) {
        std::cerr << "[Warning] " << histName << " not found in " << fname << std::endl;
        f->Close();
        continue;
      }
      if (!h2) {
        h2 = (TH2F*)htmp->Clone(histName + "_merged");
        h2->SetDirectory(0);
      } else {
        h2->Add(htmp);
      }
      f->Close();
      std::cout << "[Info] Loaded " << histName << " from " << fname << std::endl;
    }
    return h2;
  };

  // ----------------------------------------------------------
  // Helper: merge a 1D histogram from all files
  // ----------------------------------------------------------
  auto merge1D = [&](const TString& histName) -> TH1D* {
    TH1D* h1 = nullptr;
    for (const auto& fname : fileNames) {
      TFile* f = TFile::Open(fname, "READ");
      if (!f || f->IsZombie()) {
        if (f)
          f->Close();
        continue;
      }
      TH1D* htmp = (TH1D*)f->Get(taskReg + histName);
      if (!htmp) {
        std::cerr << "[Warning] 1D histogram '" << histName << "' not found in " << fname << std::endl;
        f->Close();
        continue;
      }
      if (!h1) {
        h1 = (TH1D*)htmp->Clone(histName + "_merged");
        h1->SetDirectory(0);
      } else {
        h1->Add(htmp);
      }
      f->Close();
      std::cout << "[Info] Loaded 1D " << histName << " from " << fname << std::endl;
    }
    return h1;
  };

  // ----------------------------------------------------------
  // Helper: draw a 2D histogram (COLZ) and save
  // ----------------------------------------------------------
  auto draw2D = [&](TH2F* h2,
                    const TString& xTitle,
                    const TString& yTitle,
                    const TString& title,
                    double xLo, double xHi,
                    double yLo, double yHi,
                    const TString& canName,
                    const TString& outFile) {
    TCanvas* c = new TCanvas(canName, title, 900, 700);
    c->SetLeftMargin(0.11);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.14);
    c->SetTopMargin(0.06);
    c->SetGrid();

    h2->GetXaxis()->SetTitleOffset(1.3);
    h2->GetYaxis()->SetTitleOffset(1.2);

    h2->GetXaxis()->SetRangeUser(xLo, xHi);
    h2->GetYaxis()->SetRangeUser(yLo, yHi);
    h2->GetXaxis()->SetTitle(xTitle);
    h2->GetYaxis()->SetTitle(yTitle);
    h2->Draw("COLZ");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.034);
    tex.DrawLatex(0.55, 0.89, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.55, 0.84, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");
    tex.DrawLatex(0.55, 0.79, Form("#bf{Entries: %.0f}", h2->GetEntries()));

    c->SaveAs(outDir + outFile);
    std::cout << "[Info] Saved: " << outDir + outFile << std::endl;
    delete c;
  };

  // ==========================================================
  // (A) pT 2D Response Matrix  hResponseMatrixPairPt
  //     X: pT_MC  Y: pT_reco
  // ==========================================================
  {
    TH2F* h2 = merge2D("hResponseMatrixPairPt");
    if (h2) {
      draw2D(h2,
             "p_{T,#mu#mu}^{MC} (GeV/c)",
             "p_{T,#mu#mu}^{reco} (GeV/c)",
             "Pair p_{T} Response Matrix",
             0.0, 5.0, 0.0, 5.0,
             "c_rmat_pt",
             "Step1_ResponseMatrix_PairPt.png");
      delete h2;
    }
  }

  // ==========================================================
  // (B) pT^2 2D Response Matrix  hResponseMatrixPairPt2
  //     X: pT^2_MC  Y: pT^2_reco  (|t| ≈ pT^2)
  // ==========================================================
  TH2F* h2Pt2 = nullptr; // keep for saving
  {
    h2Pt2 = merge2D("hResponseMatrixPairPt2");
    if (h2Pt2) {
      draw2D(h2Pt2,
             "p_{T,#mu#mu}^{2,MC} (GeV^{2}/c^{2})",
             "p_{T,#mu#mu}^{2,reco} (GeV^{2}/c^{2})",
             "Pair p_{T}^{2} Response (Zoom 0-0.005)",
             0.0, 0.005, 0.0, 0.005,
             "c_rmat_pt2_zoom",
             "Step1_ResponseMatrix_PairPt2_Zoom.png");

      draw2D(h2Pt2,
             "p_{T,#mu#mu}^{2,MC} (GeV^{2}/c^{2})",
             "p_{T,#mu#mu}^{2,reco} (GeV^{2}/c^{2})",
             "Pair p_{T}^{2} Response (|t| #approx p_{T}^{2})",
             0.0, 2.5, 0.0, 2.5,
             "c_rmat_pt2",
             "Step1_ResponseMatrix_PairPt2.png");

      h2Pt2->GetXaxis()->UnZoom();
      h2Pt2->GetYaxis()->UnZoom();
    }
  }

  // ==========================================================
  // (C) pT Relative Residual  hPairPtResoVsPtMC
  //     X: pT_MC  Y: (pT_reco - pT_MC) / pT_MC
  // ==========================================================
  {
    TH2F* h2 = merge2D("hPairPtResoVsPtMC");
    if (h2) {
      draw2D(h2,
             "p_{T,#mu#mu}^{MC} (GeV/c)",
             "(p_{T}^{reco} - p_{T}^{MC}) / p_{T}^{MC}",
             "Pair p_{T} Resolution",
             0.0, 5.0, -1.0, 1.0,
             "c_reso_pt",
             "Step1_Resolution_PairPt.png");
      delete h2;
    }
  }

  // ==========================================================
  // (D) pT^2 Relative Residual  hPairPt2ResoVsPt2MC
  //     X: pT^2_MC  Y: (pT^2_reco - pT^2_MC) / pT^2_MC
  // ==========================================================
  {
    TH2F* h2 = merge2D("hPairPt2ResoVsPt2MC");
    if (h2) {
      draw2D(h2,
             "p_{T,#mu#mu}^{2,MC} (GeV^{2}/c^{2})",
             "(p_{T}^{2,reco} - p_{T}^{2,MC}) / p_{T}^{2,MC}",
             "Pair p_{T}^{2} Resolution  (|t| #approx p_{T}^{2})",
             0.0, 2.5, -5.0, 30.0,
             "c_reso_pt2",
             "Step1_Resolution_PairPt2.png");
      delete h2;
    }
  }

  // ==========================================================
  // Save merged response matrix AND 1D projections for downstream use
  // ==========================================================
  TFile* fOut = TFile::Open(outDir + "Step1_merged.root", "RECREATE");

  // pT^2 projection
  if (h2Pt2) {
    h2Pt2->Write("hResponseMatrixPairPt2");

    // MC pT^2 projection (X)
    TH1D* hGenPt2 = h2Pt2->ProjectionX("hGenPt2");
    hGenPt2->SetTitle("Gen pT^2 from 2D histogram");
    hGenPt2->GetXaxis()->SetTitle("p_{T}^{2,MC} (GeV^{2}/c^{2})");
    hGenPt2->GetYaxis()->SetTitle("Counts");
    hGenPt2->Write();
    delete hGenPt2;

    // Reco pT^2 projection (Y)
    TH1D* hRecoPt2 = h2Pt2->ProjectionY("hRecoPt2");
    hRecoPt2->SetTitle("Reco pT^2 from 2D histogram");
    hRecoPt2->GetXaxis()->SetTitle("p_{T}^{2,reco} (GeV^{2}/c^{2})");
    hRecoPt2->GetYaxis()->SetTitle("Counts");
    hRecoPt2->Write();
    delete hRecoPt2;

    std::cout << "[Info] Saved 2D matrix and 1D projections for Pair pT^2." << std::endl;

    delete h2Pt2;
  }

  fOut->Close();
  std::cout << "[Info] Saved merged histograms: " << outDir + "Step1_merged.root" << std::endl;

  std::cout << "[Done] Step 1 finished." << std::endl;
}
