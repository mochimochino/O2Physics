#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"

void Step1_MergeAndView2D()
{
    // ----------------------------------------------------------
    // Input and path settings
    // ----------------------------------------------------------
    const TString dataDir  = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/";
    const TString histPath = "my-upc-muon-resolution/registry/hPtResoVsPtTrue_PostCut";

    std::vector<TString> fileNames = {
        dataDir + "jpsi-incoh.root",
        dataDir + "jpsi-coh.root",
        dataDir + "psi2s-incoh.root",
        dataDir + "psi2s-coh.root",
        dataDir + "psi2s-incoh-fd.root",
        dataDir + "psi2s-coh-fd.root",
        dataDir + "mumu-low.root",
        dataDir + "mumu-mid.root",
        dataDir + "mumu-high.root"
    };

double ptBins[] = {0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                       1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 
                       1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 
                       4.0, 5.0, 6.0};
    int nBinsX = sizeof(ptBins) / sizeof(double) - 1;

    double drawXmin = 0.0, drawXmax = 6.0;
    double drawYmin = -0.5, drawYmax = 0.5;

    // ----------------------------------------------------------
    // merge histograms from multiple files
    // ----------------------------------------------------------
    TH2F* h2Raw = nullptr;

    for (const auto& fname : fileNames) {
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[Warning] File not found or broken: " << fname << std::endl;
            if (f) f->Close();
            continue;
        }
        TH2F* h_tmp = (TH2F*)f->Get(histPath);
        if (!h_tmp) {
            std::cerr << "[Warning] Histogram not found in: " << fname << std::endl;
            f->Close();
            continue;
        }
        if (!h2Raw) {
            h2Raw = (TH2F*)h_tmp->Clone("h2Raw_total");
            h2Raw->SetDirectory(0);
        } else {
            h2Raw->Add(h_tmp);
        }
        f->Close();
        std::cout << "[Info] Loaded: " << fname << std::endl;
    }

    if (!h2Raw) {
        std::cerr << "[Error] No valid histograms could be loaded. Abort." << std::endl;
        return;
    }
    std::cout << "[Info] Total entries (merged): " << h2Raw->GetEntries() << std::endl;

    // ----------------------------------------------------------
    // Rebinning pT axis with custom bins
    // ----------------------------------------------------------
    int    nBinsY = h2Raw->GetNbinsY();
    double yMin   = h2Raw->GetYaxis()->GetXmin();
    double yMax   = h2Raw->GetYaxis()->GetXmax();

    TH2D* h2Rebinned = new TH2D("h2Rebinned",
        ";p_{T}^{true} (GeV/c);Resolution (p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}",
        nBinsX, ptBins, nBinsY, yMin, yMax);

    for (int ix = 1; ix <= h2Raw->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h2Raw->GetNbinsY(); ++iy) {
            double content = h2Raw->GetBinContent(ix, iy);
            if (content > 0) {
                double xc = h2Raw->GetXaxis()->GetBinCenter(ix);
                double yc = h2Raw->GetYaxis()->GetBinCenter(iy);
                h2Rebinned->Fill(xc, yc, content);
            }
        }
    }

    // ----------------------------------------------------------
    // Drawing the rebinned histogram
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
        gStyle->SetLabelSize(0.04,  ax);
    }

    TCanvas* c1 = new TCanvas("c1_step1", "Step1: Merged 2D Resolution", 900, 700);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.14);
    c1->SetTopMargin(0.05);
    c1->SetGrid();

    h2Rebinned->GetXaxis()->SetRangeUser(drawXmin, drawXmax);
    h2Rebinned->GetYaxis()->SetRangeUser(drawYmin, drawYmax);
    h2Rebinned->Draw("COLZ");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.035);
    tex.DrawLatex(0.40, 0.90, Form("#bf{p_{T} Resolution After Coherent J/#psi Cut}"));
    tex.DrawLatex(0.40, 0.85, "#bf{LHC26b8} #font[52]{(MC dataset)}");
    tex.DrawLatex(0.40, 0.80, Form("#bf{Entries: %.0f}", h2Raw->GetEntries()));

    c1->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_2D_Merged.png");
    std::cout << "[Info] Saved: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_2D_Merged.png" << std::endl;

    // ----------------------------------------------------------
    // Save the merged histogram for future use (e.g., fitting in Step 2)
    // ----------------------------------------------------------
    TFile* fOut = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_merged.root", "RECREATE");
    h2Rebinned->Write("h2Rebinned");
    fOut->Close();
    std::cout << "[Info] Saved merged histogram: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_merged.root" << std::endl;

    delete c1;
    delete h2Rebinned;
    delete h2Raw;

    std::cout << "[Done] Step 1 finished." << std::endl;
}
