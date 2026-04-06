#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TPad.h"

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
    // Create 1D Projections for Marginal Plots
    // ----------------------------------------------------------
    TH1D* hx = h2Rebinned->ProjectionX("hx");
    TH1D* hy = h2Rebinned->ProjectionY("hy");

    hx->SetFillColor(kAzure + 1);
    hx->SetLineColor(kBlack);
    hx->SetTitle("");
    hx->SetStats(0);

    hy->SetFillColor(kAzure + 1);
    hy->SetLineColor(kBlack);
    hy->SetTitle("");
    hy->SetStats(0);

    // ----------------------------------------------------------
    // Drawing the Marginal Plot
    // ----------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(256);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    const int fontCode = 42;

    TCanvas* c1 = new TCanvas("c1_step1", "Step1: Marginal Plot", 900, 900);

    // Define Pads (x_low, y_low, x_up, y_up)
    TPad* pad_center = new TPad("pad_center", "pad_center", 0.0,  0.0,  0.75, 0.75);
    TPad* pad_top    = new TPad("pad_top",    "pad_top",    0.0,  0.75, 0.75, 1.0);
    TPad* pad_right  = new TPad("pad_right",  "pad_right",  0.75, 0.0,  1.0,  0.75);

    // Set margins to eliminate gaps between pads
    pad_center->SetBottomMargin(0.12);
    pad_center->SetLeftMargin(0.12);
    pad_center->SetTopMargin(0.0);
    pad_center->SetRightMargin(0.0);
    pad_center->SetGrid();

    pad_top->SetBottomMargin(0.0);
    pad_top->SetLeftMargin(0.12);
    pad_top->SetTopMargin(0.1);
    pad_top->SetRightMargin(0.0);

    pad_right->SetBottomMargin(0.12);
    pad_right->SetLeftMargin(0.0);
    pad_right->SetTopMargin(0.0);
    pad_right->SetRightMargin(0.1);

    pad_center->Draw();
    pad_top->Draw();
    pad_right->Draw();

    // --- 1. Center Pad (2D Histogram) ---
    pad_center->cd();
    for (const char* ax : {"X", "Y"}) {
        h2Rebinned->GetXaxis()->SetLabelFont(fontCode);
        h2Rebinned->GetYaxis()->SetLabelFont(fontCode);
        h2Rebinned->GetXaxis()->SetTitleFont(fontCode);
        h2Rebinned->GetYaxis()->SetTitleFont(fontCode);
        h2Rebinned->GetXaxis()->SetTitleSize(0.045);
        h2Rebinned->GetYaxis()->SetTitleSize(0.045);
        h2Rebinned->GetXaxis()->SetLabelSize(0.04);
        h2Rebinned->GetYaxis()->SetLabelSize(0.04);
    }
    h2Rebinned->GetXaxis()->SetRangeUser(drawXmin, drawXmax);
    h2Rebinned->GetYaxis()->SetRangeUser(drawYmin, drawYmax);
    h2Rebinned->GetXaxis()->SetTitleOffset(1.2);
    h2Rebinned->GetYaxis()->SetTitleOffset(1.3);
    h2Rebinned->Draw("COL"); // Without "Z" to avoid color bar overlap

    // --- 2. Top Pad (1D Projection X) ---
    pad_top->cd();
    hx->GetXaxis()->SetRangeUser(drawXmin, drawXmax);
    hx->GetXaxis()->SetLabelOffset(999); // Hide X axis labels
    hx->GetXaxis()->SetTickLength(0.05);
    hx->GetYaxis()->SetLabelSize(0.08);  // Larger text due to smaller pad height
    hx->GetYaxis()->SetNdivisions(505);
    hx->Draw("HIST");

    // --- 3. Right Pad (1D Projection Y) ---
    pad_right->cd();
    hy->GetXaxis()->SetRangeUser(drawYmin, drawYmax);
    hy->GetXaxis()->SetLabelOffset(999); // Hide Y axis labels
    hy->GetXaxis()->SetTickLength(0.05);
    hy->GetYaxis()->SetLabelSize(0.08);  // Larger text due to smaller pad width
    hy->GetYaxis()->SetNdivisions(505);
    hy->Draw("hbar"); // Draw horizontally

    // --- 4. Main Canvas (Text in top-right corner) ---
    c1->cd(); 
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.025);
    tex.DrawLatex(0.77, 0.92, Form("#bf{p_{T} Resolution}"));
    tex.DrawLatex(0.77, 0.88, "#bf{LHC26b8} #font[52]{(MC)}");
    tex.DrawLatex(0.77, 0.84, Form("#bf{Entries: %.0f}", h2Raw->GetEntries()));

    c1->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_2D_Marginal.png");
    std::cout << "[Info] Saved image: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_2D_Marginal.png" << std::endl;

    // ----------------------------------------------------------
    // Save the merged histogram for future use (e.g., fitting in Step 2)
    // ----------------------------------------------------------
    TFile* fOut = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_merged.root", "RECREATE");
    h2Rebinned->Write("h2Rebinned");
    hx->Write("hx_proj");
    hy->Write("hy_proj");
    fOut->Close();
    std::cout << "[Info] Saved root file: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_merged.root" << std::endl;

    delete c1;
    delete h2Rebinned;
    delete h2Raw;

    std::cout << "[Done] Step 1 finished." << std::endl;
}