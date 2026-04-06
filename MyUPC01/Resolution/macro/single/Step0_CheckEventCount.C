// ============================================================
//  Step0_CheckEventCount.C
// ============================================================
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TColor.h"

void Step0_CheckEventCount()
{
    // ----------------------------------------------------------
    // Settings
    // ----------------------------------------------------------
    const TString dataDir  = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/";
    const TString histPath = "my-upc-muon-resolution/registry/hPtResoVsPtTrue_PostCut";

    struct Sample {
        TString file;
        TString label;
    };
    std::vector<Sample> samples = {
        { dataDir + "jpsi-incoh.root",    "J/#psi incoh." },
        { dataDir + "jpsi-coh.root",      "J/#psi coh." },
        { dataDir + "psi2s-incoh.root",   "#psi(2S) incoh." },
        { dataDir + "psi2s-coh.root",     "#psi(2S) coh." },
        { dataDir + "psi2s-incoh-fd.root","#psi(2S) incoh. fd" },
        { dataDir + "psi2s-coh-fd.root",  "#psi(2S) coh. fd" },
        { dataDir + "mumu-low.root",      "#mu#mu low" },
        { dataDir + "mumu-mid.root",      "#mu#mu mid" },
        { dataDir + "mumu-high.root",     "#mu#mu high" },
    };

    // X axis range
    const double xMin = 0.0, xMax = 6.0;

    // ----------------------------------------------------------
    // Styles
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

    // Colors for each sample
    int colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1,
                    kOrange+1, kCyan+2, kViolet+1, kTeal+2, kPink+1};

    // ----------------------------------------------------------
    // Figure1: Total pT distribution
    // ----------------------------------------------------------
    TH1D* hTotal = nullptr;

    for (size_t is = 0; is < samples.size(); ++is) {
        TFile* f = TFile::Open(samples[is].file, "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[Warning] File not found: " << samples[is].file << std::endl;
            if (f) f->Close();
            continue;
        }
        TH2F* h2 = (TH2F*)f->Get(histPath);
        if (!h2) {
            std::cerr << "[Warning] Histogram not found in: " << samples[is].file << std::endl;
            f->Close();
            continue;
        }

        // pT true distribution
        TH1D* hProjX = h2->ProjectionX(Form("hProjX_%zu", is), 1, h2->GetNbinsY());

        if (!hTotal) {
            hTotal = (TH1D*)hProjX->Clone("hTotal");
            hTotal->SetDirectory(0);
        } else {
            hTotal->Add(hProjX);
        }

        delete hProjX;
        f->Close();
        std::cout << "[Info] Loaded: " << samples[is].file << std::endl;
    }

    if (!hTotal) {
        std::cerr << "[Error] No data loaded. Abort." << std::endl;
        return;
    }

    TCanvas* c1 = new TCanvas("c1_step0", "Step0: pT True Distribution (Total)", 900, 650);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.13);
    c1->SetTopMargin(0.05);
    c1->SetRightMargin(0.05);
    c1->SetLogy();  // Log y axis
    c1->SetGrid();

    hTotal->GetXaxis()->SetRangeUser(xMin, xMax);
    hTotal->GetYaxis()->SetRangeUser(1e0, hTotal->GetMaximum()*10);
    hTotal->GetXaxis()->SetTitle("p_{T}^{true} (GeV/c)");
    hTotal->GetYaxis()->SetTitle("Counts");
    hTotal->SetLineColor(kBlack);
    hTotal->SetLineWidth(2);
    hTotal->SetFillColorAlpha(kAzure+1, 0.30);
    hTotal->SetFillStyle(3004);
    hTotal->Draw("HIST");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.45, 0.90, "#bf{p_{T}^{true} distribution after Coherent J/#psi cut}");
    tex.DrawLatex(0.45, 0.85, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");
    tex.DrawLatex(0.45, 0.80, Form("Total entries: %.0f", hTotal->GetEntries()));

    c1->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step0_EventCount_pT.png");
    std::cout << "[Info] Saved: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step0_EventCount_pT.png" << std::endl;

    // ----------------------------------------------------------
    // Figure 2: pT true distribution (per sample)
    // ----------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2_step0", "Step0: pT True Distribution (per sample)", 1000, 650);
    c2->SetLeftMargin(0.12);
    c2->SetBottomMargin(0.13);
    c2->SetTopMargin(0.05);
    c2->SetRightMargin(0.05);
    c2->SetLogy();
    c2->SetGrid();

    TLegend* leg = new TLegend(0.66, 0.55, 0.94, 0.94);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(fontCode);
    leg->SetTextSize(0.032);

    bool firstDraw = true;

    for (size_t is = 0; is < samples.size(); ++is) {
        TFile* f = TFile::Open(samples[is].file, "READ");
        if (!f || f->IsZombie()) { if (f) f->Close(); continue; }
        TH2F* h2 = (TH2F*)f->Get(histPath);
        if (!h2) { f->Close(); continue; }

        TH1D* hProj = h2->ProjectionX(Form("hProjX_each_%zu", is), 1, h2->GetNbinsY());
        hProj->SetDirectory(0);
        f->Close();

        hProj->GetXaxis()->SetRangeUser(xMin, xMax);
        hProj->GetYaxis()->SetRangeUser(1e0, hProj->GetMaximum()*100);
        hProj->GetXaxis()->SetTitle("p_{T}^{true} (GeV/c)");
        hProj->GetYaxis()->SetTitle("Counts");
        hProj->SetLineColor(colors[is % 9]);
        hProj->SetLineWidth(2);
        hProj->SetMarkerColor(colors[is % 9]);
        hProj->SetMarkerStyle(20 + (int)is);
        hProj->SetMarkerSize(0.8);

        if (firstDraw) {
            hProj->Draw("HIST");
            firstDraw = false;
        } else {
            hProj->Draw("HIST SAME");
        }

        leg->AddEntry(hProj, Form("%s  (%.0f)", samples[is].label.Data(), hProj->GetEntries()), "l");
    }

    leg->Draw();

    TLatex tex2;
    tex2.SetNDC();
    tex2.SetTextFont(fontCode);
    tex2.SetTextSize(0.038);
    tex2.DrawLatex(0.14, 0.90, "#bf{p_{T}^{true} distribution after Coherent J/#psi cut}");
    tex2.DrawLatex(0.14, 0.85, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");

    c2->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step0_EventCount_pT_each.png");
    std::cout << "[Info] Saved: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step0_EventCount_pT_each.png" << std::endl;

    // ----------------------------------------------------------
    // [コンソール出力] ビンごとのカウント数を表示
    // ----------------------------------------------------------
    std::cout << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "  pT_true bin contents (total, merged)  " << std::endl;
    std::cout << "  (X-axis bin center / width / counts)  " << std::endl;
    std::cout << "=========================================" << std::endl;
    for (int ib = 1; ib <= hTotal->GetNbinsX(); ++ib) {
        double center  = hTotal->GetBinCenter(ib);
        double width   = hTotal->GetBinWidth(ib);
        double content = hTotal->GetBinContent(ib);
        if (center > xMax) break;
        std::cout << Form("  pT = [%4.2f, %4.2f] GeV/c  ( center=%4.2f )  counts = %.0f",
                          center - width/2.0, center + width/2.0, center, content) << std::endl;
    }
    std::cout << "=========================================" << std::endl;

    delete leg;
    delete c1;
    delete c2;
    delete hTotal;

    std::cout << "[Done] Step 0 finished." << std::endl;
}
