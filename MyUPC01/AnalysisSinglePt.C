#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include <iostream>

void AnalysisSinglePt() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.04);
    
    // Using the same file path as in AnalysisMassReal.C
    TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/0305withoutpTtask.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open AnalysisResults.root" << std::endl;
        return;
    }

    TString dir = "my-upc-01/";
    
    TH1* hPt1 = dynamic_cast<TH1*>(f->Get(dir + "ptMuon1"));
    TH1* hPt2 = dynamic_cast<TH1*>(f->Get(dir + "ptMuon2"));

    if (!hPt1 || !hPt2) {
        std::cerr << "Error: Histograms ptMuon1 or ptMuon2 not found in " << dir << std::endl;
        f->ls();
        return;
    }

    // ==========================================
    // 1. Single Muon Pt Histograms
    // ==========================================
    TCanvas* c1 = new TCanvas("c1", " DiMuon Pt", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.12);
    
    gPad->SetLogy(); // y-axis to log scale
    gPad->SetGrid(); // show grid

    // Define colors from AnalysisMassReal.C
    int col1 = TColor::GetColor("#2e8b57"); // SeaGreen
    int col2 = TColor::GetColor("#d2691e"); // Chocolate

    // Muon 1 Pt settings
    hPt1->SetTitle("DiMuon p_{T}; p_{T} [GeV/c]; Counts");
    hPt1->SetLineColor(col1);
    hPt1->SetMarkerColor(col1);
    hPt1->SetLineWidth(2);
    hPt1->SetMarkerStyle(21); // Full Square
    hPt1->GetYaxis()->SetTitleOffset(1.2);
    
    // Set range to accommodate both plots
    double max = std::max(hPt1->GetMaximum(), hPt2->GetMaximum());
    hPt1->GetYaxis()->SetRangeUser(0.8, max * 5.0);

    hPt1->GetXaxis()->SetRangeUser(0.0, 5.0);
    
    hPt1->Draw("E");

    // Muon 2 Pt settings
    hPt2->SetLineColor(col2);
    hPt2->SetMarkerColor(col2);
    hPt2->SetLineWidth(2);
    hPt2->SetMarkerStyle(25); // Open Square
    hPt2->Draw("E SAME");

    // Legend settings
    TLegend* leg = new TLegend(0.65, 0.80, 0.95, 0.95);
    leg->SetBorderSize(1);
    leg->AddEntry(hPt1, "Muon 1", "lep");
    leg->AddEntry(hPt2, "Muon 2", "lep");
    leg->Draw();

    c1->SaveAs("DiMuonPt.png");

    std::cout << "Done. Saved plot to DiMuonPt.png" << std::endl;

    // f->Close(); // Keep open if running interactively, or close if batch
}
