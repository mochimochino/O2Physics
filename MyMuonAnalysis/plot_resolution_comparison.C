#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TColor.h"
#include "TF1.h" 

void compareResolution(TDirectory* dir, const char* histName1, const char* title1, const char* histName2, const char* title2, 
                       const char* canvasName, const char* canvasTitle, const char* xAxisTitle, const char* saveName)
{
    TH1F* h1 = (TH1F*)dir->Get(histName1);
    TH1F* h2 = (TH1F*)dir->Get(histName2);

    if (!h1 || !h2) {
        printf("Error: Cannot get histograms '%s' or '%s'\n", histName1, histName2);
        return;
    }

    TCanvas* c = new TCanvas(canvasName, canvasTitle, 800, 600);
    gPad->SetGrid();

    // 規格化 
    //double integral1 = h1->Integral(0, -1);
    //double integral2 = h2->Integral(0, -1);
    //if (integral1 > 0) h1->Scale(1.0 / integral1);
    //if (integral2 > 0) h2->Scale(1.0 / integral2);


    h1->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    h1->SetTitle(Form("%s;%s;Normalized Entries", canvasTitle, xAxisTitle));

    h2->SetLineColor(kRed);
    h2->SetLineWidth(2);
    h2->SetLineStyle(2); 

    float max = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    h1->GetYaxis()->SetRangeUser(0, max * 1.2);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.75, 0.75, 0.90, 0.90);
    leg->SetBorderSize(0); 
    leg->AddEntry(h1, title1, "l");
    leg->AddEntry(h2, title2, "l");
    leg->Draw();

    c->SaveAs(saveName);
    delete c;
}

void plot_resolution_comparison()
{
    gStyle->SetOptStat(0);
    // gStyle->SetOptFit(0); // Fitを使わないので不要

    TFile* f = TFile::Open("AnalysisResults.root");
    if (!f || f->IsZombie()) {
        f = TFile::Open("Analysis.root");
        if (!f || f->IsZombie()) {
            printf("Error: Cannot open AnalysisResults.root or Analysis.root\n");
            return;
        }
    }

    TDirectory* dir = (TDirectory*)f->Get("fwd_dca");
    if (!dir) {
        printf("Error: Could not find 'fwd_dca' directory.\n");
        f->Close();
        return;
    }

    // pT 
    compareResolution(dir, 
                   "TrackType/Type0/pt_resolution_Type0", "Type 0 (Global MFT+MCH+MID)", 
                   "TrackType/Type3/pt_resolution_Type3", "Type 3 (Standalone MCH+MID)",
                   "cPtRes", "pT Resolution Comparison", "(pT_reco - pT_gen) / pT_gen", 
                   "resolution_pt_compare.png");

    // eta
    compareResolution(dir, 
                   "TrackType/Type0/eta_resolution_Type0", "Type 0 (Global MFT+MCH+MID)", 
                   "TrackType/Type3/eta_resolution_Type3", "Type 3 (Standalone MCH+MID)",
                   "cEtaRes", "Eta Resolution Comparison", "#eta_reco - #eta_gen", 
                   "resolution_eta_compare.png");
    
    // phi
    compareResolution(dir, 
                   "TrackType/Type0/phi_resolution_Type0", "Type 0 (Global MFT+MCH+MID)", 
                   "TrackType/Type3/phi_resolution_Type3", "Type 3 (Standalone MCH+MID)",
                   "cPhiRes", "Phi Resolution Comparison", "#phi_reco - #phi_gen", 
                   "resolution_phi_compare.png");

    f->Close();
}