#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TDirectory.h"

void plot_2d_correlations()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(57); // kRainBow

    TFile* f = TFile::Open("AnalysisResults.root");
    if (!f || f->IsZombie()) {
        f = TFile::Open("Analysis.root");
        if (!f || f->IsZombie()) {
            printf("Error: Cannot open AnalysisResults.root or Analysis.root\n");
            return;
        }
    }

    TDirectory* dir = (TDirectory*)f->Get("fwd_dca/PDCA");
    if (!dir) {
        printf("Error: Could not find 'fwd_dca/PDCA' directory.\n");
        f->Close();
        return;
    }

    TH2F* h2_beauty = (TH2F*)dir->Get("pdca_vs_pt_beauty");
    TH2F* h2_pion   = (TH2F*)dir->Get("pdca_vs_pt_pion");
    TH2F* h2_charm  = (TH2F*)dir->Get("pdca_vs_pt_charm");

    if (h2_beauty) {
        TCanvas* c1 = new TCanvas("c_2d_beauty", "pDCA vs pT (Beauty)", 800, 600);
        c1->SetLogz(); // Z軸を対数スケール
        h2_beauty->Draw("COLZ");
        c1->SaveAs("pdca_vs_pt_beauty.png");
        delete c1;
    }

    if (h2_pion) {
        TCanvas* c2 = new TCanvas("c_2d_pion", "pDCA vs pT (Pion)", 800, 600);
        c2->SetLogz();
        h2_pion->Draw("COLZ");
        c2->SaveAs("pdca_vs_pt_pion.png");
        delete c2;
    }

    if (h2_charm) {
        TCanvas* c3 = new TCanvas("c_2d_charm", "pDCA vs pT (Charm)", 800, 600);
        c3->SetLogz();
        h2_charm->Draw("COLZ");
        c3->SaveAs("pdca_vs_pt_charm.png");
        delete c3;
    }

    f->Close();
}