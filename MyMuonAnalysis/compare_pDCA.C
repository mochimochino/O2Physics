#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"       
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TString.h"
#include "TDirectory.h"
#include "TColor.h"
#include "TMath.h"

void compare_pDCA()
{
    const char* taskname = "fwd_dca"; 
    
    gStyle->SetOptStat(0); 

    TFile* f = TFile::Open("AnalysisResults.root");
    if (!f || f->IsZombie()) {
        if (!f || f->IsZombie()) {
            printf("Error: Cannot open AnalysisResults.root\n");
            return;
        }
    }
    
    TH1F* hPDca_allgen  = (TH1F*)f->Get(Form("%s/PDCA/pdca_all", taskname));
    TH1F* hPDca_Primary = (TH1F*)f->Get(Form("%s/PDCA/pdca_primary", taskname));
    TH1F* hPDca_Pion    = (TH1F*)f->Get(Form("%s/PDCA/pdca_pion", taskname));
    TH1F* hPDca_Kaon    = (TH1F*)f->Get(Form("%s/PDCA/pdca_kaon", taskname));
    TH1F* hPDca_Charm   = (TH1F*)f->Get(Form("%s/PDCA/pdca_charm", taskname));
    TH1F* hPDca_Beauty  = (TH1F*)f->Get(Form("%s/PDCA/pdca_beauty", taskname));
    TH1F* hPDca_Other   = (TH1F*)f->Get(Form("%s/PDCA/pdca_other", taskname));

    if (!hPDca_allgen || !hPDca_Primary || !hPDca_Pion || !hPDca_Kaon || !hPDca_Charm || !hPDca_Beauty || !hPDca_Other) {
         printf("Error: Could not retrieve all required histograms from directory %s/PDCA/\n", taskname);
         f->Close();
         return;
    }

    // --- 1. Fraction Stack Plot ---
    {
        TH1F* hPDca_Total = (TH1F*)hPDca_Primary->Clone("hPDca_Total_Frac");
        hPDca_Total->Add(hPDca_Pion);
        hPDca_Total->Add(hPDca_Kaon);
        hPDca_Total->Add(hPDca_Charm);
        hPDca_Total->Add(hPDca_Beauty);
        hPDca_Total->Add(hPDca_Other);

        TH1F* hFrac_Primary = (TH1F*)hPDca_Primary->Clone("hFrac_Primary");
        TH1F* hFrac_Pion    = (TH1F*)hPDca_Pion->Clone("hFrac_Pion");
        TH1F* hFrac_Kaon    = (TH1F*)hPDca_Kaon->Clone("hFrac_Kaon");
        TH1F* hFrac_Charm   = (TH1F*)hPDca_Charm->Clone("hFrac_Charm");
        TH1F* hFrac_Beauty  = (TH1F*)hPDca_Beauty->Clone("hFrac_Beauty");
        TH1F* hFrac_Other   = (TH1F*)hPDca_Other->Clone("hFrac_Other");

        hFrac_Primary->Divide(hPDca_Total);
        hFrac_Pion->Divide(hPDca_Total);
        hFrac_Kaon->Divide(hPDca_Total);
        hFrac_Charm->Divide(hPDca_Total);
        hFrac_Beauty->Divide(hPDca_Total);
        hFrac_Other->Divide(hPDca_Total);

        THStack* hsFraction = new THStack("hsFraction", "Muon Parent Composition (Fraction);p*DCA [GeV/c * cm];Fraction");

        hFrac_Beauty->SetFillColor(TColor::GetColor("#E15759")); // 落ち着いた赤
        hFrac_Charm->SetFillColor(TColor::GetColor("#F28E2B"));  // オレンジ
        hFrac_Kaon->SetFillColor(TColor::GetColor("#59A14F"));   // 緑
        hFrac_Pion->SetFillColor(TColor::GetColor("#4E79A7"));   // 青
        hFrac_Primary->SetFillColor(TColor::GetColor("#76B7B2")); // 青緑 (シアン)
        hFrac_Other->SetFillColor(TColor::GetColor("#B07AA1"));  // 紫

        hsFraction->Add(hFrac_Other);
        hsFraction->Add(hFrac_Primary);
        hsFraction->Add(hFrac_Kaon);
        hsFraction->Add(hFrac_Pion);
        hsFraction->Add(hFrac_Charm);
        hsFraction->Add(hFrac_Beauty); 

        TCanvas* c1 = new TCanvas("c1_frac_pdca", "Muon Parent Composition (pDCA Fraction)", 800, 600);
        c1->cd();
        gPad->SetGrid();

        hsFraction->Draw("HIST");
        
        hsFraction->GetYaxis()->SetRangeUser(0.0, 1.05);
        
        TLegend* leg1 = new TLegend(0.65, 0.65, 0.88, 0.88);
        leg1->SetBorderSize(0);
        leg1->AddEntry(hFrac_Beauty, "Beauty decay", "f");
        leg1->AddEntry(hFrac_Charm, "Charm decay", "f");
        leg1->AddEntry(hFrac_Pion, "Pion decay", "f");
        leg1->AddEntry(hFrac_Kaon, "Kaon decay", "f");
        leg1->AddEntry(hFrac_Primary, "Primary", "f");
        leg1->AddEntry(hFrac_Other, "Other", "f");
        leg1->Draw();

        //c1->SaveAs("muon_pdca_stack_fraction.png");

        delete c1;
        delete hsFraction;
        delete hPDca_Total;
    }

    // --- 2. Comparison Plot ---
    { 
        hPDca_allgen->SetLineColor(TColor::GetColor("#333333")); // ダークグレー
        hPDca_Beauty->SetLineColor(TColor::GetColor("#E15759")); // 落ち着いた赤
        hPDca_Charm->SetLineColor(TColor::GetColor("#F28E2B"));  // オレンジ
        hPDca_Kaon->SetLineColor(TColor::GetColor("#59A14F"));   // 緑
        hPDca_Pion->SetLineColor(TColor::GetColor("#4E79A7"));   // 青
        hPDca_Primary->SetLineColor(TColor::GetColor("#76B7B2")); // 青緑 (シアン)
        hPDca_Other->SetLineColor(TColor::GetColor("#B07AA1"));  // 紫

        hPDca_Beauty->SetLineStyle(2);
        hPDca_Charm->SetLineStyle(2);
        hPDca_Kaon->SetLineStyle(2);
        hPDca_Pion->SetLineStyle(2);
        hPDca_Primary->SetLineStyle(2);
        hPDca_Other->SetLineStyle(2);

        hPDca_allgen->SetLineWidth(2);
        hPDca_Beauty->SetLineWidth(2);
        hPDca_Charm->SetLineWidth(2);
        hPDca_Kaon->SetLineWidth(2);
        hPDca_Pion->SetLineWidth(2);
        hPDca_Primary->SetLineWidth(2);
        hPDca_Other->SetLineWidth(2);

        double yMax = 0;
        yMax = TMath::Max(yMax, hPDca_allgen->GetMaximum());
        yMax = TMath::Max(yMax, hPDca_Primary->GetMaximum());
        yMax = TMath::Max(yMax, hPDca_Pion->GetMaximum());
        yMax = TMath::Max(yMax, hPDca_Kaon->GetMaximum());
        yMax = TMath::Max(yMax, hPDca_Charm->GetMaximum());
        yMax = TMath::Max(yMax, hPDca_Beauty->GetMaximum());
        yMax = TMath::Max(yMax, hPDca_Other->GetMaximum());

        TCanvas* c2 = new TCanvas("c2_counts_pdca", "Muon Parent Composition (pDCA Comparison)", 800, 600);
        c2->cd();
        c2->SetLogy();
        // c2->SetLogx();
        gPad->SetGrid();

        hPDca_Pion->SetTitle("Muon Parent Composition (pDCA);p*DCA [GeV/c * cm];Counts");
        hPDca_Pion->GetYaxis()->SetRangeUser(0.5, yMax * 5.0);
        hPDca_Pion->Draw("HIST");
        hPDca_allgen->Draw("HIST SAME");
        hPDca_Kaon->Draw("HIST SAME");
        hPDca_Primary->Draw("HIST SAME");
        hPDca_Charm->Draw("HIST SAME");
        hPDca_Beauty->Draw("HIST SAME");
        hPDca_Other->Draw("HIST SAME");

        TLegend* leg2 = new TLegend(0.70, 0.70, 0.90, 0.90);
        leg2->SetBorderSize(0);
        leg2->AddEntry(hPDca_allgen, "All matched muons", "l");
        leg2->AddEntry(hPDca_Beauty, "Beauty decay", "l");
        leg2->AddEntry(hPDca_Charm, "Charm decay", "l");
        leg2->AddEntry(hPDca_Pion, "Pion decay", "l");
        leg2->AddEntry(hPDca_Kaon, "Kaon decay", "l");
        leg2->AddEntry(hPDca_Primary, "Primary", "l");
        leg2->AddEntry(hPDca_Other, "Other", "l");
        leg2->Draw();

        c2->SaveAs("muon_pdca_compare.pdf");
        
        delete c2;
    }

    f->Close();
}