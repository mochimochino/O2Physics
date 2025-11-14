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

void muon_mother()
{
    const char* taskname = "mc_muon_resto_track"; 
    
    gStyle->SetOptStat(0); 

    TFile* f = TFile::Open("AnalysisResults.root");
    if (!f || f->IsZombie()) {
        printf("Error: Cannot open AnalysisResults.root\n");
        return;
    }
    TH1F* hPt_Reco_allgen  = (TH1F*)f->Get(Form("%s/pt_reco", taskname));
    TH1F* hPt_Reco_Primary = (TH1F*)f->Get(Form("%s/pt_reco_primary", taskname));
    TH1F* hPt_Reco_Pion    = (TH1F*)f->Get(Form("%s/pt_reco_pion", taskname));
    TH1F* hPt_Reco_Kaon    = (TH1F*)f->Get(Form("%s/pt_reco_kaon", taskname));
    TH1F* hPt_Reco_Charm   = (TH1F*)f->Get(Form("%s/pt_reco_charm", taskname));
    TH1F* hPt_Reco_Beauty  = (TH1F*)f->Get(Form("%s/pt_reco_beauty", taskname));
    TH1F* hPt_Reco_Other   = (TH1F*)f->Get(Form("%s/pt_reco_other", taskname));

    if (!hPt_Reco_allgen || !hPt_Reco_Primary || !hPt_Reco_Pion || !hPt_Reco_Kaon || !hPt_Reco_Charm || !hPt_Reco_Beauty || !hPt_Reco_Other) {
         printf("Error: Could not retrieve all required histograms from directory %s\n", taskname);
         f->Close();
         return;
    }

    {
        TH1F* hPt_Reco_Total = (TH1F*)hPt_Reco_Primary->Clone("hPt_Reco_Total_Frac");
        hPt_Reco_Total->Add(hPt_Reco_Pion);
        hPt_Reco_Total->Add(hPt_Reco_Kaon);
        hPt_Reco_Total->Add(hPt_Reco_Charm);
        hPt_Reco_Total->Add(hPt_Reco_Beauty);
        hPt_Reco_Total->Add(hPt_Reco_Other);

        TH1F* hFrac_Primary = (TH1F*)hPt_Reco_Primary->Clone("hFrac_Primary");
        TH1F* hFrac_Pion    = (TH1F*)hPt_Reco_Pion->Clone("hFrac_Pion");
        TH1F* hFrac_Kaon    = (TH1F*)hPt_Reco_Kaon->Clone("hFrac_Kaon");
        TH1F* hFrac_Charm   = (TH1F*)hPt_Reco_Charm->Clone("hFrac_Charm");
        TH1F* hFrac_Beauty  = (TH1F*)hPt_Reco_Beauty->Clone("hFrac_Beauty");
        TH1F* hFrac_Other   = (TH1F*)hPt_Reco_Other->Clone("hFrac_Other");

        hFrac_Primary->Divide(hPt_Reco_Total);
        hFrac_Pion->Divide(hPt_Reco_Total);
        hFrac_Kaon->Divide(hPt_Reco_Total);
        hFrac_Charm->Divide(hPt_Reco_Total);
        hFrac_Beauty->Divide(hPt_Reco_Total);
        hFrac_Other->Divide(hPt_Reco_Total);

        THStack* hsFraction = new THStack("hsFraction", "Muon Parent Composition (Fraction);p_{T}^{reco} [GeV/c];Fraction");

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

        TCanvas* c1 = new TCanvas("c1_frac", "Muon Parent Composition (Fraction)", 800, 600);
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

        c1->SaveAs("muon_composition_stack_fraction.png");

        delete c1;
        delete hsFraction;
        delete hPt_Reco_Total;
    }


    { 
        hPt_Reco_allgen->SetLineColor(TColor::GetColor("#333333")); // ダークグレー
        hPt_Reco_Beauty->SetLineColor(TColor::GetColor("#E15759")); // 落ち着いた赤
        hPt_Reco_Charm->SetLineColor(TColor::GetColor("#F28E2B"));  // オレンジ
        hPt_Reco_Kaon->SetLineColor(TColor::GetColor("#59A14F"));   // 緑
        hPt_Reco_Pion->SetLineColor(TColor::GetColor("#4E79A7"));   // 青
        hPt_Reco_Primary->SetLineColor(TColor::GetColor("#76B7B2")); // 青緑 (シアン)
        hPt_Reco_Other->SetLineColor(TColor::GetColor("#B07AA1"));  // 紫

        hPt_Reco_Beauty->SetLineStyle(2); // 破線
        hPt_Reco_Charm->SetLineStyle(2);  // 破線
        hPt_Reco_Kaon->SetLineStyle(2);   // 破線
        hPt_Reco_Pion->SetLineStyle(2);   // 破線
        hPt_Reco_Primary->SetLineStyle(2); // 破線
        hPt_Reco_Other->SetLineStyle(2);  // 破線

        hPt_Reco_allgen->SetLineWidth(2);
        hPt_Reco_Beauty->SetLineWidth(2);
        hPt_Reco_Charm->SetLineWidth(2);
        hPt_Reco_Kaon->SetLineWidth(2);
        hPt_Reco_Pion->SetLineWidth(2);
        hPt_Reco_Primary->SetLineWidth(2);
        hPt_Reco_Other->SetLineWidth(2);

        double yMax = 0;
        yMax = TMath::Max(yMax, hPt_Reco_allgen->GetMaximum());
        yMax = TMath::Max(yMax, hPt_Reco_Primary->GetMaximum());
        yMax = TMath::Max(yMax, hPt_Reco_Pion->GetMaximum());
        yMax = TMath::Max(yMax, hPt_Reco_Kaon->GetMaximum());
        yMax = TMath::Max(yMax, hPt_Reco_Charm->GetMaximum());
        yMax = TMath::Max(yMax, hPt_Reco_Beauty->GetMaximum());
        yMax = TMath::Max(yMax, hPt_Reco_Other->GetMaximum());

        TCanvas* c2 = new TCanvas("c2_counts", "Muon Parent Composition (Comparison)", 800, 600);
        c2->cd();
        c2->SetLogy();
        gPad->SetGrid();

        hPt_Reco_Pion->SetTitle("Muon Parent Composition (Comparison);p_{T}^{reco} [GeV/c];Counts");
        hPt_Reco_Pion->GetYaxis()->SetRangeUser(0.5, yMax * 5.0);
        hPt_Reco_Pion->Draw("HIST");
        hPt_Reco_allgen->Draw("HIST SAME");
        hPt_Reco_Kaon->Draw("HIST SAME");
        hPt_Reco_Primary->Draw("HIST SAME");
        hPt_Reco_Charm->Draw("HIST SAME");
        hPt_Reco_Beauty->Draw("HIST SAME");
        hPt_Reco_Other->Draw("HIST SAME");

        TLegend* leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg2->SetBorderSize(0);
        leg2->AddEntry(hPt_Reco_allgen, "All generated reco muons", "l");
        leg2->AddEntry(hPt_Reco_Beauty, "Beauty decay", "l");
        leg2->AddEntry(hPt_Reco_Charm, "Charm decay", "l");
        leg2->AddEntry(hPt_Reco_Pion, "Pion decay", "l");
        leg2->AddEntry(hPt_Reco_Kaon, "Kaon decay", "l");
        leg2->AddEntry(hPt_Reco_Primary, "Primary", "l");
        leg2->AddEntry(hPt_Reco_Other, "Other", "l");
        leg2->Draw();

        c2->SaveAs("muon_composition_compare.png");
        
        delete c2;
    }

    f->Close();
}