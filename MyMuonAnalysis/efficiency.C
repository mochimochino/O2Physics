#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TAxis.h"
// my_muon_tracking_efficiency

void efficiency()
{
    TFile* f = TFile::Open("AnalysisResults.root");
    if (!f || f->IsZombie()) {
        printf("Error: Cannot open AnalysisResults.root\n");
        return;
    }

    TH1F* hDenPt = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Gen");

    TH1F* hNumPt_Type0 = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Reco_Matched_Type0");
    TH1F* hNumPt_Type3 = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Reco_Matched_Type3");

    if (!hDenPt || !hNumPt_Type0 || !hNumPt_Type3) {
        printf("Error: Cannot find histograms in file.\n");
        if (!hDenPt) printf("... hPt_Gen not found.\n");
        if (!hNumPt_Type0) printf("... hPt_Reco_Matched_Type0 not found.\n");
        if (!hNumPt_Type3) printf("... hPt_Reco_Matched_Type3 not found.\n");
        return;
    }

    TEfficiency* effPt_Type0 = new TEfficiency(*hNumPt_Type0, *hDenPt);
    effPt_Type0->SetName("effPt_Type0");

    TEfficiency* effPt_Type3 = new TEfficiency(*hNumPt_Type3, *hDenPt);
    effPt_Type3->SetName("effPt_Type3");

    TCanvas* c1 = new TCanvas("c1", "Tracking Efficiency vs pT", 800, 600);
    c1->cd();
    gPad->SetGrid(); 

    effPt_Type0->SetMarkerStyle(21);
    effPt_Type0->SetMarkerColor(kRed);
    effPt_Type0->SetLineColor(kRed);
    effPt_Type0->SetTitle("Muon Tracking Efficiency vs p_{T};p_{T} [GeV/c];Efficiency #epsilon");
    effPt_Type0->Draw("AP");

    gPad->Update();
    TGraphAsymmErrors* graph_Type0 = effPt_Type0->GetPaintedGraph();
    if (graph_Type0) {
        graph_Type0->GetYaxis()->SetRangeUser(0.0, 1.5);
    }

    effPt_Type3->SetMarkerStyle(22);
    effPt_Type3->SetMarkerColor(kBlue);
    effPt_Type3->SetLineColor(kBlue);
    //effPt_Type3->Draw("P SAME");

    TLegend* leg = new TLegend(0.75, 0.75, 0.90, 0.90);
    leg->AddEntry(effPt_Type0, "Type 0 (MFT+MCH+MID)", "ep");
    //leg->AddEntry(effPt_Type3, "Type 3 (MCH+MID)", "ep");
    leg->Draw();

    c1->SaveAs("tracking_efficiency_pt.png");
}

