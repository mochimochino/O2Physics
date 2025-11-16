#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TColor.h>

void setHistStyle(TH1* hist, int color, int markerStyle, TString title) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetLineWidth(2);
    hist->SetStats(0);
    hist->SetTitle(title);
    hist->GetYaxis()->SetTitle("counts");
}

void normalizeHist(TH1* hist) {
    if (!hist) return;
    double integral = hist->Integral("width");
    if (integral > 0) {
        hist->Scale(1.0 / integral, "width");
    }
}

void fake_muon(const char* inFileName = "AnalysisResults.root") {
    gStyle->SetOptStat(0); 

    TFile* file = TFile::Open(inFileName);
    if (!file || file->IsZombie()) {
        printf("Error: Cannot open file %s\n", inFileName);
        return;
    }

    TDirectory* dir = (TDirectory*)file->Get("fake-muon-comparator");
    if (!dir) {
        printf("Error: Cannot find directory 'histos' in file %s\n", inFileName);
        file->Close();
        return;
    }

    TH1F* hPt_All = (TH1F*)dir->Get("hPt_matchedQualityCuts");
    TH1F* hPt_True = (TH1F*)dir->Get("hPt_GlobalMuon_True");
    TH1F* hPt_Fake = (TH1F*)dir->Get("hPt_GlobalMuon_Fake");

    TH1F* hEta_All = (TH1F*)dir->Get("hEta_matchedQualityCuts");
    TH1F* hEta_True = (TH1F*)dir->Get("hEta_GlobalMuon_True");
    TH1F* hEta_Fake = (TH1F*)dir->Get("hEta_GlobalMuon_Fake");

    TH1F* hPhi_All = (TH1F*)dir->Get("hPhi_matchedQualityCuts");
    TH1F* hPhi_True = (TH1F*)dir->Get("hPhi_GlobalMuon_True");
    TH1F* hPhi_Fake = (TH1F*)dir->Get("hPhi_GlobalMuon_Fake");

    int colorTrue = TColor::GetColor("#76B7B2");
    int colorFake = TColor::GetColor("#F28E2B"); 
    
    setHistStyle(hPt_All, kBlack, kFullCircle, "p_{T} Comparison (matchedQualityCuts)");
    setHistStyle(hPt_True, colorTrue, kFullSquare, "p_{T} Comparison (matchedQualityCuts)");
    setHistStyle(hPt_Fake, colorFake, kFullTriangleUp, "p_{T} Comparison (matchedQualityCuts)");

    setHistStyle(hEta_All, kBlack, kFullCircle, "#eta Comparison (matchedQualityCuts)");
    setHistStyle(hEta_True, colorTrue, kFullSquare, "#eta Comparison (matchedQualityCuts)");
    setHistStyle(hEta_Fake, colorFake, kFullTriangleUp, "#eta Comparison (matchedQualityCuts)");

    setHistStyle(hPhi_All, kBlack, kFullCircle, "#phi Comparison (matchedQualityCuts)");
    setHistStyle(hPhi_True, colorTrue, kFullSquare, "#phi Comparison (matchedQualityCuts)");
    setHistStyle(hPhi_Fake, colorFake, kFullTriangleUp, "#phi Comparison (matchedQualityCuts)");

    TCanvas* cPt = new TCanvas("cPt", "pT Comparison", 800, 600);
    cPt->SetLogy();
    hPt_All->Draw("HIST P");
    hPt_True->Draw("HIST P SAME");
    hPt_Fake->Draw("HIST P SAME");
    hPt_All->GetYaxis()->SetRangeUser(hPt_Fake->GetMinimum(1e-6) * 0.1, hPt_All->GetMaximum() * 10);

    TLegend* legPt = new TLegend(0.70, 0.70, 0.90, 0.90);
    legPt->SetHeader("Global Muons (Type 0)");
    legPt->AddEntry(hPt_All, "All (matchedQualityCuts)", "p");
    legPt->AddEntry(hPt_True, "True Match", "p");
    legPt->AddEntry(hPt_Fake, "Fake Match", "p");
    legPt->Draw();
    //cPt->SaveAs("compare_pT.png");
    cPt->SaveAs("compare_pT.pdf");
    delete cPt;

    TCanvas* cEta = new TCanvas("cEta", "Eta Comparison", 800, 600);
    hEta_All->Draw("HIST P");
    hEta_True->Draw("HIST P SAME");
    hEta_Fake->Draw("HIST P SAME");
    hEta_All->GetYaxis()->SetRangeUser(0, hEta_All->GetMaximum() * 1.5);

    TLegend* legEta = new TLegend(0.70, 0.70, 0.90, 0.90);
    legEta->SetHeader("Global Muons (Type 0)");
    legEta->AddEntry(hEta_All, "All (matchedQualityCuts)", "p");
    legEta->AddEntry(hEta_True, "True Match", "p");
    legEta->AddEntry(hEta_Fake, "Fake Match", "p");
    legEta->Draw();
    //cEta->SaveAs("compare_Eta.png");
    cEta->SaveAs("compare_Eta.pdf");
    delete cEta;

    TCanvas* cPhi = new TCanvas("cPhi", "Phi Comparison", 800, 600);
    hPhi_All->Draw("HIST P");
    hPhi_True->Draw("HIST P SAME");
    hPhi_Fake->Draw("HIST P SAME");
    hPhi_All->GetYaxis()->SetRangeUser(0, hPhi_All->GetMaximum() * 1.5);

    TLegend* legPhi = new TLegend(0.70, 0.70, 0.90, 0.90);
    legPhi->SetHeader("Global Muons (Type 0)");
    legPhi->AddEntry(hPhi_All, "All (matchedQualityCuts)", "p");
    legPhi->AddEntry(hPhi_True, "True Match", "p");
    legPhi->AddEntry(hPhi_Fake, "Fake Match", "p");
    legPhi->Draw();
    //cPhi->SaveAs("compare_Phi.png");
    cPhi->SaveAs("compare_Phi.pdf");
    delete cPhi;


    file->Close();
}