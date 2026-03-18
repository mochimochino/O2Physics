#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"

// 1. Expo3
Double_t Draw_Expo3(Double_t* x, Double_t* par) {
    double m = x[0];
    return par[0] * TMath::Exp(par[1] * m + par[2] * m * m);
}

// 2. ExpoPol4
Double_t Draw_ExpoPol4(Double_t* x, Double_t* par) {
    double m = x[0];
    double poly = par[2] + par[3] * m + par[4] * m * m + par[5] * m * m * m + par[6] * m * m * m * m;
    return par[0] * TMath::Exp(par[1] * m) * poly;
}

// 3. VMG (Variable Width Gaussian)
Double_t Draw_VMG(Double_t* x, Double_t* par) {
    double m = x[0];
    double mbar = par[1];
    double sigma = par[2] + par[3] * (m - mbar) / mbar;
    if (sigma <= 0) return 0;
    return par[0] * TMath::Exp(-(m - mbar) * (m - mbar) / (2.0 * sigma * sigma));
}

void DrawBkgFunctions() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.04);

    double minMass = 1.8;
    double maxMass = 7.0;

    TCanvas* c1 = new TCanvas("c1", "Background Functions", 800, 600);
    // c1->SetLogy();
    c1->SetGrid();

    // Expo3の設定
    TF1* fExpo3 = new TF1("fExpo3", Draw_Expo3, minMass, maxMass, 3);
    fExpo3->SetParameters(5000, -1.2, 0.05);
    fExpo3->SetLineColor(kBlue);
    fExpo3->SetLineWidth(3);
    fExpo3->SetTitle("Background Function Shapes; M_{#mu#mu} [GeV/c^{2}]; Arbitrary Counts");

    // ExpoPol4の設定
    TF1* fExpoPol4 = new TF1("fExpoPol4", Draw_ExpoPol4, minMass, maxMass, 7);
    fExpoPol4->SetParameters(4000, -1.5, 1.0, 0.5, 0.1, -0.01, 0.001);
    fExpoPol4->SetLineColor(kRed);
    fExpoPol4->SetLineStyle(2);
    fExpoPol4->SetLineWidth(3);

    // VMGの設定
    TF1* fVMG = new TF1("fVMG", Draw_VMG, minMass, maxMass, 4);
    fVMG->SetParameters(10000, 1.0, 1.5, 0.3);
    fVMG->SetLineColor(kGreen+2);
    fVMG->SetLineStyle(7);
    fVMG->SetLineWidth(3);

    // 描画
    fExpo3->GetYaxis()->SetRangeUser(10, 20000);
    fExpo3->Draw();
    fExpoPol4->Draw("SAME");
    fVMG->Draw("SAME");

    // 凡例
    TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->AddEntry(fExpo3, "Expo3 (Nominal)", "l");
    leg->AddEntry(fExpoPol4, "ExpoPol4 (Systematics 1)", "l");
    leg->AddEntry(fVMG, "VMG (Systematics 2)", "l");
    leg->Draw();

    c1->SaveAs("Bkg_Functions_Comparison.png");
}