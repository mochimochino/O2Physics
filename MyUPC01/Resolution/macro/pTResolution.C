#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

// 共通処理を行うヘルパー関数
void processAndDraw(const std::vector<TString>& fileNames, const TString& histPath, const TString& title, const TString& outName, const TString& suffix) {
    
    double drawXmin = 0.0;
    double drawXmax = 5.0;
    double drawYmin = -0.5;
    double drawYmax =  0.5;

    double ptBins[] = {0.0, 0.5, 0.75, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 10.0};
    int nBinsX = sizeof(ptBins)/sizeof(double) - 1;

    TH2F *h2Pre = nullptr;

    for (const auto& fileName : fileNames) {
        TFile *f = TFile::Open(fileName, "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[Warning] File not found or broken: " << fileName << std::endl;
            if (f) f->Close();
            continue;
        }

        TH2F *h_temp = (TH2F*)f->Get(histPath);
        if (!h_temp) {
            std::cerr << "[Warning] 2D histogram not found in: " << fileName << " (" << histPath << ")" << std::endl;
            f->Close();
            continue;
        }

        if (!h2Pre) {
            h2Pre = (TH2F*)h_temp->Clone(Form("h2Pre_total_%s", suffix.Data()));
            h2Pre->SetDirectory(0);
        } else {
            h2Pre->Add(h_temp);
        }
        f->Close();
    }

    if (!h2Pre) {
        std::cerr << "[Error] No valid histograms could be loaded for " << title << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    int fontCode = 42;
    gStyle->SetLabelFont(fontCode, "XYZ");
    gStyle->SetTitleFont(fontCode, "XYZ");
    gStyle->SetTextFont(fontCode);
    gStyle->SetLegendFont(fontCode);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetLabelSize(0.03, "XYZ");
    gStyle->SetPalette(kBird);

    int nBinsY = h2Pre->GetNbinsY();
    double yMin = h2Pre->GetYaxis()->GetXmin();
    double yMax = h2Pre->GetYaxis()->GetXmax();

    TH2D *h2Rebinned = new TH2D(Form("h2Rebinned_%s", suffix.Data()), 
                                Form("%s; p_{T}^{true} (GeV/c);Resolution (p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}", title.Data()), 
                                nBinsX, ptBins, nBinsY, yMin, yMax);

    for (int ix = 1; ix <= h2Pre->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h2Pre->GetNbinsY(); ++iy) {
            double content = h2Pre->GetBinContent(ix, iy);
            if (content > 0) {
                double xCenter = h2Pre->GetXaxis()->GetBinCenter(ix);
                double yCenter = h2Pre->GetYaxis()->GetBinCenter(iy);
                h2Rebinned->Fill(xCenter, yCenter, content);
            }
        }
    }

    TGraphErrors *grFitResult = new TGraphErrors(nBinsX);

    for(int i = 0; i < nBinsX; ++i) {
        TH1D *hProj = h2Rebinned->ProjectionY(Form("hProj_%s_%d", suffix.Data(), i), i+1, i+1);

        if (hProj->GetEntries() < 10) continue;

        int maxBin = hProj->GetMaximumBin();
        double peakX = hProj->GetBinCenter(maxBin);
        double peakY = hProj->GetMaximum();
        
        TF1 *fitFunc = new TF1(Form("fit_%s_%d", suffix.Data(), i), "gaus", peakX - 0.1, peakX + 0.1);
        fitFunc->SetParameters(peakY, peakX, 0.02);
        hProj->Fit(fitFunc, "Q0R");

        double mean1  = fitFunc->GetParameter(1);
        double sigma1 = fitFunc->GetParameter(2);
        
        if(sigma1 > 0.0 && sigma1 < 1.0) {
            fitFunc->SetRange(mean1 - 1.5 * sigma1, mean1 + 1.5 * sigma1);
            hProj->Fit(fitFunc, "Q0R");
        }

        if(fitFunc) {
            double ptCenter = h2Rebinned->GetXaxis()->GetBinCenter(i+1);
            double ptErr    = h2Rebinned->GetXaxis()->GetBinWidth(i+1) / 2.0;
            double mean  = fitFunc->GetParameter(1);
            double sigma = fitFunc->GetParameter(2);

            grFitResult->SetPoint(i, ptCenter, mean);
            grFitResult->SetPointError(i, ptErr, sigma);
        }
        delete fitFunc;
    }

    TCanvas *c1 = new TCanvas(Form("c1_%s", suffix.Data()), Form("pT Resolution Fit Overlay - %s", suffix.Data()), 800, 600);
    c1->SetTopMargin(0.08);
    c1->SetRightMargin(0.05);
    c1->SetGrid();

    h2Rebinned->GetXaxis()->SetRangeUser(drawXmin, drawXmax);
    h2Rebinned->GetYaxis()->SetRangeUser(drawYmin, drawYmax);
    h2Rebinned->Draw("COL");

    grFitResult->SetMarkerStyle(20);
    grFitResult->SetMarkerSize(1.2);
    grFitResult->SetMarkerColor(kBlack);
    grFitResult->SetLineColor(kBlack);
    grFitResult->SetLineWidth(2);
    grFitResult->Draw("P SAME");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.65, 0.88, Form("Entries: %.0f", h2Pre->GetEntries()));
    
    c1->SaveAs(outName);

    delete c1;
    delete grFitResult;
    delete h2Rebinned;
    delete h2Pre;
}

void pTResolution() {
// Input files
    std::vector<TString> fileNames = {
        "jpsi-incoh.root",
        "jpsi-coh.root",
        "psi2s-incoh.root",
        "psi2s-coh.root",
        "psi2s-incoh-fd.root",
        "psi2s-coh-fd.root",
        "mumu-low.root",
        "mumu-mid.root",
        "mumu-high.root"
    }; 

    std::cout << "Processing Pre-Cut Data..." << std::endl;
    processAndDraw(
        fileNames, 
        "my-upc-muon-resolution/registry/hPtResoVsPtTrue", 
        "p_{T} Resolution Pre Cut", 
        "pTResolution_PreCut.png", 
        "pre"
    );

    std::cout << "Processing After-Cut Data..." << std::endl;
    processAndDraw(
        fileNames, 
        "my-upc-muon-resolution/registry/hPtResoVsPtTrue_PostCut", 
        "p_{T} Resolution After Cut", 
        "pTResolution_AfterCut.png", 
        "post"
    );

    std::cout << "All done!" << std::endl;
}