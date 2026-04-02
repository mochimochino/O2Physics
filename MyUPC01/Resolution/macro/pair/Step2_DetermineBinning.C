// ============================================================
//  Step2_DetermineBinning.C
// ============================================================
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include <algorithm>

// ==========================================================
// Top-down binning
// ==========================================================
std::vector<double> AutoDetermineBinningTopDown(TH2F* hMatrixFine, double targetPurity, double targetStability, double maxPt2) {
    std::vector<double> edges;
    edges.push_back(maxPt2);

    int maxBinFine = hMatrixFine->GetXaxis()->FindBin(maxPt2 - 1e-6);
    int upperBin = maxBinFine;

    for (int lowerBin = maxBinFine; lowerBin >= 1; --lowerBin) {
        
        // ★ 進捗の表示 (100ビンごとに更新して負荷を下げる)
        if (lowerBin % 100 == 0 || lowerBin == 1) {
            int progress = 100 - (int)(100.0 * lowerBin / maxBinFine);
            std::cout << "\r[1/4] Auto Binning... " << progress << "% completed" << std::flush;
        }

        double diag = hMatrixFine->Integral(lowerBin, upperBin, lowerBin, upperBin);
        double sumTrue = hMatrixFine->Integral(lowerBin, upperBin, 1, hMatrixFine->GetNbinsY());
        double sumReco = hMatrixFine->Integral(1, hMatrixFine->GetNbinsX(), lowerBin, upperBin);

        if (sumTrue <= 0 || sumReco <= 0) continue;

        double purity = diag / sumReco;
        double stability = diag / sumTrue;

        if (purity >= targetPurity && stability >= targetStability) {
            double lowerEdge = hMatrixFine->GetXaxis()->GetBinLowEdge(lowerBin);
            if (lowerBin > 1 && lowerEdge > 0.0) {
                edges.push_back(lowerEdge);
                upperBin = lowerBin - 1; 
            }
        }
    }
    std::cout << "\r[1/4] Auto Binning... 100% completed!       " << std::endl; // ★ 完了

    if (edges.back() != 0.0) edges.push_back(0.0);
    std::reverse(edges.begin(), edges.end());
    return edges;
}

void Step2_DetermineBinning()
{
    // ----------------------------------------------------------
    // Settings
    // ----------------------------------------------------------
    const TString inDir   = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0402/withoutCut000001/";
    const TString inFile  = inDir + "Step1_merged.root";
    const TString outDir  = inDir;

    gStyle->SetOptStat(0);

    TFile* fIn = TFile::Open(inFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Error: Cannot open " << inFile << std::endl;
        return;
    }
    TH2F* hMatrixFine = (TH2F*)fIn->Get("hResponseMatrixPairPt2"); // X:MC, Y:Reco
    if (!hMatrixFine) {
        std::cerr << "Error: hResponseMatrixPairPt2 not found." << std::endl;
        return;
    }

    double targetPurity = 0.50;    // 50%
    double targetStability = 0.50; // 50%
    double maxPt2 = 2.50;
    
    std::vector<double> pt2Bins = AutoDetermineBinningTopDown(hMatrixFine, targetPurity, targetStability, maxPt2);
    const int nBins = pt2Bins.size() - 1;

    std::cout << "[Auto Binning Top-Down] Generated " << nBins << " bins: {";
    for (size_t i = 0; i < pt2Bins.size(); ++i) {
        std::cout << pt2Bins[i] << (i == pt2Bins.size() - 1 ? "" : ", ");
    }
    std::cout << "}" << std::endl;
    // ==========================================================
    // Resolution evaluation
    // ==========================================================
    TGraphErrors* grResolution = new TGraphErrors();
    grResolution->SetName("grResolution");
    grResolution->SetTitle("Absolute p_{T}^{2} Resolution vs True p_{T}^{2};True p_{T}^{2} (GeV^{2}/c^{2});#sigma(p_{T}^{2}) (GeV^{2}/c^{2})");
    grResolution->SetMarkerStyle(20);
    grResolution->SetMarkerColor(kBlue+1);
    grResolution->SetLineColor(kBlue+1);

    int nSlices = 20; // slice number
    //double maxPt2 = 2.50;
    for (int i = 0; i < nSlices; ++i) {
        std::cout << "\r[2/4] Resolution Fitting... " << (int)(100.0 * i / nSlices) << "%" << std::flush;
        double pt2Min = i * (maxPt2 / nSlices);
        double pt2Max = (i + 1) * (maxPt2 / nSlices);
        double pt2Center = (pt2Min + pt2Max) / 2.0;

        int binMin = hMatrixFine->GetXaxis()->FindBin(pt2Min + 1e-6);
        int binMax = hMatrixFine->GetXaxis()->FindBin(pt2Max - 1e-6);
        
        // True pT^2 slice -> Reco pT^2 distribution
        TH1D* hSlice = hMatrixFine->ProjectionY(Form("slice_%d", i), binMin, binMax);
        if (hSlice->GetEntries() > 50) {
            // Gaussian fit --------------------------------(Crystall Ball is better?)
            hSlice->Fit("gaus", "Q0");
            TF1* fitFunc = hSlice->GetFunction("gaus");
            if (fitFunc) {
                double sigma = fitFunc->GetParameter(2);
                double sigmaErr = fitFunc->GetParError(2);
                int np = grResolution->GetN();
                grResolution->SetPoint(np, pt2Center, sigma);
                grResolution->SetPointError(np, 0.0, sigmaErr);
            }
        }
        delete hSlice;
    }

    TCanvas* cRes = new TCanvas("cRes", "Resolution", 800, 600);
    cRes->SetGrid();
    grResolution->Draw("APE");
    cRes->SaveAs(outDir + "Step2_Resolution_Abs.png");

    // ==========================================================
    // Rebinning
    // ==========================================================
    TH2D* hMatrixRebinned = new TH2D("hMatrixRebinned", 
                                     "Rebinned Response Matrix;True p_{T}^{2};Reco p_{T}^{2}", 
                                     nBins, pt2Bins.data(), nBins, pt2Bins.data());

    for (int ix = 1; ix <= hMatrixFine->GetNbinsX(); ++ix) {
        if (ix % 100 == 0 || ix == hMatrixFine->GetNbinsX()) {
            std::cout << "\r[3/4] Rebinning Matrix... " << (int)(100.0 * ix / hMatrixFine->GetNbinsX()) << "%" << std::flush;
        }
        double truePt2 = hMatrixFine->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hMatrixFine->GetNbinsY(); ++iy) {
            double recoPt2 = hMatrixFine->GetYaxis()->GetBinCenter(iy);
            double content = hMatrixFine->GetBinContent(ix, iy);
            if (content > 0) {
                hMatrixRebinned->Fill(truePt2, recoPt2, content);
            }
        }
    }
    std::cout << "\r[3/4] Rebinning Matrix... 100% completed!" << std::endl; // ★ 完了

// ==========================================================
    // 描画用のパーセント行列 (Percentage Matrix) の作成
    // ==========================================================
    TH2D* hMatrixPercentage = (TH2D*)hMatrixRebinned->Clone("hMatrixPercentage");
    hMatrixPercentage->SetTitle("Rebinned Response Matrix (Reco + Gen)/Gen in %);True p_{T}^{2};Reco p_{T}^{2}");
    
    // TEXTで表示する数字のフォーマットを指定 (例: 小数点第1位まで "95.2")
    gStyle->SetPaintTextFormat(".1f");

    for (int ix = 1; ix <= hMatrixRebinned->GetNbinsX(); ++ix) {
        std::cout << "\r[4/4] Purity & Stability... " << (int)(100.0 * ix / hMatrixRebinned->GetNbinsX()) << "%" << std::flush;
        // Gen (Trueビンの合計) を計算
        double sumGen = hMatrixRebinned->Integral(ix, ix, 1, hMatrixRebinned->GetNbinsY());
        
        for (int iy = 1; iy <= hMatrixRebinned->GetNbinsY(); ++iy) {
            double content = hMatrixRebinned->GetBinContent(ix, iy);
            if (sumGen > 0) {
                // (Rec+Gen)/Gen * 100 をセット
                hMatrixPercentage->SetBinContent(ix, iy, 100.0 * content / sumGen);
            } else {
                hMatrixPercentage->SetBinContent(ix, iy, 0.0);
            }
        }
    }
    std::cout << "\r[4/4] Purity & Stability... 100% completed!" << std::endl; // ★ 完了

    TCanvas* cMat = new TCanvas("cMat", "Rebinned Matrix", 800, 600);
    cMat->SetLogz();
    // 描画には hMatrixPercentage を使う
    hMatrixPercentage->Draw("COLZ TEXT");
    cMat->SaveAs(outDir + "Step2_RebinnedMatrix.png");

    // ==========================================================
    // 0.0 - 0.01 Zoom
    // ==========================================================
    TCanvas* cMatZoom = new TCanvas("cMatZoom", "Rebinned Matrix Zoom", 800, 600);
    cMatZoom->SetLogz();
    hMatrixPercentage->GetXaxis()->SetRangeUser(0.0, 0.01);
    hMatrixPercentage->GetYaxis()->SetRangeUser(0.0, 0.01);
    hMatrixPercentage->SetMarkerSize(1.5); 
    hMatrixPercentage->Draw("COLZ TEXT");
    cMatZoom->SaveAs(outDir + "Step2_RebinnedMatrix_Zoom.png");

    // Reset zoom
    hMatrixPercentage->GetXaxis()->UnZoom();
    hMatrixPercentage->GetYaxis()->UnZoom();
    hMatrixPercentage->SetMarkerSize(1.0);

    // ==========================================================
    // Purity and Stability
    // ==========================================================
    TH1D* hPurity = new TH1D("hPurity", "Purity and Stability;p_{T}^{2} Bin;Percentage", nBins, pt2Bins.data());
    TH1D* hStability = new TH1D("hStability", "Stability", nBins, pt2Bins.data());

    hPurity->SetLineColor(kRed);     hPurity->SetLineWidth(2);
    hStability->SetLineColor(kBlue); hStability->SetLineWidth(2);

    for (int i = 1; i <= nBins; ++i) {
        double diag = hMatrixRebinned->GetBinContent(i, i);
        
        // Purity = Diag / sumReco
        double sumReco = hMatrixRebinned->Integral(1, nBins, i, i); 
        double purity = (sumReco > 0) ? (diag / sumReco) : 0;
        hPurity->SetBinContent(i, purity);

        // Stability = Diag / sumTrue
        double sumTrue = hMatrixRebinned->Integral(i, i, 1, nBins);
        double stability = (sumTrue > 0) ? (diag / sumTrue) : 0;
        hStability->SetBinContent(i, stability);
    }

    TCanvas* cPurStab = new TCanvas("cPurStab", "Purity and Stability", 800, 600);
    cPurStab->SetGrid();
    hPurity->GetYaxis()->SetRangeUser(0.0, 1.1);
    hPurity->Draw("HIST");
    hStability->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.2, 0.85, 0.35);
    leg->AddEntry(hPurity, "Purity", "l");
    leg->AddEntry(hStability, "Stability", "l");
    leg->Draw();

    TF1* line50 = new TF1("line50", "0.5", 0, maxPt2);
    line50->SetLineStyle(2); line50->SetLineColor(kBlack);
    line50->Draw("SAME");

    cPurStab->SaveAs(outDir + "Step2_PurityStability.png");

    // ==========================================================
    // Purity & Stability Zoom
    // ==========================================================
    TCanvas* cPurStabZoom = new TCanvas("cPurStabZoom", "Purity and Stability Zoom", 800, 600);
    cPurStabZoom->SetGrid();
    hPurity->GetXaxis()->SetRangeUser(0.0, 0.01);
    hPurity->Draw("HIST");
    hStability->Draw("HIST SAME");
    leg->Draw();
    line50->Draw("SAME");
    cPurStabZoom->SaveAs(outDir + "Step2_PurityStability_Zoom.png");

    hPurity->GetXaxis()->UnZoom();

    // ==========================================================
    // Save output for Unfolding macro
    // ==========================================================
    TFile* fOut = TFile::Open(outDir + "Response_for_Unfolding.root", "RECREATE");
    hMatrixRebinned->Write("hResponseMatrix");
    hPurity->Write();
    hStability->Write();
    fOut->Close();

    std::cout << "[Done] Step 2 finished. Binning optimized and Response Matrix ready." << std::endl;
}