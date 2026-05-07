// ============================================================
//  RebinTopDownTUnfold.C
// ============================================================

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------
// ヘルパー関数: Reco用 (Flat Stats になるようなビン境界を計算)
// ---------------------------------------------------------
std::vector<double> CalculateFlatStatsBins(TH1D* h, int nTargetBins, double xMin, double xMax)
{
    int binMin = h->FindBin(xMin);
    int binMax = h->FindBin(xMax);

    std::vector<double> cumSum;
    cumSum.push_back(0.0);
    for (int i = binMin; i <= binMax; ++i) {
        cumSum.push_back(cumSum.back() + h->GetBinContent(i));
    }

    double totalIntegral = cumSum.back();
    double step = totalIntegral / nTargetBins;

    std::vector<double> bins;
    bins.push_back(xMin);

    for (int k = 1; k < nTargetBins; ++k) {
        double target = k * step;
        int bestIdx = 1;
        double minDiff = std::abs(cumSum[1] - target);

        for (size_t i = 2; i < cumSum.size(); ++i) {
            double diff = std::abs(cumSum[i] - target);
            if (diff < minDiff) {
                minDiff = diff;
                bestIdx = i;
            }
        }

        double edge = h->GetBinLowEdge(binMin + bestIdx);
        if (edge > bins.back() && edge < xMax) {
            bins.push_back(edge);
        }
    }

    if (bins.back() < xMax) {
        bins.push_back(xMax);
    }

    return bins;
}

// ---------------------------------------------------------
// ヘルパー関数: Gen用 (TopDownでPurity/Stabilityを満たすビン境界を計算)
// ---------------------------------------------------------
std::vector<double> TopDownBinEdges(TH2D* hMatrixFine, double targetPurity, double targetStability, double minStats, double maxPt2)
{
    std::vector<double> edges;
    edges.push_back(maxPt2);

    int maxBinFine = hMatrixFine->GetXaxis()->FindBin(maxPt2 - 1e-6);
    int upperBin = maxBinFine;

    for (int lowerBin = maxBinFine; lowerBin >= 1; --lowerBin) {
        double diag = hMatrixFine->Integral(lowerBin, upperBin, lowerBin, upperBin);
        double sumTrue = hMatrixFine->Integral(lowerBin, upperBin, 0, hMatrixFine->GetNbinsY() + 1);
        double sumReco = hMatrixFine->Integral(0, hMatrixFine->GetNbinsX() + 1, lowerBin, upperBin);

        if (sumTrue <= 0 || sumReco <= 0) continue;

        double purity = diag / sumReco;
        double stability = diag / sumTrue;

        if (purity >= targetPurity && stability >= targetStability && sumReco >= minStats && sumTrue >= minStats) {
            double lowerEdge = hMatrixFine->GetXaxis()->GetBinLowEdge(lowerBin);
            edges.push_back(lowerEdge);
            upperBin = lowerBin - 1;
        }
    }

    if (upperBin >= 1) {
        if (edges.size() > 1) {
            std::cout << "[TopDown] Final segment [0.0, " << edges.back() << "] failed criteria. Merging with adjacent upper bin." << std::endl;
            edges.pop_back();
        }
    }

    if (edges.back() != 0.0) {
        edges.push_back(0.0);
    }

    std::reverse(edges.begin(), edges.end());
    return edges;
}

// ---------------------------------------------------------
// メイン関数
// ---------------------------------------------------------
void RebinTopDownTUnfold(
    const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Step1_merged.root",
    const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/TopDown/30/Step2_Rebinned.root",
    double targetPurity = 0.30,
    double targetStability = 0.30,
    double minStats = 500,
    double xMin = 0.0,
    double xMax = 0.020)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // --- 1. 元データの読み込み ---
    TFile* fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "[Error] Cannot open input file: " << inputFile << std::endl;
        return;
    }

    TH1D* hGenFine = (TH1D*)fIn->Get("hGenPt2");
    TH1D* hRecoFine = (TH1D*)fIn->Get("hRecoPt2");
    // TH2Fで保存されている場合も考慮し、汎用的なTH2として取得後、TH2Dにキャストまたは新規作成が安全ですが、
    // ここではご提示のコードに合わせて取得します（元の型がTH2Fの場合はキャストに注意してください）
    TH2D* hMatrixFine = (TH2D*)fIn->Get("hResponseMatrixPairPt2"); 

    if (!hGenFine || !hRecoFine || !hMatrixFine) {
        std::cerr << "[Error] Required histograms not found in input file!" << std::endl;
        fIn->Close();
        return;
    }

    // --- 2. Genビンの計算 (TopDown) ---
    std::vector<double> binsGen = TopDownBinEdges(hMatrixFine, targetPurity, targetStability, minStats, xMax);
    int nBinsGen = binsGen.size() - 1;

    if (nBinsGen <= 0) {
        std::cerr << "[Error] No Gen bins determined. Try loosening the criteria." << std::endl;
        fIn->Close();
        return;
    }

    std::cout << "========================================" << std::endl;
    std::cout << "[Info] Determined Gen Bins (TopDown): " << nBinsGen << std::endl;
    for (int i = 0; i <= nBinsGen; ++i) {
        std::cout << "  Gen Edge[" << i << "] = " << binsGen[i] << std::endl;
    }

    // --- 3. Recoビンの計算 (Genビン数の2倍をFlat Statsで) ---
    int nRecoTargetBins = nBinsGen * 2;
    std::vector<double> binsReco = CalculateFlatStatsBins(hRecoFine, nRecoTargetBins, xMin, xMax);
    int nBinsReco = binsReco.size() - 1;

    std::cout << "[Info] Determined Reco Bins (Flat Stats): " << nBinsReco << std::endl;
    std::cout << "========================================" << std::endl;

    // --- 4. 1DヒストグラムのRebin ---
    TH1D* hGenRebin = (TH1D*)hGenFine->Rebin(nBinsGen, "hGenPt2_rebin", binsGen.data());
    TH1D* hRecoRebin = (TH1D*)hRecoFine->Rebin(nBinsReco, "hRecoPt2_rebin", binsReco.data());

    // --- 5. 非正方(Rectangular)な応答行列(TH2D)の作成 ---
    TH2D* hMatrixRebin = new TH2D("hResponseMatrix", "Response Matrix (Gen:X, Reco:Y)",
                                  nBinsGen, binsGen.data(),
                                  nBinsReco, binsReco.data());

    for (int i = 1; i <= hMatrixFine->GetNbinsX(); ++i) {
        for (int j = 1; j <= hMatrixFine->GetNbinsY(); ++j) {
            double content = hMatrixFine->GetBinContent(i, j);
            if (content > 0) {
                double xCenter = hMatrixFine->GetXaxis()->GetBinCenter(i); // Gen
                double yCenter = hMatrixFine->GetYaxis()->GetBinCenter(j); // Reco
                hMatrixRebin->Fill(xCenter, yCenter, content);
            }
        }
    }

    // --- 6. 結果の保存 ---
    TFile* fOut = TFile::Open(outputFile, "RECREATE");
    hGenRebin->Write();
    hRecoRebin->Write();
    hMatrixRebin->Write();
    fOut->Close();
    fIn->Close();

    std::cout << "[Info] Saved ROOT file: " << outputFile << std::endl;
}