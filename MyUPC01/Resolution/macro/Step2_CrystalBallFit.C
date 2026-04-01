// ============================================================
//  Step2_CrystalBallFit.C
//  Step 2: Step1 で保存したリビン済み2Dヒストグラムを読み込み、
//          各 pT_true ビンでの 1D ヒストグラム（分解能分布）を取り出し、
//          Crystal Ball 関数でフィッティングして各ビンの結果を
//          個別にキャンバス表示・保存する。
//
//  ▼ Crystal Ball 関数 (片側 tail 版):
//   f(x) = N * exp(-(x-mu)^2/(2*sigma^2))        for (x-mu)/sigma > -alpha
//          N * (n/|alpha|)^n * exp(-|alpha|^2/2)
//              * (n/|alpha| - |alpha| - (x-mu)/sigma)^(-n)    otherwise
//   パラメータ: [0]=N(規格化), [1]=mu(平均), [2]=sigma, [3]=alpha, [4]=n
//
//  実行例:
//    root -l -b -q 'Step2_CrystalBallFit.C'
//
//  出力:
//    ../Step2_CrystalBallFit_bin{i}.png  … 各 pT ビンのフィット結果
//    ../Step2_fitresults.root            … フィット結果 (μ, σ) を TGraphErrors に保存
// ============================================================
#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMath.h"

// ------------------------------------------------------------------
// Crystal Ball 関数（ROOT TF1 用）
//   par[0] = N      (normalization)
//   par[1] = mu     (mean / peak position)
//   par[2] = sigma  (core width → resolution)
//   par[3] = alpha  (tail 開始点、> 0)
//   par[4] = n      (tail べき指数、> 1)
// ------------------------------------------------------------------
double CrystalBall(double* x, double* par)
{
    double N     = par[0];
    double mu    = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double n     = par[4];

    double t = (x[0] - mu) / sigma;
    if (t >= -TMath::Abs(alpha)) {
        return N * TMath::Exp(-0.5 * t * t);
    } else {
        double A = TMath::Power(n / TMath::Abs(alpha), n)
                 * TMath::Exp(-0.5 * alpha * alpha);
        double B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
        return N * A * TMath::Power(B - t, -n);
    }
}

void Step2_CrystalBallFit()
{
    // ----------------------------------------------------------
    // [設定]
    // ----------------------------------------------------------
    const TString inFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step1_merged.root";
    const TString histName = "h2Rebinned";

    const int minEntries = 30;   // フィットを行う最低エントリー数
    const double fitRange = 0.12; // 初期フィット幅 ±fitRange

    // ----------------------------------------------------------
    // [読み込み]
    // ----------------------------------------------------------
    TFile* fIn = TFile::Open(inFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "[Error] Cannot open: " << inFile << std::endl;
        return;
    }
    TH2D* h2 = (TH2D*)fIn->Get(histName);
    if (!h2) {
        std::cerr << "[Error] h2Rebinned not found in " << inFile << std::endl;
        return;
    }
    h2->SetDirectory(0);
    fIn->Close();

    int nBinsX = h2->GetNbinsX();

    // ----------------------------------------------------------
    // [スタイル]
    // ----------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    const int fontCode = 42;
    for (const char* ax : {"X", "Y", "Z"}) {
        gStyle->SetLabelFont(fontCode, ax);
        gStyle->SetTitleFont(fontCode, ax);
        gStyle->SetTitleSize(0.048, ax);
        gStyle->SetLabelSize(0.042, ax);
    }

    // ----------------------------------------------------------
    // [フィット結果格納グラフ]
    // ----------------------------------------------------------
    // gr_mu  : μ (平均 → バイアス) vs pT_true
    // gr_sig : σ (分解能)         vs pT_true
    TGraphErrors* gr_mu  = new TGraphErrors();
    TGraphErrors* gr_sig = new TGraphErrors();
    gr_mu ->SetName("gr_mu");
    gr_sig->SetName("gr_sig");
    int nPts = 0;

    // ----------------------------------------------------------
    // [各ビンのフィット]
    // ----------------------------------------------------------
    for (int i = 0; i < nBinsX; ++i) {
        TH1D* hProj = h2->ProjectionY(Form("hProj_step2_%d", i), i+1, i+1);

        double ptLo  = h2->GetXaxis()->GetBinLowEdge(i+1);
        double ptHi  = h2->GetXaxis()->GetBinUpEdge(i+1);
        double ptCen = h2->GetXaxis()->GetBinCenter(i+1);
        double ptErr = h2->GetXaxis()->GetBinWidth(i+1) / 2.0;

        if (hProj->GetEntries() < minEntries) {
            std::cout << Form("[Skip] pT bin %d (%.2f-%.2f GeV/c): entries=%.0f < %d",
                i, ptLo, ptHi, hProj->GetEntries(), minEntries) << std::endl;
            delete hProj;
            continue;
        }

        // --- 初期パラメータ推定 ---
        int    maxBin = hProj->GetMaximumBin();
        double peakX  = hProj->GetBinCenter(maxBin);
        double peakY  = hProj->GetMaximum();

        // 1st pass: Gauss で大まかな peak/sigma を推定
        TF1* fGaus = new TF1(Form("fGaus_%d", i), "gaus",
                              peakX - fitRange, peakX + fitRange);
        fGaus->SetParameters(peakY, peakX, 0.02);
        hProj->Fit(fGaus, "Q0R");
        double mu0    = fGaus->GetParameter(1);
        double sigma0 = std::abs(fGaus->GetParameter(2));
        delete fGaus;

        // 収束していなければ幅を広げて再試行
        if (sigma0 < 1e-4 || sigma0 > 0.5) { mu0 = peakX; sigma0 = 0.03; }

        // 2nd pass: Crystal Ball フィット
        // フィット範囲: ±3σ (推定値)
        double fitLo = mu0 - 3.0 * sigma0;
        double fitHi = mu0 + 3.0 * sigma0;

        TF1* fCB = new TF1(Form("fCB_%d", i), CrystalBall, fitLo, fitHi, 5);
        fCB->SetParNames("N", "#mu", "#sigma", "#alpha", "n");
        fCB->SetParameters(peakY, mu0, sigma0, 1.2, 5.0);
        // パラメータ範囲を設定して Physical な解を保証
        fCB->SetParLimits(0, 0, 1e9);        // N > 0
        fCB->SetParLimits(2, 1e-5, 1.0);     // sigma > 0
        fCB->SetParLimits(3, 0.1, 5.0);      // alpha > 0
        fCB->SetParLimits(4, 1.1, 50.0);     // n > 1

        int fitStatus = hProj->Fit(fCB, "Q0RS");  // R: 関数範囲でfit, S: 結果取得

        double mu_fit    = fCB->GetParameter(1);
        double sigma_fit = std::abs(fCB->GetParameter(2));
        double mu_err    = fCB->GetParError(1);
        double sigma_err = fCB->GetParError(2);
        double chi2ndf   = (fCB->GetNDF() > 0) ? fCB->GetChisquare()/fCB->GetNDF() : -1;

        // --- 描画 ---
        TCanvas* c2 = new TCanvas(Form("c2_bin%d", i),
                                  Form("Step2: Crystal Ball Fit, pT bin %d", i),
                                  700, 550);
        c2->SetLeftMargin(0.13);
        c2->SetBottomMargin(0.13);
        c2->SetTopMargin(0.10);
        c2->SetRightMargin(0.05);

        hProj->SetLineColor(kAzure+1);
        hProj->SetLineWidth(2);
        hProj->SetFillColorAlpha(kAzure+1, 0.25);
        hProj->SetFillStyle(1001);
        hProj->GetXaxis()->SetTitle("(p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}");
        hProj->GetYaxis()->SetTitle("Counts");
        hProj->GetXaxis()->SetRangeUser(-0.4, 0.4);
        hProj->Draw("HIST");

        fCB->SetLineColor(kRed+1);
        fCB->SetLineWidth(2);
        fCB->Draw("SAME");

        // 凡例
        TLegend* leg = new TLegend(0.55, 0.68, 0.92, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(fontCode);
        leg->SetTextSize(0.038);
        leg->AddEntry(hProj, "Data (projection)", "f");
        leg->AddEntry(fCB,   "Crystal Ball fit", "l");
        leg->Draw();

        // ラベル
        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(fontCode);
        tex.SetTextSize(0.042);
        tex.DrawLatex(0.15, 0.90, "#bf{ALICE Simulation}");
        tex.SetTextSize(0.036);
        tex.DrawLatex(0.15, 0.85, Form("p_{T}^{true} #in [%.2f, %.2f] GeV/c", ptLo, ptHi));
        tex.DrawLatex(0.15, 0.80, Form("#mu = %.4f #pm %.4f", mu_fit, mu_err));
        tex.DrawLatex(0.15, 0.75, Form("#sigma = %.4f #pm %.4f", sigma_fit, sigma_err));
        if (chi2ndf >= 0)
            tex.DrawLatex(0.15, 0.70, Form("#chi^{2}/ndf = %.2f", chi2ndf));
        tex.DrawLatex(0.15, 0.65, Form("Entries = %.0f", hProj->GetEntries()));

        TString outPng = Form("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_CrystalBallFit_bin%02d.png", i);
        c2->SaveAs(outPng);
        std::cout << "[Info] Saved: " << outPng
                  << Form("  mu=%.4f  sigma=%.4f  chi2/ndf=%.2f", mu_fit, sigma_fit, chi2ndf)
                  << std::endl;

        // フィット結果をグラフに追加
        gr_mu ->SetPoint(nPts, ptCen, mu_fit);
        gr_mu ->SetPointError(nPts, ptErr, mu_err);
        gr_sig->SetPoint(nPts, ptCen, sigma_fit);
        gr_sig->SetPointError(nPts, ptErr, sigma_err);
        ++nPts;

        delete leg;
        delete fCB;
        delete c2;
        delete hProj;
    }

    // ----------------------------------------------------------
    // [フィット結果を ROOT ファイルに保存]
    // ----------------------------------------------------------
    TFile* fOut = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_fitresults.root", "RECREATE");
    gr_mu ->Write("gr_mu");
    gr_sig->Write("gr_sigma");
    fOut->Close();
    std::cout << "[Info] Saved fit results: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_fitresults.root" << std::endl;

    delete gr_mu;
    delete gr_sig;
    delete h2;

    std::cout << "[Done] Step 2 finished. " << nPts << " bins fitted." << std::endl;
}
