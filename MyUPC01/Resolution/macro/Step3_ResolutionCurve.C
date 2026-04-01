// ============================================================
//  Step3_ResolutionCurve.C
//  Step 3: Step2 のフィット結果を読み込み、
//          分解能カーブ（σ vs pT_true）と
//          バイアスカーブ（μ vs pT_true）を描画・保存する。
//
//  実行例:
//    root -l -b -q 'Step3_ResolutionCurve.C'
//
//  出力:
//    ../Step3_ResolutionCurve.png   … σ (%) vs pT_true カーブ
//    ../Step3_BiasCurve.png         … μ     vs pT_true カーブ
//    ../Step3_Combined.png          … σ と μ を上下 2 パッドに並べた複合図
// ============================================================
#include <iostream>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"

void Step3_ResolutionCurve()
{
    // ----------------------------------------------------------
    // [読み込み]
    // ----------------------------------------------------------
    const TString inFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step2_fitresults.root";
    TFile* fIn = TFile::Open(inFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "[Error] Cannot open: " << inFile << std::endl;
        return;
    }

    TGraphErrors* gr_sig = (TGraphErrors*)fIn->Get("gr_sigma");
    TGraphErrors* gr_mu  = (TGraphErrors*)fIn->Get("gr_mu");
    if (!gr_sig || !gr_mu) {
        std::cerr << "[Error] Graphs not found in " << inFile << std::endl;
        return;
    }
    gr_sig->SetDirectory(nullptr);
    gr_mu ->SetDirectory(nullptr);
    fIn->Close();

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
        gStyle->SetTitleSize(0.050, ax);
        gStyle->SetLabelSize(0.042, ax);
    }

    double xMin = 0.0, xMax = 5.5;

    // ----------------------------------------------------------
    // [グラフの体裁設定]
    // ----------------------------------------------------------
    auto styleGraph = [&](TGraphErrors* gr, int color, int mStyle) {
        gr->SetMarkerStyle(mStyle);
        gr->SetMarkerSize(1.3);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);
    };

    // σ: この値が「分解能」 (resolution)
    //    百分率で見たい場合は y *= 100 の変換が必要だが、
    //    ここでは元のΔpT/pT 比のまま表示する（横に %軸を補足ラベルで示す）
    styleGraph(gr_sig, kRed+1,    20);
    styleGraph(gr_mu,  kAzure+1,  21);

    // ----------------------------------------------------------
    // [図1] 分解能カーブ (σ vs pT_true)
    // ----------------------------------------------------------
    TCanvas* cSig = new TCanvas("cSig", "Resolution Curve", 800, 600);
    cSig->SetLeftMargin(0.13);
    cSig->SetBottomMargin(0.13);
    cSig->SetTopMargin(0.10);
    cSig->SetRightMargin(0.05);
    cSig->SetGrid();

    gr_sig->GetXaxis()->SetLimits(xMin, xMax);
    gr_sig->GetYaxis()->SetRangeUser(0.0, 0.15);
    gr_sig->GetXaxis()->SetTitle("p_{T}^{true} (GeV/c)");
    gr_sig->GetYaxis()->SetTitle("#sigma_{CB}  [= (#Deltap_{T}/p_{T})]");
    gr_sig->Draw("APE");

    // Zero 線
    TLine* l0sig = new TLine(xMin, 0, xMax, 0);
    l0sig->SetLineStyle(2);
    l0sig->SetLineColor(kGray+1);
    l0sig->Draw();

    TLegend* legSig = new TLegend(0.15, 0.78, 0.55, 0.88);
    legSig->SetBorderSize(0);
    legSig->SetFillStyle(0);
    legSig->SetTextFont(fontCode);
    legSig->SetTextSize(0.040);
    legSig->AddEntry(gr_sig, "CB #sigma (resolution core)", "pe");
    legSig->Draw();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(fontCode);
    tex.SetTextSize(0.045);
    tex.DrawLatex(0.14, 0.91, "#bf{ALICE Simulation}  #font[52]{(MC, UPC #mu pair)}");
    tex.SetTextSize(0.037);
    tex.DrawLatex(0.14, 0.86, "Crystal Ball fit to (p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}");

    cSig->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_ResolutionCurve.png");
    std::cout << "[Info] Saved: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_ResolutionCurve.png" << std::endl;

    // ----------------------------------------------------------
    // [図2] バイアスカーブ (μ vs pT_true)
    // ----------------------------------------------------------
    TCanvas* cMu = new TCanvas("cMu", "Bias Curve", 800, 600);
    cMu->SetLeftMargin(0.13);
    cMu->SetBottomMargin(0.13);
    cMu->SetTopMargin(0.10);
    cMu->SetRightMargin(0.05);
    cMu->SetGrid();

    gr_mu->GetXaxis()->SetLimits(xMin, xMax);
    gr_mu->GetYaxis()->SetRangeUser(-0.05, 0.05);
    gr_mu->GetXaxis()->SetTitle("p_{T}^{true} (GeV/c)");
    gr_mu->GetYaxis()->SetTitle("#mu_{CB}  [= (#Deltap_{T}/p_{T}) bias]");
    gr_mu->Draw("APE");

    TLine* l0mu = new TLine(xMin, 0, xMax, 0);
    l0mu->SetLineStyle(2);
    l0mu->SetLineColor(kGray+1);
    l0mu->Draw();

    TLegend* legMu = new TLegend(0.15, 0.78, 0.55, 0.88);
    legMu->SetBorderSize(0);
    legMu->SetFillStyle(0);
    legMu->SetTextFont(fontCode);
    legMu->SetTextSize(0.040);
    legMu->AddEntry(gr_mu, "CB #mu (momentum bias)", "pe");
    legMu->Draw();

    TLatex tex2;
    tex2.SetNDC();
    tex2.SetTextFont(fontCode);
    tex2.SetTextSize(0.045);
    tex2.DrawLatex(0.14, 0.91, "#bf{ALICE Simulation}  #font[52]{(MC, UPC #mu pair)}");
    tex2.SetTextSize(0.037);
    tex2.DrawLatex(0.14, 0.86, "Crystal Ball fit to (p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}");

    cMu->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_BiasCurve.png");
    std::cout << "[Info] Saved: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_BiasCurve.png" << std::endl;

    // ----------------------------------------------------------
    // [図3] σ と μ を上下2パッドに並べた複合図
    // ----------------------------------------------------------
    TCanvas* cComb = new TCanvas("cComb", "Resolution + Bias Combined", 800, 900);
    cComb->SetFillColor(0);

    // 上段パッド: σ(分解能)
    TPad* padTop = new TPad("padTop", "", 0.0, 0.43, 1.0, 1.0);
    padTop->SetLeftMargin(0.13);
    padTop->SetRightMargin(0.05);
    padTop->SetBottomMargin(0.02);
    padTop->SetTopMargin(0.14);
    padTop->SetGrid();
    padTop->Draw();
    padTop->cd();

    gr_sig->GetXaxis()->SetLimits(xMin, xMax);
    gr_sig->GetXaxis()->SetLabelSize(0);  // 上パッドはX軸ラベル非表示
    gr_sig->GetXaxis()->SetTitleSize(0);
    gr_sig->GetYaxis()->SetRangeUser(0.0, 0.15);
    gr_sig->GetYaxis()->SetTitle("#sigma_{CB}  (#Deltap_{T}/p_{T})");
    gr_sig->GetYaxis()->SetTitleSize(0.055);
    gr_sig->GetYaxis()->SetLabelSize(0.050);
    gr_sig->GetYaxis()->SetTitleOffset(1.1);
    gr_sig->Draw("APE");
    l0sig->Draw();

    TLegend* leg1 = new TLegend(0.15, 0.76, 0.60, 0.86);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->SetTextFont(fontCode); leg1->SetTextSize(0.053);
    leg1->AddEntry(gr_sig, "CB #sigma (resolution)","pe");
    leg1->Draw();

    TLatex texTop;
    texTop.SetNDC();
    texTop.SetTextFont(fontCode);
    texTop.SetTextSize(0.058);
    texTop.DrawLatex(0.14, 0.88, "#bf{ALICE Simulation}  #font[52]{(MC, UPC #mu pair)}");

    // 下段パッド: μ(バイアス)
    cComb->cd();
    TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.43);
    padBot->SetLeftMargin(0.13);
    padBot->SetRightMargin(0.05);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.20);
    padBot->SetGrid();
    padBot->Draw();
    padBot->cd();

    gr_mu->GetXaxis()->SetLimits(xMin, xMax);
    gr_mu->GetYaxis()->SetRangeUser(-0.05, 0.05);
    gr_mu->GetXaxis()->SetTitle("p_{T}^{true} (GeV/c)");
    gr_mu->GetYaxis()->SetTitle("#mu_{CB}  (bias)");
    gr_mu->GetXaxis()->SetTitleSize(0.075);
    gr_mu->GetXaxis()->SetLabelSize(0.065);
    gr_mu->GetYaxis()->SetTitleSize(0.070);
    gr_mu->GetYaxis()->SetLabelSize(0.065);
    gr_mu->GetYaxis()->SetTitleOffset(0.85);
    gr_mu->Draw("APE");
    l0mu->Draw();

    TLegend* leg2 = new TLegend(0.15, 0.82, 0.60, 0.96);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->SetTextFont(fontCode); leg2->SetTextSize(0.068);
    leg2->AddEntry(gr_mu, "CB #mu (bias)","pe");
    leg2->Draw();

    cComb->cd();
    cComb->SaveAs("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_Combined.png");
    std::cout << "[Info] Saved: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_Combined.png" << std::endl;

    // ----------------------------------------------------------
    // [グラフを ROOT ファイルにも保存]
    // ----------------------------------------------------------
    TFile* fOut = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_ResolutionCurve.root", "RECREATE");
    gr_sig->Write("gr_sigma");
    gr_mu ->Write("gr_mu");
    fOut->Close();
    std::cout << "[Info] Saved graphs: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/Step3_ResolutionCurve.root" << std::endl;

    delete leg1; delete leg2; delete legSig; delete legMu;
    delete l0sig; delete l0mu;
    delete cSig; delete cMu; delete cComb;
    delete gr_sig; delete gr_mu;

    std::cout << "[Done] Step 3 finished." << std::endl;
}
