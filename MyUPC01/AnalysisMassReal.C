#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>

void AnalysisMassReal() {
    // 統計ボックスを非表示
    gStyle->SetOptStat(0);
    
    // ROOTファイルを開く
    TFile* f = TFile::Open("/media/takuma/ESD-EAWA/Data/UPCcandMuon/AnalysisResults.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open AnalysisResults.root" << std::endl;
        return;
    }

    // ヒストグラムが保存されているディレクトリのパス（MyUPCTask.cxxの定義に対応）
    TString dir = "my-upc-01/";
    
    // 不変質量分布のヒストグラムを取得
    TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dir + "MMuonUnlike"));
    TH1* hLike = dynamic_cast<TH1*>(f->Get(dir + "MMuonLike"));

    if (!hUnlike || !hLike) {
        std::cerr << "Error: Histograms not found. Check the directory structure inside the root file." << std::endl;
        f->ls();
        return;
    }

    // ==========================================
    // 1. Raw Histograms (UnlikeとLikeを重ねて描画)
    // ==========================================
    TCanvas* c1 = new TCanvas("c1", "Raw Histograms", 800, 600);
    gPad->SetLogy(); // y軸をログスケールに
    gPad->SetGrid(); // グリッドを表示

    // Unlike Sign の描画設定
    hUnlike->SetTitle("Dimuon Invariant Mass; M_{#mu#mu} [GeV/c^{2}]; Counts");
    hUnlike->SetLineColor(kGreen+2);
    hUnlike->SetMarkerColor(kGreen+2);
    hUnlike->SetLineWidth(2);
    hUnlike->SetMarkerStyle(21); // 21: Full Square
    hUnlike->GetYaxis()->SetRangeUser(0.8, hUnlike->GetMaximum() * 5.0);
    hUnlike->Draw("E");

    // Like Sign の描画設定
    hLike->SetLineColor(kOrange+7);
    hLike->SetMarkerColor(kOrange+7);
    hLike->SetLineWidth(2);
    hLike->SetMarkerStyle(25); // 25: Open Square
    hLike->Draw("E SAME");

    // 凡例の設定
    TLegend* leg = new TLegend(0.65, 0.75, 0.90, 0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(hUnlike, "Unlike Sign (+-)", "lep");
    leg->AddEntry(hLike, "Like Sign (++ & --)", "lep");
    leg->Draw();

    c1->SaveAs("RawMass.png");

    // ==================================================
    // 2. Signal Extraction (Unlike - Like のバックグラウンド引き去り)
    // ==================================================
    TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
    hSignal->SetTitle("Signal (Unlike - Like); M_{#mu#mu} [GeV/c^{2}]; Counts");
    
    // Backgroundを引く (Like Signそのものをバックグラウンドとして減算)
    hSignal->Add(hLike, -1.0);

    TCanvas* c2 = new TCanvas("c2", "Signal Extraction", 800, 600);
    gPad->SetGrid();

    // 負の値があるかもしれないのでy軸の範囲を調整
    double ymin = hSignal->GetMinimum();
    if (ymin > 0) ymin = 0;
    
    hSignal->SetLineColor(kGreen+2);
    hSignal->SetMarkerColor(kGreen+2);
    hSignal->SetLineWidth(2);
    hSignal->SetMarkerStyle(21); // 21: Full Square
    hSignal->GetYaxis()->SetRangeUser(ymin * 1.2, hSignal->GetMaximum() * 1.5);
    hSignal->Draw("E");

    TLegend* leg2 = new TLegend(0.70, 0.80, 0.90, 0.89);
    leg2->SetBorderSize(0);
    leg2->AddEntry(hSignal, "Signal Extracted", "lep");
    leg2->Draw();

    c2->SaveAs("Unlike_bkg.png");

    // ファイルを閉じる
    f->Close();
}
