#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>

void AnalysisMassReal() {
    // グラフのグローバルスタイル設定
    gStyle->SetOptStat(0); // 統計ボックスは非表示
    gStyle->SetTitleFontSize(0.04); // タイトルの文字サイズを小さくして被りを防ぐ
    
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
    // キャンバスのマージンを調整し、グラフが占める割合を増やす
    TCanvas* c1 = new TCanvas("c1", "Raw Histograms", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.12);
    
    gPad->SetLogy(); // y軸をログスケールに
    gPad->SetGrid(); // グリッドを表示

    // 落ち着いた色合いの定義（HEXコード）
    int colUnlike = TColor::GetColor("#2e8b57"); // SeaGreen（落ち着いた緑）
    int colLike   = TColor::GetColor("#d2691e"); // Chocolate（落ち着いたオレンジ）

    // Unlike Sign の描画設定
    hUnlike->SetTitle("Dimuon Invariant Mass; M_{#mu#mu} [GeV/c^{2}]; Counts");
    hUnlike->SetLineColor(colUnlike);
    hUnlike->SetMarkerColor(colUnlike);
    hUnlike->SetLineWidth(2);
    hUnlike->SetMarkerStyle(21); // 21: Full Square
    hUnlike->GetYaxis()->SetRangeUser(0.8, hUnlike->GetMaximum() * 5.0);
    hUnlike->GetYaxis()->SetTitleOffset(1.2);
    hUnlike->Draw("E"); // 最初に描画したものの統計ボックスが表示される

    // Like Sign の描画設定
    hLike->SetLineColor(colLike);
    hLike->SetMarkerColor(colLike);
    hLike->SetLineWidth(2);
    hLike->SetMarkerStyle(25); // 25: Open Square
    hLike->Draw("E SAME");

    // 凡例の設定（右上にピッタリ配置）
    TLegend* leg = new TLegend(0.65, 0.80, 0.95, 0.95);
    leg->SetBorderSize(1);
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

    // キャンバス2の設定
    TCanvas* c2 = new TCanvas("c2", "Signal Extraction", 800, 600);
    c2->SetLeftMargin(0.12);
    c2->SetRightMargin(0.05);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.12);
    
    gPad->SetGrid();

    // 負の値があるかもしれないのでy軸の範囲を調整
    double ymin = hSignal->GetMinimum();
    if (ymin > 0) ymin = 0;
    
    int colSignal = TColor::GetColor("#2e8b57"); // シグナルも落ち着いた緑に
    
    hSignal->SetLineColor(colSignal);
    hSignal->SetMarkerColor(colSignal);
    hSignal->SetLineWidth(2);
    hSignal->SetMarkerStyle(21); // 21: Full Square
    hSignal->GetYaxis()->SetRangeUser(ymin * 1.2, hSignal->GetMaximum() * 1.5);
    hSignal->GetYaxis()->SetTitleOffset(1.2);
    hSignal->Draw("E");

    // 凡例の設定（右上にピッタリ配置）
    TLegend* leg2 = new TLegend(0.65, 0.85, 0.95, 0.95);
    leg2->SetBorderSize(1);
    leg2->AddEntry(hSignal, "Signal Extracted", "lep");
    leg2->Draw();

    c2->SaveAs("Unlike_bkg.png");

    // ファイルを閉じる
    f->Close();
}
