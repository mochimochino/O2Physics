void efficiency()
{
    // 1. 解析結果ファイルを開く
    TFile* f = TFile::Open("AnalysisResults.root");
    if (!f || f->IsZombie()) {
        printf("Error: Cannot open AnalysisResults.root\n");
        return;
    }

    // 2. ヒストグラムを取得
    // "my_muon_tracking_efficiency" は HistogramRegistry で指定した名前 (TBrowserで確認したもの)
    
    // TH1F* hNumPt_All = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Reco_Matched_All"); // pT_ALL は描画しない
    TH1F* hDenPt = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Gen");

    TH1F* hNumPt_Type0 = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Reco_Matched_Type0");
    TH1F* hNumPt_Type2 = (TH1F*)f->Get("my_muon_tracking_efficiency/hPt_Reco_Matched_Type2");

    // 3. ヒストグラムの存在をチェック
    if (/*!hNumPt_All ||*/ !hDenPt || !hNumPt_Type0 || !hNumPt_Type2) {
        printf("Error: Cannot find histograms in file.\n");
        if (!hDenPt) printf("... hPt_Gen not found.\n");
        if (!hNumPt_Type0) printf("... hPt_Reco_Matched_Type0 not found.\n");
        if (!hNumPt_Type2) printf("... hPt_Reco_Matched_Type2 not found.\n");
        return;
    }

    // --- TH1F::Divide を "B" オプション付きで使用 ---

    // 1. All Types (コメントアウト)
    // TH1F* hEffPt_All = (TH1F*)hNumPt_All->Clone("hEffPt_All");
    // hEffPt_All->Divide(hEffPt_All, hDenPt, 1.0, 1.0, "B"); 
    // hEffPt_All->SetTitle("Tracking Efficiency (All Types); p_{T}^{gen} [GeV/c]; Efficiency #epsilon");
    // hEffPt_All->SetStats(kFALSE); 

    // 2. Type 0
    TH1F* hEffPt_Type0 = (TH1F*)hNumPt_Type0->Clone("hEffPt_Type0");
    // 正しい二項分布誤差で割り算
    hEffPt_Type0->Divide(hEffPt_Type0, hDenPt, 1.0, 1.0, "B");
    hEffPt_Type0->SetTitle("Tracking Efficiency (Type 0); p_{T}^{gen} [GeV/c]; Efficiency #epsilon");
    hEffPt_Type0->SetStats(kFALSE);

    // 3. Type 2
    TH1F* hEffPt_Type2 = (TH1F*)hNumPt_Type2->Clone("hEffPt_Type2");
    // 正しい二項分布誤差で割り算
    hEffPt_Type2->Divide(hEffPt_Type2, hDenPt, 1.0, 1.0, "B");
    hEffPt_Type2->SetTitle("Tracking Efficiency (Type 2); p_{T}^{gen} [GeV/c]; Efficiency #epsilon");
    hEffPt_Type2->SetStats(kFALSE);


    // --- 描画 ---
    TCanvas* c1 = new TCanvas("c1", "Tracking Efficiency vs pT", 800, 600);
    c1->cd();
    gPad->SetGrid(); 

    // 最初に描画するヒストグラム (Type0) で軸を設定
    hEffPt_Type0->SetMarkerStyle(21);
    hEffPt_Type0->SetMarkerColor(kRed);
    hEffPt_Type0->SetLineColor(kRed);
    hEffPt_Type0->Draw("E1 P"); // エラーバー付きで描画

    // Y軸の範囲を 0.0 から 1.2 (120%) に設定
    hEffPt_Type0->GetYaxis()->SetRangeUser(0.0, 1.2); 

    // Type 2 を重ね描き
    hEffPt_Type2->SetMarkerStyle(22);
    hEffPt_Type2->SetMarkerColor(kBlue);
    hEffPt_Type2->SetLineColor(kBlue);
    hEffPt_Type2->Draw("E1 P SAME"); // "SAME" で重ね描き

    // 凡例 (凡例) の追加
    TLegend* leg = new TLegend(0.6, 0.2, 0.88, 0.4);
    // leg->AddEntry(hEffPt_All, "All Types", "ep"); // コメントアウト
    leg->AddEntry(hEffPt_Type0, "Type 0 (MFT+MCH+MID)", "ep");
    leg->AddEntry(hEffPt_Type2, "Type 2 (MFT+MCH)", "ep");
    leg->Draw();

    c1->SaveAs("tracking_efficiency_pt.png");
}

