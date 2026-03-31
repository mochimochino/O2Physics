import ROOT
import array
import sys

def fit_and_overlay_resolution():
    # ==========================================
    # Setting
    # ==========================================
    file_name = "AnalysisResults.root" 
    hist_path_pre = "my-upc-muon-resolution/registry/hPtResoVsPtTrue" 

    draw_xmin, draw_xmax = 0.0, 5.0
    draw_ymin, draw_ymax = -0.5, 0.5

    # PyROOTで可変長ビンを扱うための配列変換
    pt_bins_list = [0.0, 0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 5.0, 10.0]
    pt_bins = array.array('d', pt_bins_list)
    n_bins_x = len(pt_bins_list) - 1
    # ==========================================

    # --- 1. ファイルとヒストグラムの読み込み ---
    f = ROOT.TFile.Open(file_name, "READ")
    if not f or f.IsZombie():
        print(f"[Error] File not found: {file_name}", file=sys.stderr)
        return

    h2_pre = f.Get(hist_path_pre)
    if not h2_pre:
        print(f"[Error] 2D histogram not found: {hist_path_pre}", file=sys.stderr)
        return

    # --- 2. グラフの体裁（かっこよくする）設定 ---
    ROOT.gStyle.SetPalette(ROOT.kViridis) # モダンで美しい紫〜緑〜黄のカラーマップ
    ROOT.gStyle.SetNumberContours(256)    # グラデーションを滑らかに
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)            # デフォルトのタイトル枠を消す
    ROOT.gStyle.SetPadTickX(1)            # 上部の目盛りをオン
    ROOT.gStyle.SetPadTickY(1)            # 右側の目盛りをオン
    
    # フォントをHelvetica(42)に統一
    ROOT.gStyle.SetTextFont(42)
    for axis in ["X", "Y", "Z"]:
        ROOT.gStyle.SetLabelFont(42, axis)
        ROOT.gStyle.SetTitleFont(42, axis)
        ROOT.gStyle.SetTitleSize(0.045, axis)
        ROOT.gStyle.SetLabelSize(0.04, axis)
    ROOT.gStyle.SetTitleOffset(1.2, "Y")

    # --- 3. 2DヒストグラムのRebin ---
    n_bins_y = h2_pre.GetNbinsY()
    y_min = h2_pre.GetYaxis().GetXmin()
    y_max = h2_pre.GetYaxis().GetXmax()

    h2_rebinned = ROOT.TH2D("h2Rebinned", 
                            ";p_{T}^{true} (GeV/c);Resolution (p_{T}^{reco} - p_{T}^{true}) / p_{T}^{true}", 
                            n_bins_x, pt_bins, n_bins_y, y_min, y_max)

    for ix in range(1, h2_pre.GetNbinsX() + 1):
        for iy in range(1, h2_pre.GetNbinsY() + 1):
            content = h2_pre.GetBinContent(ix, iy)
            if content > 0:
                x_center = h2_pre.GetXaxis().GetBinCenter(ix)
                y_center = h2_pre.GetYaxis().GetBinCenter(iy)
                h2_rebinned.Fill(x_center, y_center, content)

    # --- 4. スライスと反復(Iterative)フィット ---
    gr_fit_result = ROOT.TGraphErrors(n_bins_x)

    for i in range(n_bins_x):
        h_proj = h2_rebinned.ProjectionY(f"hProj_{i}", i + 1, i + 1)

        # エントリー数が少なすぎる場合はスキップ
        if h_proj.GetEntries() < 10:
            continue

        max_bin = h_proj.GetMaximumBin()
        peak_x = h_proj.GetBinCenter(max_bin)
        peak_y = h_proj.GetMaximum()
        
        fit_func = ROOT.TF1(f"fit_{i}", "gaus", peak_x - 0.1, peak_x + 0.1)
        fit_func.SetParameters(peak_y, peak_x, 0.02)

        h_proj.Fit(fit_func, "Q0R")

        mean1  = fit_func.GetParameter(1)
        sigma1 = fit_func.GetParameter(2)
        
        if 0.0 < sigma1 < 1.0:
            fit_func.SetRange(mean1 - 1.5 * sigma1, mean1 + 1.5 * sigma1)
            h_proj.Fit(fit_func, "Q0R")

        if fit_func:
            pt_center = h2_rebinned.GetXaxis().GetBinCenter(i + 1)
            pt_err    = h2_rebinned.GetXaxis().GetBinWidth(i + 1) / 2.0
            
            mean  = fit_func.GetParameter(1)
            sigma = fit_func.GetParameter(2)

            gr_fit_result.SetPoint(i, pt_center, mean)
            gr_fit_result.SetPointError(i, pt_err, sigma)

    # --- 5. 描画と保存 ---
    c1 = ROOT.TCanvas("c1", "pT Resolution Fit Overlay", 800, 600)
    
    # プロット領域のマージンを綺麗に設定
    c1.SetRightMargin(0.05)
    c1.SetTopMargin(0.08)
    c1.SetLeftMargin(0.12)
    c1.SetBottomMargin(0.12)
    c1.SetGrid()

    h2_rebinned.GetXaxis().SetRangeUser(draw_xmin, draw_xmax)
    h2_rebinned.GetYaxis().SetRangeUser(draw_ymin, draw_ymax)

    # カラーバーなし(COL)で描画
    h2_rebinned.Draw("COL")
    # h2_rebinned.Draw("kThermometer")

    gr_fit_result.SetMarkerStyle(20)
    gr_fit_result.SetMarkerSize(1.2)
    gr_fit_result.SetMarkerColor(ROOT.kBlack)
    gr_fit_result.SetLineColor(ROOT.kBlack)
    gr_fit_result.SetLineWidth(2)
    gr_fit_result.Draw("P SAME")

    # 凡例 (透明＆枠なしにしてプロの仕上がりに)
    leg = ROOT.TLegend(0.65, 0.80, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) # 背景を透明にする
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry(gr_fit_result, "Fit: #mu #pm 1#sigma", "pe")
    leg.Draw()

    # TLatexを使ってタイトルやメタデータをカッコよく配置
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.045)
    tex.DrawLatex(0.15, 0.85, "#bf{ALICE Simulation}") # 太字
    tex.SetTextSize(0.035)
    tex.DrawLatex(0.15, 0.80, "p_{T} Resolution (Pre-Cut)")
    
    c1.SaveAs("Resolution_FitOverlay_PreCut_Rebinnedpy.png")

if __name__ == "__main__":
    fit_and_overlay_resolution()