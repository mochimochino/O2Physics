#!/usr/bin/env python3
import os
import ROOT

# ============================================================
# 1. 全体設定
# ============================================================
INPUT_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0403/without/"
OUTPUT_FILE = "mass_distribution_rebin_colored.png"

HISTO_PATH = "my-upc-mass-02/registry/hMassUnlike"

REBIN_FACTOR = 20

SAMPLES = [
    ("jpsi-incoh.root",     "J/#psi incoh.",         ROOT.kRed + 1),
    ("jpsi-coh.root",       "J/#psi coh.",           ROOT.kBlue + 1),
    ("psi2s-incoh.root",    "#psi(2S) incoh.",       ROOT.kGreen + 2),
    ("psi2s-coh.root",      "#psi(2S) coh.",         ROOT.kMagenta + 1),
    ("psi2s-incoh-fd.root", "#psi(2S) incoh. fd",    ROOT.kOrange + 1),
    ("psi2s-coh-fd.root",   "#psi(2S) coh. fd",      ROOT.kCyan + 2),
    ("mumu-low.root",       "#mu#mu low",            ROOT.kViolet + 1),
    ("mumu-mid.root",       "#mu#mu mid",            ROOT.kTeal + 2),
    ("mumu-high.root",      "#mu#mu high",           ROOT.kPink + 1)
]

def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    
    fontCode = 42
    ROOT.gStyle.SetLabelFont(fontCode, "X")
    ROOT.gStyle.SetLabelFont(fontCode, "Y")
    ROOT.gStyle.SetTitleFont(fontCode, "X")
    ROOT.gStyle.SetTitleFont(fontCode, "Y")
    ROOT.gStyle.SetTitleSize(0.048, "X")
    ROOT.gStyle.SetTitleSize(0.048, "Y")
    ROOT.gStyle.SetLabelSize(0.042, "X")
    ROOT.gStyle.SetLabelSize(0.042, "Y")

    canvas = ROOT.TCanvas("c1", "Mass Distribution", 1000, 650)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.06)
    canvas.SetRightMargin(0.05)
    canvas.SetGrid()
    canvas.SetLogy()

    legend = ROOT.TLegend(0.64, 0.55, 0.94, 0.94)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(fontCode)
    legend.SetTextSize(0.030)

    histograms = []
    first_draw = True
    max_y = 0.0

    print(f"--- 不変質量分布の描画を開始 (Rebin Factor: {REBIN_FACTOR}) ---")

    for filename, label, color in SAMPLES:
        filepath = os.path.join(INPUT_DIR, filename)
        if not os.path.exists(filepath):
            print(f"  [SKIP] ファイルが見つかりません: {filename}")
            continue

        root_file = ROOT.TFile.Open(filepath, "READ")
        if not root_file or root_file.IsZombie():
            print(f"  [ERROR] ファイルを開けません: {filename}")
            continue

        hist = root_file.Get(HISTO_PATH)
        if not hist:
            print(f"  [WARNING] ヒストグラムが見つかりません -> {HISTO_PATH} in {filename}")
            root_file.Close()
            continue

        h_clone = hist.Clone(f"h_{filename}")
        h_clone.SetDirectory(0) 
        
        # Rebin
        if REBIN_FACTOR > 1:
            h_clone.Rebin(REBIN_FACTOR)
            
        h_clone.SetLineColor(color)
        h_clone.SetLineWidth(2)
        
        entries = int(h_clone.GetEntries())
        label_with_entries = f"{label} ({entries})"
        
        histograms.append(h_clone)
        legend.AddEntry(h_clone, label_with_entries, "l")

        current_max = h_clone.GetMaximum()
        if current_max > max_y:
            max_y = current_max

        if first_draw:
            bin_width_gev = h_clone.GetBinWidth(1)
            bin_width_mev = bin_width_gev * 1000.0
            y_title = f"Counts / {bin_width_mev:.0f} MeV/#it{{c}}^{{2}}"
            
            h_clone.SetTitle(f";m_{{#mu#mu}} (GeV/#it{{c}}^{{2}});{y_title}")
            h_clone.GetXaxis().SetRangeUser(1.0, 10.0)
            h_clone.Draw("HIST")
            first_draw = False
        else:
            h_clone.Draw("HIST SAME")

        root_file.Close()
        print(f"  [OK] 読み込み完了: {filename} (Entries: {entries})")

    if not histograms:
        print("\n[ERROR] 描画するヒストグラムが一つもありません。")
        return

    histograms[0].SetMinimum(1e0) 
    histograms[0].SetMaximum(max_y * 100.0) 

    legend.Draw()
    
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(fontCode)
    tex.SetTextSize(0.036)
    tex.DrawLatex(0.16, 0.90, "#bf{Invariant Mass (Unlike)}")
    tex.DrawLatex(0.16, 0.85, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}")
    tex.DrawLatex(0.16, 0.80, "with Coherent Cut")
    
    out_path = os.path.join(INPUT_DIR, OUTPUT_FILE)
    canvas.SaveAs(out_path)
    print(f"\n--- 完了: {out_path} を生成しました ---")

if __name__ == "__main__":
    main()