import array
import glob
import math
import os

import ROOT

# ============================================================
#  RebinFlatStatsUnbinned.py
#  UPCMuonAnalysisMC.cxx
# "upc-muon-analysis": {"processExportOmniFoldTuple": "true"})
# ============================================================


def RebinFlatStatsUnbinned(
    inputFiles="AnalysisResults.root",
    outputFile="Rebinned_Unbinned.root",
    treeName="omniFoldPairTuple",
    nTargetBins=5,
    xMin=0.0,
    xMax=1.40,
    # 等統計(Flat-Stats)分割を実際に行う範囲。[xMin, xMax] 全体ではなく
    # [splitMin, splitMax] の内側だけを nTargetBins 個に等統計分割する。
    # Coherent が支配的な低 p_{T}^{2} 領域(今回は splitMin 未満)は分割対象
    # から外し、[xMin, splitMin) をひとまとめの 1 bin として残す(結果の
    # 2D Migration 確率プロットの最小値は xMin=0 のまま)。
    # splitMax=None なら xMax と同じ(上側には lump bin を作らない)。
    splitMin=0.3,
    splitMax=None,
    nBinsPtQA=40,
    # 重ね書き(Overlay QA)専用の表示レンジ。None ならビン決定用の xMin/xMax
    # (pT は sqrt(xMin)/sqrt(xMax))をそのまま使う。Flat-Stats ビン決定と
    # 独立に、QAプロットだけ広い/狭いレンジで見たい場合に指定する。
    overlayPt2Min=None,
    overlayPt2Max=None,
    overlayPtMin=None,
    overlayPtMax=None,
):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    # --------------------------------------------------------
    # 0. 入力ファイル展開(単一パス/リスト/glob パターンいずれも可)し、
    #    TTree を直接 TChain する。
    #    (.cxx タスクの生出力に直接つなぐので、Step1 相当のマージ
    #    macro/ファイルを経由しない)
    # --------------------------------------------------------
    patterns = [inputFiles] if isinstance(inputFiles, str) else list(inputFiles)
    files = []
    for pattern in patterns:
        matched = sorted(glob.glob(pattern))
        files.extend(matched if matched else [pattern])

    if not files:
        print("[Error] No input files resolved from inputFiles.")
        return

    resolvedTreePath = None
    for fname in files:
        f = ROOT.TFile.Open(fname, "READ")
        if not f or f.IsZombie():
            print(f"[Warning] Cannot open {fname}, skipping.")
            continue
        for cand in (treeName, f"upc-muon-analysis/{treeName}"):
            if f.Get(cand):
                resolvedTreePath = cand
                break
        f.Close()
        if resolvedTreePath:
            break

    if resolvedTreePath is None:
        print(f"[Error] Tree '{treeName}' not found in any of {len(files)} input file(s) "
              f"(tried top-level and 'upc-muon-analysis/{treeName}'). "
              "Make sure processExportOmniFoldTuple was enabled for the upc-muon-analysis task.")
        return

    chain = ROOT.TChain(resolvedTreePath)
    for fname in files:
        chain.Add(fname)

    df = ROOT.RDataFrame(chain)
    nEntries = df.Count().GetValue()
    if nEntries == 0:
        print(f"[Error] Chained tree '{resolvedTreePath}' has 0 entries across {len(files)} file(s).")
        return

    print(f"[Info] Chained {len(files)} file(s) directly from the upc-muon-analysis task output, "
          f"tree path '{resolvedTreePath}', {nEntries} entries total.")

    # --------------------------------------------------------
    # 1. Gen 側の生イベント値(pairPt2MC)を splitMin-splitMax の範囲だけで
    #    ソートし、等統計になる境界値をそのまま quantile として取り出す。
    #    (ヒストグラムの bin を介さないので離散化誤差が入らない)
    #    Coherent が支配的な [xMin, splitMin) は quantile 計算に使わない
    #    (この領域の統計量に等統計分割が引きずられないようにする)。
    # --------------------------------------------------------
    splitMaxResolved = xMax if splitMax is None else splitMax

    genValues = list(df.Filter(f"pairPt2MC >= {splitMin} && pairPt2MC < {splitMaxResolved}")
                       .Take["float"]("pairPt2MC").GetValue())
    genValues.sort()

    n = len(genValues)
    if n == 0:
        print(f"[Error] No Gen entries with pairPt2MC in [{splitMin}, {splitMaxResolved}).")
        return

    step = n / float(nTargetBins)
    splitEdges = [splitMin]
    for k in range(1, nTargetBins):
        idx = min(max(int(round(k * step)), 1), n - 1)
        val = genValues[idx]
        if val > splitEdges[-1] and val < splitMaxResolved:
            splitEdges.append(val)
    if splitEdges[-1] < splitMaxResolved:
        splitEdges.append(splitMaxResolved)

    # [xMin, splitMin) / (splitMaxResolved, xMax] を、分割せず単一 bin として
    # 前後に付け足す(Coherent 領域などを 2D プロットの表示範囲には残しつつ、
    # 等統計分割の対象からは外すため)。
    edges = list(splitEdges)
    if splitMin > xMin:
        edges.insert(0, xMin)
    if splitMaxResolved < xMax:
        edges.append(xMax)

    nBins = len(edges) - 1
    binArray = array.array("d", edges)

    # --------------------------------------------------------
    # 2. Gen / Reco の 1D ヒストグラムを、決定した境界で
    #    TTree のイベントから直接 Fill(疑似ビン center 近似を挟まない)
    # --------------------------------------------------------
    hGenModel = ROOT.RDF.TH1DModel(
        "hGenPt2_flatUnbinned",
        "Gen p_{T}^{2} (Flat Stats, Unbinned);p_{T}^{2} (GeV^{2}/c^{2});Events",
        nBins, binArray)
    hRecoModel = ROOT.RDF.TH1DModel(
        "hRecoPt2_flatUnbinned",
        "Reco p_{T}^{2} (Flat Stats, Unbinned);p_{T}^{2} (GeV^{2}/c^{2});Events",
        nBins, binArray)

    hGenFlat = df.Filter(f"pairPt2MC >= {xMin} && pairPt2MC < {xMax}") \
                 .Histo1D(hGenModel, "pairPt2MC").GetValue()
    hRecoFlat = df.Filter(f"passRecoFlag && pairPt2Reco >= {xMin} && pairPt2Reco < {xMax}") \
                  .Histo1D(hRecoModel, "pairPt2Reco").GetValue()

    print("=" * 40)
    print(f"[Info] Target events per bin (within split range [{splitMin}, {splitMaxResolved})): {step:.1f}")
    print("[Info] Gen Bin Contents after rebinning:")
    for i in range(1, nBins + 1):
        lo, hi = hGenFlat.GetBinLowEdge(i), hGenFlat.GetBinLowEdge(i + 1)
        tag = " (lump, outside split range)" if (lo < splitMin or hi > splitMaxResolved) else ""
        print(f"  Bin {i} [{lo:.6f}, {hi:.6f}] : {hGenFlat.GetBinContent(i):.0f} events{tag}")
    print("=" * 40)

    # --------------------------------------------------------
    # 3. 2D レスポンス行列も、あらかじめ詰め込まれた fine 2D hist を
    #    リビンするのではなく、(pairPt2MC, pairPt2Reco) の event pair を
    #    直接 Fill する(passRecoFlag が立っている = マッチした pair のみ)
    # --------------------------------------------------------
    hMatModel = ROOT.RDF.TH2DModel(
        "hResponseMatrixPairPt2_flatUnbinned",
        "Response Matrix (Flat Stats, Unbinned);Gen p_{T}^{2};Reco p_{T}^{2}",
        nBins, binArray, nBins, binArray)
    hMatFlat = df.Filter(f"passRecoFlag && pairPt2MC >= {xMin} && pairPt2MC < {xMax} "
                          f"&& pairPt2Reco >= {xMin} && pairPt2Reco < {xMax}") \
                 .Histo2D(hMatModel, "pairPt2MC", "pairPt2Reco").GetValue()

    # --------------------------------------------------------
    # 4. Purity / Stability / Migration probability
    #
    #    Purity/Stability は Make2DHistFromBinnedHists.C と同じ定義のまま
    #    (efficiency 込みの Gen/Reco 総数で正規化)。
    #
    #    Migration probability M_ij = P(reco が bin j | gen が bin i, かつ
    #    reco された)は、行 i の合計が 100% になる「真の確率」として定義
    #    し直す。以前は分母に hGenFlat(未検出分も含む全 Gen)を使っていた
    #    ため、行の合計が efficiency 分しかなく確率になっていなかった
    #    (合計が100%にならない)。正しくは、行 i の分母は同じ行 i の
    #    hMatFlat の合計(= reco された分だけ)にする。
    # --------------------------------------------------------
    hPurity = ROOT.TH1D("hPurity_flatUnbinned", "Purity;Gen bin;Purity", nBins, binArray)
    hStability = ROOT.TH1D("hStability_flatUnbinned", "Stability;Reco bin;Stability", nBins, binArray)
    hMigrationProb = ROOT.TH2D(
        "hMigrationProb_flatUnbinned", "Migration Probability P(reco bin | gen bin, reconstructed) [%];Gen p_{T}^{2};Reco p_{T}^{2}",
        nBins, binArray, nBins, binArray)

    print("\n" + "=" * 76)
    print("  Rebinning Result Table (Unbinned)")
    print("=" * 76)
    print(f"{'Bin':<5} | {'X Range (Gen)':<16} | {'Y Range (Reco)':<16} | "
          f"{'Purity%':<8} | {'Stab%':<8} | {'Diag Migration (%)':<18}")
    print("-" * 76)
    for i in range(1, nBins + 1):
        diag = hMatFlat.GetBinContent(i, i)
        sumTrue = hGenFlat.GetBinContent(i)
        sumReco = hRecoFlat.GetBinContent(i)
        rowSumMatched = sum(hMatFlat.GetBinContent(i, j) for j in range(1, nBins + 1))

        purity = diag / sumReco if sumReco > 0 else 0.0
        stability = diag / sumTrue if sumTrue > 0 else 0.0
        hPurity.SetBinContent(i, purity)
        hStability.SetBinContent(i, stability)

        for j in range(1, nBins + 1):
            cellVal = hMatFlat.GetBinContent(i, j)
            prob = (cellVal / rowSumMatched * 100.0) if rowSumMatched > 0 else 0.0
            hMigrationProb.SetBinContent(i, j, prob)

        diagProb = hMigrationProb.GetBinContent(i, i)
        print(f"{i:<5} | {edges[i - 1]:6.4f}-{edges[i]:6.4f}   | "
              f"{edges[i - 1]:6.4f}-{edges[i]:6.4f}   | "
              f"{purity * 100:<8.1f} | {stability * 100:<8.1f} | {diagProb:<18.1f}")
    print("=" * 76 + "\n")

    # --------------------------------------------------------
    # 5. QA: MC True と Reco の重ね書き用ヒストグラム
    #
    #    表示レンジは overlayPt2Min/Max, overlayPtMin/Max で xMin/xMax
    #    (Flat-Stats ビン決定用レンジ)とは独立に指定できる(未指定なら
    #    xMin/xMax およびその sqrt を使う)。ビン決定用の hGenFlat/hRecoFlat
    #    とは別に、この QA レンジ専用の等幅ヒストグラムをここで作る。
    #
    #    Reco は passRecoFlag が立っている pair だけなので、Gen 側も同じ
    #    条件(Reco-matched)で filter した版を作り、そちらを重ね書きに使う。
    #    "全 Gen" と比較すると再構成効率の分だけ形が変わって見えてしまい、
    #    分解能(smearing)だけを見たい QA としては公平な比較にならないため。
    # --------------------------------------------------------
    qaPt2Min = xMin if overlayPt2Min is None else overlayPt2Min
    qaPt2Max = xMax if overlayPt2Max is None else overlayPt2Max
    qaPtMin = math.sqrt(max(qaPt2Min, 0.0)) if overlayPtMin is None else overlayPtMin
    qaPtMax = math.sqrt(qaPt2Max) if overlayPtMax is None else overlayPtMax

    hGenPt2QA = df.Filter(f"pairPt2MC >= {qaPt2Min} && pairPt2MC < {qaPt2Max}") \
                  .Histo1D(("hGenPt2_QA", "Pair p_{T}^{2} MC True;p_{T}^{2} (GeV^{2}/c^{2});Events",
                            nBinsPtQA, qaPt2Min, qaPt2Max), "pairPt2MC").GetValue()
    hRecoPt2QA = df.Filter(f"passRecoFlag && pairPt2Reco >= {qaPt2Min} && pairPt2Reco < {qaPt2Max}") \
                   .Histo1D(("hRecoPt2_QA", "Pair p_{T}^{2} Reco;p_{T}^{2} (GeV^{2}/c^{2});Events",
                             nBinsPtQA, qaPt2Min, qaPt2Max), "pairPt2Reco").GetValue()
    hGenPt2QAMatched = df.Filter(f"passRecoFlag && pairPt2MC >= {qaPt2Min} && pairPt2MC < {qaPt2Max} "
                                 f"&& pairPt2Reco >= {qaPt2Min} && pairPt2Reco < {qaPt2Max}") \
                        .Histo1D(("hGenPt2_QA_matched", "Pair p_{T}^{2} MC True, Reco-matched pairs only;p_{T}^{2} (GeV^{2}/c^{2});Events",
                                  nBinsPtQA, qaPt2Min, qaPt2Max), "pairPt2MC").GetValue()

    hGenPtQA = df.Filter(f"pairPt2MC >= {qaPt2Min} && pairPt2MC < {qaPt2Max}") \
                 .Define("pairPtMC_QA", "sqrt(pairPt2MC)") \
                 .Histo1D(("hGenPt_QA", "Pair p_{T} MC True;p_{T} (GeV/c);Events",
                           nBinsPtQA, qaPtMin, qaPtMax), "pairPtMC_QA").GetValue()
    hRecoPtQA = df.Filter(f"passRecoFlag && pairPt2Reco >= {qaPt2Min} && pairPt2Reco < {qaPt2Max}") \
                  .Define("pairPtReco_QA", "sqrt(pairPt2Reco)") \
                  .Histo1D(("hRecoPt_QA", "Pair p_{T} Reco;p_{T} (GeV/c);Events",
                            nBinsPtQA, qaPtMin, qaPtMax), "pairPtReco_QA").GetValue()
    hGenPtQAMatched = df.Filter(f"passRecoFlag && pairPt2MC >= {qaPt2Min} && pairPt2MC < {qaPt2Max} "
                                f"&& pairPt2Reco >= {qaPt2Min} && pairPt2Reco < {qaPt2Max}") \
                        .Define("pairPtMC_QA_matched", "sqrt(pairPt2MC)") \
                        .Histo1D(("hGenPt_QA_matched", "Pair p_{T} MC True, Reco-matched pairs only;p_{T} (GeV/c);Events",
                                  nBinsPtQA, qaPtMin, qaPtMax), "pairPtMC_QA_matched").GetValue()

    # --------------------------------------------------------
    # Save ROOT
    # --------------------------------------------------------
    outDir = os.path.dirname(outputFile)
    if outDir:
        os.makedirs(outDir, exist_ok=True)
    fOut = ROOT.TFile.Open(outputFile, "RECREATE")
    hGenFlat.Write()
    hRecoFlat.Write()
    hMatFlat.Write()
    hPurity.Write()
    hStability.Write()
    hMigrationProb.Write()
    hGenPt2QA.Write()
    hRecoPt2QA.Write()
    hGenPt2QAMatched.Write()
    hGenPtQA.Write()
    hRecoPtQA.Write()
    hGenPtQAMatched.Write()
    fOut.Close()
    print(f"[Info] Saved ROOT: {outputFile}")

    # --------------------------------------------------------
    # Draw / Save PNG
    # --------------------------------------------------------
    outBase = outputFile.replace(".root", "")

    def drawAndSave1D(h, title, suffix, color):
        c = ROOT.TCanvas("c_" + suffix, title, 800, 600)
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.12)
        c.SetGrid()
        h.SetMinimum(0.0)
        h.SetMaximum(h.GetMaximum() * 1.2)
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.Draw("HIST")

        tex = ROOT.TLatex()
        tex.SetNDC()
        tex.SetTextSize(0.03)
        tex.DrawLatex(0.55, 0.79, "Method: Flat Distributions (Unbinned/TTree)")
        tex.DrawLatex(0.55, 0.74, f"Bins: {nBins}")
        c.SaveAs(outBase + "_" + suffix + ".png")

    drawAndSave1D(hGenFlat, "Gen p_{T}^{2} (Flat Stats, Unbinned)", "Gen", ROOT.kRed + 1)
    drawAndSave1D(hRecoFlat, "Reco p_{T}^{2} (Flat Stats, Unbinned)", "Reco", ROOT.kBlue + 1)

    def drawOverlay(hGen, hReco, title, suffix):
        c = ROOT.TCanvas("c_overlay_" + suffix, title, 800, 600)
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.12)
        c.SetGrid()

        hGen.SetTitle(title)
        hGen.SetLineColor(ROOT.kBlack)
        hGen.SetLineWidth(2)
        hGen.SetMinimum(0.0)
        hGen.SetMaximum(max(hGen.GetMaximum(), hReco.GetMaximum()) * 1.2)
        hGen.Draw("HIST")

        hReco.SetLineColor(ROOT.kRed + 1)
        hReco.SetMarkerColor(ROOT.kRed + 1)
        hReco.SetMarkerStyle(20)
        hReco.SetMarkerSize(0.7)
        hReco.Draw("SAME P E")

        leg = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
        leg.AddEntry(hGen, "MC True (Reco-matched)", "l")
        leg.AddEntry(hReco, "Reco", "pe")
        leg.Draw()

        c.SaveAs(outBase + "_Overlay_" + suffix + ".png")

    # Reco と同じ pair 集合の Gen (passRecoFlag 一致)で重ねる。
    # 全 Gen (hGenPt2QA 等、passRecoFlag 不問版)と比べると効率損失分の
    # 形状差が混ざるため、smearing だけを見たいこの QA には使わない。
    drawOverlay(hGenPt2QAMatched, hRecoPt2QA, "Pair p_{T}^{2}: MC True vs Reco (Reco-matched)", "Pt2")
    drawOverlay(hGenPtQAMatched, hRecoPtQA, "Pair p_{T}: MC True vs Reco (Reco-matched)", "Pt")

    ROOT.gStyle.SetPaintTextFormat(".1f")
    cMat = ROOT.TCanvas("c_Matrix", "Response Matrix (Unbinned)", 800, 700)
    cMat.SetRightMargin(0.15)
    cMat.SetLogz()
    hMatFlat.Draw("COLZ")
    cMat.SaveAs(outBase + "_Matrix.png")

    cMigration = ROOT.TCanvas("c_Migration", "Migration Probability [%]", 800, 700)
    cMigration.SetRightMargin(0.15)
    cMigration.SetLeftMargin(0.15)
    hMigrationProb.Draw("COLZ TEXT")
    cMigration.SaveAs(outBase + "_MigrationProbability.png")

    cPS = ROOT.TCanvas("c_PS", "Purity & Stability", 900, 450)
    cPS.Divide(2, 1)
    cPS.cd(1)
    hPurity.GetYaxis().SetRangeUser(0, 1)
    hPurity.Draw("HIST")
    cPS.cd(2)
    hStability.GetYaxis().SetRangeUser(0, 1)
    hStability.Draw("HIST")
    cPS.SaveAs(outBase + "_PurityStability.png")

    print(f"[Info] Saved PNGs: {outBase}_Gen.png / _Reco.png / _Matrix.png / "
          f"_MigrationProbability.png / _PurityStability.png / "
          f"_Overlay_Pt2.png / _Overlay_Pt.png")


if __name__ == "__main__":
    RebinFlatStatsUnbinned()
