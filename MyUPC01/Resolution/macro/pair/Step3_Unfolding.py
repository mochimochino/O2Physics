import ROOT
import sys
import math

try:
    ROOT.gSystem.Load("libRooUnfold")
except Exception as e:
    print(f"[Error] Failed to load RooUnfold: {e}")
    sys.exit(1)

def Step3_Unfolding():
    # --------------------------------------------
    # Setting
    # -------------------------------------------
    inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0408/Resolution/"
    inFile = inDir + "Step2_Response_for_Unfolding.root"
    outDir = inDir

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFontSize(0.045)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)

    fIn = ROOT.TFile.Open(inFile, "READ")
    if not fIn or fIn.IsZombie():
        print(f"[Error] Cannot open {inFile}")
        return

    hGen = fIn.Get("hMCGen_Rebinned")
    hReco = fIn.Get("hMCReco_Rebinned")
    hMat = fIn.Get("hResponseMatrix")
    hData = fIn.Get("hDataReco_Rebinned")

    if not all([hGen, hReco, hMat, hData]):
        print(f"[Error] Required histograms not found in {inFile}")
        return

    print("[Info] Constructing RooUnfoldResponse...")
    response = ROOT.RooUnfoldResponse(hReco, hGen, hMat)

    # ----------------------------------------------------------
    # Unfolding (Bayes & SVD)
    # ----------------------------------------------------------
    nIterations = 4
    unfoldBayes = ROOT.RooUnfoldBayes(response, hData, nIterations)
    hUnfoldedBayes = unfoldBayes.Hreco().Clone("hUnfoldedBayes")
    hUnfoldedBayes.SetLineColor(ROOT.kRed)
    hUnfoldedBayes.SetMarkerColor(ROOT.kRed)
    hUnfoldedBayes.SetMarkerStyle(20)

    kTerm = int(hGen.GetNbinsX() / 2)
    if kTerm < 2: kTerm = 2
    unfoldSvd = ROOT.RooUnfoldSvd(response, hData, kTerm)
    hUnfoldedSvd = unfoldSvd.Hreco().Clone("hUnfoldedSvd")
    hUnfoldedSvd.SetLineColor(ROOT.kBlue)
    hUnfoldedSvd.SetMarkerColor(ROOT.kBlue)
    hUnfoldedSvd.SetMarkerStyle(21)

    # ----------------------------------------------------------
    # 1. Closure Test Plot
    # ----------------------------------------------------------
    print("[Info] Plotting Closure Test...")
    hGen.SetLineColor(ROOT.kBlack)
    hGen.SetLineWidth(2)

    c1 = ROOT.TCanvas("c1", "Closure Test", 800, 800)
    c1.Divide(1, 2)

    pad1 = c1.GetPad(1)
    pad1.SetPad(0.0, 0.3, 1.0, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    hGen.SetTitle("Closure Test: Unfolded vs True MC")
    hGen.GetXaxis().SetLabelSize(0)
    hGen.GetXaxis().SetTitleSize(0)
    hGen.Draw("HIST")
    hUnfoldedBayes.Draw("PE SAME")
    hUnfoldedSvd.Draw("PE SAME")

    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.85)
    leg.AddEntry(hGen, "True (MC Gen)", "l")
    leg.AddEntry(hUnfoldedBayes, f"Bayes (iter={nIterations})", "pe")
    leg.AddEntry(hUnfoldedSvd, f"SVD (k={kTerm})", "pe")
    leg.Draw()

    c1.cd()
    pad2 = c1.GetPad(2)
    pad2.SetPad(0.0, 0.0, 1.0, 0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()

    hRatioBayes = hUnfoldedBayes.Clone("hRatioBayes")
    hRatioBayes.Divide(hGen)
    hRatioBayes.SetTitle("")
    hRatioBayes.GetYaxis().SetTitle("Unfolded / True")
    hRatioBayes.GetYaxis().SetRangeUser(0.1, 1.9)
    hRatioBayes.GetYaxis().SetNdivisions(505)
    hRatioBayes.GetYaxis().SetLabelSize(0.1)
    hRatioBayes.GetYaxis().SetTitleSize(0.12)
    hRatioBayes.GetYaxis().SetTitleOffset(0.4)
    hRatioBayes.GetXaxis().SetTitle("p_{T}^{2} (GeV^{2}/c^{2})")
    hRatioBayes.GetXaxis().SetLabelSize(0.1)
    hRatioBayes.GetXaxis().SetTitleSize(0.12)
    hRatioBayes.Draw("PE")

    hRatioSvd = hUnfoldedSvd.Clone("hRatioSvd")
    hRatioSvd.Divide(hGen)
    hRatioSvd.Draw("PE SAME")

    line = ROOT.TLine(hGen.GetXaxis().GetXmin(), 1.0, hGen.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw("SAME")
    c1.SaveAs(outDir + "Step3_ClosureTest_PyROOT.png")

    # ----------------------------------------------------------
    # 2. Bayesian Iteration Study
    # ----------------------------------------------------------
    print("[Info] Running Bayesian Iteration Study...")
    c2 = ROOT.TCanvas("c2", "Bayes Iteration Study", 800, 800)
    c2.Divide(1, 2)

    pad2_1 = c2.GetPad(1)
    pad2_1.SetPad(0.0, 0.3, 1.0, 1.0)
    pad2_1.SetBottomMargin(0.02)
    pad2_1.SetLogy()
    pad2_1.Draw()
    pad2_1.cd()

    hGenIter = hGen.Clone("hGenIter")
    hGenIter.SetTitle("Bayesian Unfolding: Convergence Study")
    hGenIter.Draw("HIST")

    legIter = ROOT.TLegend(0.65, 0.55, 0.85, 0.85)
    legIter.AddEntry(hGenIter, "True (MC Gen)", "l")

    maxIters = 5
    colors = [ROOT.kRed, ROOT.kOrange + 1, ROOT.kGreen + 2, ROOT.kAzure + 1, ROOT.kMagenta]
    markers = [20, 21, 22, 23, 33]
    vecIterHists = []
    vecRatHists = []

    for i in range(1, maxIters + 1):
        unfoldTest = ROOT.RooUnfoldBayes(response, hData, i)
        hUnf = unfoldTest.Hreco().Clone(f"hUnf_iter{i}")
        hUnf.SetLineColor(colors[i - 1])
        hUnf.SetMarkerColor(colors[i - 1])
        hUnf.SetMarkerStyle(markers[i - 1])
        hUnf.Draw("PE SAME")
        legIter.AddEntry(hUnf, f"Iter = {i}", "pe")
        vecIterHists.append(hUnf)

    legIter.Draw()

    c2.cd()
    pad2_2 = c2.GetPad(2)
    pad2_2.SetPad(0.0, 0.0, 1.0, 0.3)
    pad2_2.SetTopMargin(0.02)
    pad2_2.SetBottomMargin(0.3)
    pad2_2.Draw()
    pad2_2.cd()

    for i in range(maxIters):
        hRatIter = vecIterHists[i].Clone(f"hRatIter_{i}")
        hRatIter.Divide(hGen)
        if i == 0:
            hRatIter.SetTitle("")
            hRatIter.GetYaxis().SetTitle("Unfolded / True")
            hRatIter.GetYaxis().SetRangeUser(0.5, 1.5)
            hRatIter.GetYaxis().SetNdivisions(505)
            hRatIter.GetYaxis().SetLabelSize(0.1)
            hRatIter.GetYaxis().SetTitleSize(0.12)
            hRatIter.GetYaxis().SetTitleOffset(0.4)
            hRatIter.GetXaxis().SetTitle("p_{T}^{2} (GeV^{2}/c^{2})")
            hRatIter.GetXaxis().SetLabelSize(0.1)
            hRatIter.GetXaxis().SetTitleSize(0.12)
            hRatIter.Draw("PE")
        else:
            hRatIter.Draw("PE SAME")
        vecRatHists.append(hRatIter) # Keep objects alive
        
    line.Draw("SAME")
    c2.SaveAs(outDir + "Step3_IterationStudy_PyROOT.png")

    # ----------------------------------------------------------
    # 3. Binning Appropriateness Check: Correlation Matrix
    # ----------------------------------------------------------
    print("[Info] Extracting Covariance Matrix to evaluate binning...")


    covMatrix = unfoldBayes.Ereco(ROOT.RooUnfold.kCovariance)
    nBins = hGen.GetNbinsX()
    
    hCorr = ROOT.TH2D("hCorr", "Correlation Matrix (Bayesian);Bin i;Bin j", nBins, 0.5, nBins+0.5, nBins, 0.5, nBins+0.5)
    
    for i in range(nBins):
        for j in range(nBins):
            var_i = covMatrix(i, i)
            var_j = covMatrix(j, j)
            if var_i > 0 and var_j > 0:
                # Correlation Matrix = Cov(i,j) / sqrt(V(i)*V(j))
                corr = covMatrix(i, j) / math.sqrt(var_i * var_j)
                hCorr.SetBinContent(i + 1, j + 1, corr)
            else:
                hCorr.SetBinContent(i + 1, j + 1, 0)

    # Draw
    c3 = ROOT.TCanvas("c3", "Correlation Matrix", 800, 600)
    c3.SetRightMargin(0.15)
    hCorr.GetZaxis().SetRangeUser(-1.0, 1.0)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    hCorr.Draw("COLZ TEXT")
    c3.SaveAs(outDir + "Step3_CorrelationMatrix_PyROOT.png")
    
    print("[Done] All plots generated successfully.")

if __name__ == "__main__":
    Step3_Unfolding()
