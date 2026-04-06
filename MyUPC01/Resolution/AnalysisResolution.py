#!/usr/bin/env python3
"""
AnalysisResolution.py
---------------------
pT Resolution analysis macro for forward muons in UPC events.

Reads AnalysisResults.root (output of o2-analysis-my-upc-muon-resolution) and:
  1. Slices hPtResoVsPtTrue along the X-axis (pT_true) in configurable bin widths
  2. Fits each Y-axis projection with a Gaussian to extract sigma
  3. Plots sigma vs pT_true (resolution curve)
  4. Also draws the response matrix hResponseMatrix

Usage:
  python3 AnalysisResolution.py [AnalysisResults.root]

If no argument is given, it looks for AnalysisResults.root in the current directory.
"""

import sys
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetPalette(ROOT.kBird)

# ── Configuration ──────────────────────────────────────────────────────────────
INPUT_FILE   = sys.argv[1] if len(sys.argv) > 1 else "AnalysisResults.root"
TASK_DIR     = "registry"          # HistogramRegistry name in the cxx task
HIST_RESO    = "hPtResoVsPtTrue"
HIST_RESP    = "hResponseMatrix"

# Slice width in GeV/c along the pT_true axis
SLICE_WIDTH  = 0.5   # GeV

# Gaussian fit range in units of sigma (±N*sigma after first rough fit)
FIT_NSIGMA   = 2.0

OUT_RESO     = "resolution_plot.pdf"
OUT_RESP     = "response_matrix.pdf"
# ──────────────────────────────────────────────────────────────────────────────


def open_histogram(f, task_dir, hist_name):
    """Retrieve a histogram from the ROOT file, navigating into task sub-directories."""
    # Try direct path first
    h = f.Get(f"{task_dir}/{hist_name}")
    if h:
        h.SetDirectory(0)
        return h
    # Some O2 outputs wrap things in an extra directory
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if hasattr(obj, "Get"):
            h = obj.Get(hist_name)
            if h:
                h.SetDirectory(0)
                return h
    # Last resort: recursive search
    h = f.Get(hist_name)
    if h:
        h.SetDirectory(0)
        return h
    return None


def fit_projection(h1, ptLow, ptHigh):
    """
    Fit a 1D histogram (Y-projection) with a Gaussian.
    Returns (mean, sigma, sigma_err) or None on failure.
    """
    if h1.GetEntries() < 10:
        return None

    # First rough fit over full range
    fit0 = ROOT.TF1("gaus0", "gaus", h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
    h1.Fit(fit0, "Q0")

    mean0  = fit0.GetParameter(1)
    sigma0 = fit0.GetParameter(2)
    if sigma0 <= 0:
        return None

    # Refined fit within ±FIT_NSIGMA * sigma
    lo = mean0 - FIT_NSIGMA * sigma0
    hi = mean0 + FIT_NSIGMA * sigma0
    fit1 = ROOT.TF1(f"gaus_{ptLow:.2f}_{ptHigh:.2f}", "gaus", lo, hi)
    result = h1.Fit(fit1, "QRS")

    sigma     = fit1.GetParameter(2)
    sigma_err = fit1.GetParError(2)
    mean      = fit1.GetParameter(1)

    if sigma <= 0:
        return None
    return (mean, sigma, sigma_err)


def make_resolution_plot(h2d):
    """
    Slice h2d along the X axis and fit each Y projection.
    Returns a TGraphErrors(sigma_vs_ptTrue).
    """
    xAxis  = h2d.GetXaxis()
    xMin   = xAxis.GetXmin()
    xMax   = xAxis.GetXmax()
    nBinsX = xAxis.GetNbins()

    # Build a list of (binLo, binHi) groups corresponding to SLICE_WIDTH
    pt_centers = []
    pt_errors  = []
    sigmas     = []
    sigma_errs = []

    pt = xMin + SLICE_WIDTH / 2.0
    while pt < xMax:
        ptLow  = pt - SLICE_WIDTH / 2.0
        ptHigh = pt + SLICE_WIDTH / 2.0

        binLo = xAxis.FindBin(ptLow  + 1e-6)
        binHi = xAxis.FindBin(ptHigh - 1e-6)
        if binLo > nBinsX or binHi < 1:
            pt += SLICE_WIDTH
            continue

        binLo = max(binLo, 1)
        binHi = min(binHi, nBinsX)

        h1 = h2d.ProjectionY(f"proj_{ptLow:.2f}", binLo, binHi)
        h1.SetDirectory(0)

        res = fit_projection(h1, ptLow, ptHigh)
        if res is not None:
            _, sigma, sigma_err = res
            pt_centers.append((ptLow + ptHigh) / 2.0)
            pt_errors.append(SLICE_WIDTH / 2.0)
            sigmas.append(sigma)
            sigma_errs.append(sigma_err)

        pt += SLICE_WIDTH

    n = len(pt_centers)
    if n == 0:
        return None

    gr = ROOT.TGraphErrors(n)
    gr.SetName("grResolution")
    gr.SetTitle("pT Resolution;#it{p}_{T}^{true} (GeV/#it{c});#sigma [(#it{p}_{T}^{reco}-#it{p}_{T}^{true})/#it{p}_{T}^{true}]")
    for i in range(n):
        gr.SetPoint(i, pt_centers[i], sigmas[i])
        gr.SetPointError(i, pt_errors[i], sigma_errs[i])

    return gr


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    if not os.path.isfile(INPUT_FILE):
        print(f"[ERROR] Input file not found: {INPUT_FILE}")
        sys.exit(1)

    f = ROOT.TFile.Open(INPUT_FILE, "READ")
    if not f or f.IsZombie():
        print(f"[ERROR] Cannot open ROOT file: {INPUT_FILE}")
        sys.exit(1)

    # ── Load histograms ────────────────────────────────────────────────────────
    h_reso = open_histogram(f, TASK_DIR, HIST_RESO)
    h_resp = open_histogram(f, TASK_DIR, HIST_RESP)

    if not h_reso:
        print(f"[ERROR] Histogram '{HIST_RESO}' not found in {INPUT_FILE}")
        f.Close()
        sys.exit(1)
    if not h_resp:
        print(f"[WARNING] Histogram '{HIST_RESP}' not found — skipping response matrix plot")

    # ── Resolution plot ────────────────────────────────────────────────────────
    print(f"[INFO] Building resolution curve from '{HIST_RESO}' ...")
    gr = make_resolution_plot(h_reso)

    c1 = ROOT.TCanvas("c1", "pT Resolution", 900, 700)
    c1.SetLeftMargin(0.14)
    c1.SetBottomMargin(0.13)

    # Draw 2D first (as color map in top pad)
    c1.Divide(1, 2)

    pad1 = c1.cd(1)
    pad1.SetLeftMargin(0.14)
    h_reso_draw = h_reso.Clone("hResoForDraw")
    h_reso_draw.SetTitle("Relative pT Error vs pT^{true}")
    h_reso_draw.Draw("COLZ")

    pad2 = c1.cd(2)
    pad2.SetLeftMargin(0.14)
    if gr and gr.GetN() > 0:
        gr.SetMarkerStyle(20)
        gr.SetMarkerColor(ROOT.kBlue + 1)
        gr.SetLineColor(ROOT.kBlue + 1)
        gr.SetMarkerSize(1.1)
        gr.Draw("APE")
        # Overlay a smooth line fit (pol2)
        fit_pol = ROOT.TF1("fitPol", "pol2", gr.GetXaxis().GetXmin(), gr.GetXaxis().GetXmax())
        fit_pol.SetLineColor(ROOT.kRed)
        fit_pol.SetLineWidth(2)
        gr.Fit(fit_pol, "Q")
    else:
        print("[WARNING] No valid resolution points — check histogram statistics")

    c1.SaveAs(OUT_RESO)
    print(f"[INFO] Saved resolution plot → {OUT_RESO}")

    # ── Response matrix ────────────────────────────────────────────────────────
    if h_resp:
        print(f"[INFO] Drawing response matrix '{HIST_RESP}' ...")
        c2 = ROOT.TCanvas("c2", "Response Matrix", 800, 700)
        c2.SetLeftMargin(0.14)
        c2.SetBottomMargin(0.13)
        c2.SetRightMargin(0.15)
        h_resp.SetTitle("Response Matrix")
        h_resp.Draw("COLZ")
        # Draw diagonal reference line
        xmin = h_resp.GetXaxis().GetXmin()
        xmax = h_resp.GetXaxis().GetXmax()
        diag = ROOT.TLine(xmin, xmin, xmax, xmax)
        diag.SetLineColor(ROOT.kRed)
        diag.SetLineWidth(2)
        diag.SetLineStyle(2)
        diag.Draw("SAME")
        c2.SaveAs(OUT_RESP)
        print(f"[INFO] Saved response matrix  → {OUT_RESP}")

    f.Close()
    print("[INFO] Done.")


if __name__ == "__main__":
    main()
