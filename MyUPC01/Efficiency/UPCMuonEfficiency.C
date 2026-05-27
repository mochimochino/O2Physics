/// \file   UPCMuonEfficiency.C
/// \brief  Efficiency analysis macro for UPC dimuon photoproduction.
///         Reads the output of UPCMuonEfficiencyTask (AnalysisResults.root),
///         computes Eff = Reco / MC for single muons and dimuon pairs
///         (Phi, Eta/y, pT), and saves plots and histograms to a ROOT file.
///
/// Output ROOT file structure:
///   MC distributions    (6 TCanvas): cMC_TrkPhi/Eta/Pt, cMC_PairPhi/Rap/Pt
///   Reco distributions  (6 TCanvas): cReco_TrkPhi/Eta/Pt, cReco_PairPhi/Rap/Pt
///   MC+Reco overlays    (6 TCanvas): cOvl_TrkPhi/Eta/Pt, cOvl_PairPhi/Rap/Pt
///   Efficiency          (6 TCanvas): cEff_TrkPhi/Eta/Pt, cEff_PairPhi/Rap/Pt
///   Source histograms  (12 TH1)
///   Efficiency hists    (6 TH1D)
///
/// \usage  root -l -b -q 'UPCMuonEfficiency.C("AnalysisResults.root")'
/// \author Takuma Matsumoto

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

// ============================================================================
// Global style constants (matching Step1/Step2 macro style)
// ============================================================================
static const int kFont = 42;        // Helvetica
static const double kTitSz = 0.045; // axis title size
static const double kLabSz = 0.040; // axis label size
static const double kTitOff = 1.3;  // X-axis title offset
static const double kTitOffY = 1.2; // Y-axis title offset

// ============================================================================
// Apply global gStyle (call once in main)
// ============================================================================
void ApplyStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  for (const char* ax : {"X", "Y", "Z"}) {
    gStyle->SetLabelFont(kFont, ax);
    gStyle->SetTitleFont(kFont, ax);
    gStyle->SetTitleSize(kTitSz, ax);
    gStyle->SetLabelSize(kLabSz, ax);
  }
}

// ============================================================================
// Helper: apply standard canvas margins (matching Step1 draw2D)
// ============================================================================
void SetCanvasStyle(TCanvas* c, bool hasRightMargin = false)
{
  c->SetLeftMargin(0.11);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(hasRightMargin ? 0.14 : 0.05);
  c->SetTopMargin(0.08);
  c->SetGrid();
}

// ============================================================================
// Helper: apply standard axis style to a TH1
// ============================================================================
void StyleAxis(TH1* h, const char* xTitle, const char* yTitle)
{
  if (!h)
    return;
  h->GetXaxis()->SetTitle(xTitle);
  h->GetXaxis()->SetTitleFont(kFont);
  h->GetXaxis()->SetLabelFont(kFont);
  h->GetXaxis()->SetTitleSize(kTitSz);
  h->GetXaxis()->SetLabelSize(kLabSz);
  h->GetXaxis()->SetTitleOffset(kTitOff);

  h->GetYaxis()->SetTitle(yTitle);
  h->GetYaxis()->SetTitleFont(kFont);
  h->GetYaxis()->SetLabelFont(kFont);
  h->GetYaxis()->SetTitleSize(kTitSz);
  h->GetYaxis()->SetLabelSize(kLabSz);
  h->GetYaxis()->SetTitleOffset(kTitOffY);
}

// ============================================================================
// Helper: draw a TLatex label block (upper-right corner, like Step1)
// ============================================================================
void DrawInfoBlock(const char* titleLine,
                   double x = 0.55, double yStart = 0.86,
                   double dy = 0.055, double sz = 0.034,
                   const char* entriesLine = "")
{
  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(kFont);
  tex.SetTextSize(sz);
  tex.DrawLatex(x, yStart, Form("#bf{%s}", titleLine));
  tex.DrawLatex(x, yStart - dy, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");
  tex.DrawLatex(x, yStart - 2 * dy, "Global Muon (MCH-MID-MFT)");
  if (entriesLine && strlen(entriesLine) > 0)
    tex.DrawLatex(x, yStart - 3 * dy, entriesLine);
}

// ============================================================================
// Helper: fetch a histogram from the O2 task output directory.
// ============================================================================
TH1* GetHist(TFile* f, const char* histName,
             const char* task = "upc-muon-efficiency",
             const char* reg = "registry")
{
  TString path = TString::Format("%s/%s/%s", task, reg, histName);
  TH1* h = dynamic_cast<TH1*>(f->Get(path));
  if (!h) {
    ::Warning("GetHist", "Histogram not found: %s", path.Data());
  }
  return h;
}

// ============================================================================
// Helper: compute Eff = hReco / hMC (bin-by-bin, binomial error).
// ============================================================================
TH1D* MakeEfficiency(TH1* hReco, TH1* hMC, const char* name, const char* title)
{
  if (!hReco || !hMC)
    return nullptr;

  TH1D* hEff = new TH1D(name, title,
                        hReco->GetNbinsX(),
                        hReco->GetXaxis()->GetXmin(),
                        hReco->GetXaxis()->GetXmax());

  for (int b = 1; b <= hReco->GetNbinsX(); ++b) {
    double numer = hReco->GetBinContent(b);
    double denom = hMC->GetBinContent(b);
    if (denom <= 0.) {
      hEff->SetBinContent(b, 0.);
      hEff->SetBinError(b, 0.);
      continue;
    }
    double eff = numer / denom;
    double err = (eff <= 1.) ? TMath::Sqrt(eff * (1. - eff) / denom) : 0.;
    hEff->SetBinContent(b, eff);
    hEff->SetBinError(b, err);
  }
  return hEff;
}

// ============================================================================
// Build a single-histogram canvas (MC only or Reco only)
// ============================================================================
TCanvas* MakeSingleCanvas(const char* name, const char* canTitle,
                          TH1* h, int color,
                          const char* xTitle, const char* yTitle,
                          const char* plotLabel)
{
  TCanvas* c = new TCanvas(name, canTitle, 900, 700);
  SetCanvasStyle(c);
  if (!h)
    return c;

  h->SetLineColor(color);
  h->SetLineWidth(2);
  StyleAxis(h, xTitle, yTitle);
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum() * 1.5);
  h->Draw("HIST");

  TString entriesStr = Form("Entries: #bf{%.0f}", h->GetEntries());
  DrawInfoBlock(plotLabel, 0.55, 0.86, 0.055, 0.034, entriesStr.Data());
  return c;
}

// ============================================================================
// Build an overlay canvas (MC + Reco)
// ============================================================================
TCanvas* MakeOverlayCanvas(const char* name, const char* canTitle,
                           TH1* hMC, TH1* hReco,
                           const char* xTitle, const char* yTitle,
                           const char* plotLabel)
{
  TCanvas* c = new TCanvas(name, canTitle, 900, 700);
  SetCanvasStyle(c);

  if (hMC) {
    hMC->SetLineColor(kBlue + 1);
    hMC->SetLineWidth(2);
    StyleAxis(hMC, xTitle, yTitle);
    if (hReco) {
      double ymax = TMath::Max(hMC->GetMaximum(), hReco->GetMaximum()) * 1.5;
      hMC->GetYaxis()->SetRangeUser(0., ymax);
    } else {
      hMC->GetYaxis()->SetRangeUser(0., hMC->GetMaximum() * 1.5);
    }
    hMC->Draw("HIST");
  }
  if (hReco) {
    hReco->SetLineColor(kRed + 1);
    hReco->SetLineWidth(2);
    StyleAxis(hReco, xTitle, yTitle);
    hReco->Draw("HIST SAME");
  }

  // Legend
  TLegend* leg = new TLegend(0.13, 0.72, 0.42, 0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(kFont);
  leg->SetTextSize(kTitSz);
  if (hMC)
    leg->AddEntry(hMC, "#bf{MC truth}", "l");
  if (hReco)
    leg->AddEntry(hReco, "#bf{Reco}", "l");
  leg->Draw();

  DrawInfoBlock(plotLabel);

  // Event counts for MC and Reco (drawn separately in respective colors)
  {
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(kFont);
    tex.SetTextSize(0.034);
    double xE = 0.55, yE = 0.86 - 3 * 0.055;
    if (hMC)
      tex.DrawLatex(xE, yE,
                    Form("#color[%d]{MC Entries: #bf{%.0f}}",
                         kBlue + 1, hMC->GetEntries()));
    if (hReco)
      tex.DrawLatex(xE, yE - 0.055,
                    Form("#color[%d]{Reco Entries: #bf{%.0f}}",
                         kRed + 1, hReco->GetEntries()));
  }
  return c;
}

// ============================================================================
// Build an efficiency canvas (single panel)
// ============================================================================
TCanvas* MakeEffCanvas(const char* name, const char* canTitle,
                       TH1D* hEff,
                       const char* xTitle,
                       const char* plotLabel)
{
  TCanvas* c = new TCanvas(name, canTitle, 900, 700);
  SetCanvasStyle(c);
  if (!hEff)
    return c;

  hEff->SetMarkerStyle(20);
  hEff->SetMarkerSize(0.9);
  hEff->SetMarkerColor(kBlack);
  hEff->SetLineColor(kBlack);
  hEff->SetLineWidth(2);
  StyleAxis(hEff, xTitle, "Efficiency  (Reco / MC)");
  double effYmax = (hEff->GetMaximum() > 0.) ? hEff->GetMaximum() * 1.5 : 1.5;
  hEff->GetYaxis()->SetRangeUser(0., effYmax);
  hEff->GetYaxis()->SetNdivisions(505);
  hEff->Draw("E1");

  // Unity line
  TLine* line = new TLine(hEff->GetXaxis()->GetXmin(), 1.,
                          hEff->GetXaxis()->GetXmax(), 1.);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 1);
  line->SetLineWidth(2);
  line->Draw("SAME");

  DrawInfoBlock(plotLabel, 0.55, 0.86, 0.055, 0.034);
  return c;
}

// ============================================================================
// Helper: write a canvas to file and PDF
// ============================================================================
void SaveCanvas(TCanvas* c, TFile* fout, const char* pdfFile,
                bool isFirst = false, bool isLast = false)
{
  if (!c)
    return;
  if (fout && fout->IsOpen()) {
    fout->cd();
    c->Write();
  }
  TString pdfArg = pdfFile;
  if (isFirst)
    pdfArg = TString::Format("%s(", pdfFile);
  else if (isLast)
    pdfArg = TString::Format("%s)", pdfFile);
  c->Print(pdfArg);
}

// ============================================================================
// Main macro
// ============================================================================
void UPCMuonEfficiency(const char* inFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/efficiency/0527/Efficiency/040/Incoherent/jpsi-incoh.root",
                       const char* outFile = "UPCMuonEfficiencyIncoherent.root",
                       const char* pdfFile = "UPCMuonEfficiencyIncoherent.pdf")
{
  ApplyStyle();

  // --- Open input file ---
  TFile* fin = TFile::Open(inFile, "READ");
  if (!fin || fin->IsZombie()) {
    ::Error("UPCMuonEfficiency", "Cannot open input file: %s", inFile);
    return;
  }
  ::Info("UPCMuonEfficiency", "Opened: %s", inFile);

  // =========================================================================
  // Fetch source histograms
  // =========================================================================
  TH1* hTrkPhiMC = GetHist(fin, "hTrkPhiMC");
  TH1* hTrkPhiReco = GetHist(fin, "hTrkPhiReco");
  TH1* hTrkEtaMC = GetHist(fin, "hTrkEtaMC");
  TH1* hTrkEtaReco = GetHist(fin, "hTrkEtaReco");
  TH1* hTrkPtMC = GetHist(fin, "hTrkPtMC");
  TH1* hTrkPtReco = GetHist(fin, "hTrkPtReco");

  TH1* hPairPhiMC = GetHist(fin, "hPairPhiMC");
  TH1* hPairPhiReco = GetHist(fin, "hPairPhiReco");
  TH1* hPairRapidityMC = GetHist(fin, "hPairRapidityMC");
  TH1* hPairRapidityReco = GetHist(fin, "hPairRapidityReco");
  TH1* hPairPtMC = GetHist(fin, "hPairPtMC");
  TH1* hPairPtReco = GetHist(fin, "hPairPtReco");

  TH1* hCutFlowMC = GetHist(fin, "hCutFlowMC");
  TH1* hCutFlowReco = GetHist(fin, "hCutFlowReco");

  // =========================================================================
  // Compute efficiencies (Reco / MC)
  // =========================================================================
  TH1D* hEffTrkPhi = MakeEfficiency(hTrkPhiReco, hTrkPhiMC, "hEffTrkPhi", "Single Muon Eff vs #phi");
  TH1D* hEffTrkEta = MakeEfficiency(hTrkEtaReco, hTrkEtaMC, "hEffTrkEta", "Single Muon Eff vs #eta");
  TH1D* hEffTrkPt = MakeEfficiency(hTrkPtReco, hTrkPtMC, "hEffTrkPt", "Single Muon Eff vs p_{T}");
  TH1D* hEffPairPhi = MakeEfficiency(hPairPhiReco, hPairPhiMC, "hEffPairPhi", "Dimuon Pair Eff vs #phi");
  TH1D* hEffPairRap = MakeEfficiency(hPairRapidityReco, hPairRapidityMC, "hEffPairRap", "Dimuon Pair Eff vs y");
  TH1D* hEffPairPt = MakeEfficiency(hPairPtReco, hPairPtMC, "hEffPairPt", "Dimuon Pair Eff vs p_{T}");

  // =========================================================================
  // Open output ROOT file
  // =========================================================================
  TFile* fout = TFile::Open(outFile, "RECREATE");

  // =========================================================================
  // Variable table for loop-driven canvas creation
  // =========================================================================
  struct VarEntry {
    TH1* hMC;
    TH1* hReco;
    TH1D* hEff;
    const char* xTitle;    // axis X label
    const char* yTitle;    // axis Y label ("Counts")
    const char* plotLabel; // label shown via TLatex on the plot
    const char* tag;       // short name for canvas naming
  };

  VarEntry vars[] = {
    {hTrkPhiMC, hTrkPhiReco, hEffTrkPhi,
     "#phi (rad)", "Counts", "Single Muon: #phi", "TrkPhi"},
    {hTrkEtaMC, hTrkEtaReco, hEffTrkEta,
     "#eta", "Counts", "Single Muon: #eta", "TrkEta"},
    {hTrkPtMC, hTrkPtReco, hEffTrkPt,
     "p_{T} (GeV/c)", "Counts", "Single Muon: p_{T}", "TrkPt"},
    {hPairPhiMC, hPairPhiReco, hEffPairPhi,
     "#phi (rad)", "Counts", "Dimuon Pair: #phi", "PairPhi"},
    {hPairRapidityMC, hPairRapidityReco, hEffPairRap,
     "y", "Counts", "Dimuon Pair: y", "PairRap"},
    {hPairPtMC, hPairPtReco, hEffPairPt,
     "p_{T} (GeV/c)", "Counts", "Dimuon Pair: p_{T}", "PairPt"},
  };
  const int nVars = 6;

  // =========================================================================
  // Page 1: CutFlow (opens PDF)
  // =========================================================================
  {
    TCanvas* cCF = new TCanvas("cCutFlow", "CutFlow", 1200, 550);
    cCF->Divide(2, 1);

    auto drawCF = [&](int pad, TH1* h, const char* lbl, int color) {
      cCF->cd(pad);
      gPad->SetLeftMargin(0.11);
      gPad->SetBottomMargin(0.26);
      gPad->SetRightMargin(0.05);
      gPad->SetTopMargin(0.08);
      gPad->SetGrid();
      if (!h)
        return;
      h->SetLineColor(color);
      h->SetFillColorAlpha(color, 0.18);
      h->SetLineWidth(2);
      h->GetXaxis()->SetLabelFont(kFont);
      h->GetXaxis()->SetLabelSize(0.036);
      h->GetXaxis()->LabelsOption("v");
      h->GetYaxis()->SetTitle("Counts");
      h->GetYaxis()->SetTitleFont(kFont);
      h->GetYaxis()->SetTitleSize(kTitSz);
      h->GetYaxis()->SetLabelFont(kFont);
      h->GetYaxis()->SetLabelSize(kLabSz);
      h->GetYaxis()->SetTitleOffset(kTitOffY);
      h->GetYaxis()->SetRangeUser(0., h->GetMaximum() * 1.5);
      h->Draw("HIST");
      TLatex tex;
      tex.SetNDC();
      tex.SetTextFont(kFont);
      tex.SetTextSize(0.040);
      tex.DrawLatex(0.13, 0.92, Form("#bf{%s}", lbl));
    };

    drawCF(1, hCutFlowMC, "CutFlow  MC truth", kBlue + 1);
    drawCF(2, hCutFlowReco, "CutFlow  Reco", kRed + 1);

    SaveCanvas(cCF, fout, pdfFile, /*isFirst=*/true, /*isLast=*/false);
    delete cCF;
  }

  // =========================================================================
  // Pages 2-7: MC-only distributions
  // =========================================================================
  for (int i = 0; i < nVars; ++i) {
    auto& v = vars[i];
    TCanvas* c = MakeSingleCanvas(
      Form("cMC_%s", v.tag),
      Form("MC: %s", v.plotLabel),
      v.hMC, kBlue + 1,
      v.xTitle, v.yTitle,
      Form("MC truth:  %s", v.plotLabel));
    SaveCanvas(c, fout, pdfFile, false, false);
    delete c;
  }

  // =========================================================================
  // Pages 8-13: Reco-only distributions
  // =========================================================================
  for (int i = 0; i < nVars; ++i) {
    auto& v = vars[i];
    TCanvas* c = MakeSingleCanvas(
      Form("cReco_%s", v.tag),
      Form("Reco: %s", v.plotLabel),
      v.hReco, kRed + 1,
      v.xTitle, v.yTitle,
      Form("Reco:  %s", v.plotLabel));
    SaveCanvas(c, fout, pdfFile, false, false);
    delete c;
  }

  // =========================================================================
  // Pages 14-19: Overlay (MC + Reco)
  // =========================================================================
  for (int i = 0; i < nVars; ++i) {
    auto& v = vars[i];
    TCanvas* c = MakeOverlayCanvas(
      Form("cOvl_%s", v.tag),
      Form("MC + Reco: %s", v.plotLabel),
      v.hMC, v.hReco,
      v.xTitle, v.yTitle,
      v.plotLabel);
    SaveCanvas(c, fout, pdfFile, false, false);
    delete c;
  }

  // =========================================================================
  // Pages 20-25: Efficiency (closes PDF on last page)
  // =========================================================================
  for (int i = 0; i < nVars; ++i) {
    auto& v = vars[i];
    bool isLast = (i == nVars - 1);
    TCanvas* c = MakeEffCanvas(
      Form("cEff_%s", v.tag),
      Form("Efficiency: %s", v.plotLabel),
      v.hEff,
      v.xTitle,
      Form("Efficiency:  %s", v.plotLabel));
    SaveCanvas(c, fout, pdfFile, false, isLast);
    delete c;
  }

  // =========================================================================
  // Save all histograms to the output ROOT file
  // =========================================================================
  if (fout && fout->IsOpen()) {
    fout->cd();

    // Source histograms (12)
    if (hTrkPhiMC)
      hTrkPhiMC->Write("hTrkPhiMC");
    if (hTrkPhiReco)
      hTrkPhiReco->Write("hTrkPhiReco");
    if (hTrkEtaMC)
      hTrkEtaMC->Write("hTrkEtaMC");
    if (hTrkEtaReco)
      hTrkEtaReco->Write("hTrkEtaReco");
    if (hTrkPtMC)
      hTrkPtMC->Write("hTrkPtMC");
    if (hTrkPtReco)
      hTrkPtReco->Write("hTrkPtReco");
    if (hPairPhiMC)
      hPairPhiMC->Write("hPairPhiMC");
    if (hPairPhiReco)
      hPairPhiReco->Write("hPairPhiReco");
    if (hPairRapidityMC)
      hPairRapidityMC->Write("hPairRapidityMC");
    if (hPairRapidityReco)
      hPairRapidityReco->Write("hPairRapidityReco");
    if (hPairPtMC)
      hPairPtMC->Write("hPairPtMC");
    if (hPairPtReco)
      hPairPtReco->Write("hPairPtReco");

    // Efficiency histograms (6)
    if (hEffTrkPhi)
      hEffTrkPhi->Write();
    if (hEffTrkEta)
      hEffTrkEta->Write();
    if (hEffTrkPt)
      hEffTrkPt->Write();
    if (hEffPairPhi)
      hEffPairPhi->Write();
    if (hEffPairRap)
      hEffPairRap->Write();
    if (hEffPairPt)
      hEffPairPt->Write();

    // CutFlow histograms (2)
    if (hCutFlowMC)
      hCutFlowMC->Write("hCutFlowMC");
    if (hCutFlowReco)
      hCutFlowReco->Write("hCutFlowReco");

    fout->Close();
    ::Info("UPCMuonEfficiency", "Saved output to: %s", outFile);
  }

  fin->Close();
  ::Info("UPCMuonEfficiency", "PDF saved to: %s  (%d pages)", pdfFile, 1 + 3 * nVars * 2 + 1);
}
