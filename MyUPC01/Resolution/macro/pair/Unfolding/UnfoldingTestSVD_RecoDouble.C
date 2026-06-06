// ===========================================================
// UnfoldingTestSVD_RecoDouble.C
// SVD Unfolding for non-square response matrix
//   Truth (Gen) : N bins
//   Measured (Reco) : 2N bins
//
// Outputs:
//   Step3_2_SVD_RecoDouble_SelectedK.png  -- selected k values
//   Step3_2_SVD_RecoDouble_AllK_1toN.png  -- all k from 1 to nGenBins
//   Step3_2_SVD_RecoDouble_Dvector.png    -- d-vector diagnostic
// ===========================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include <TH1.h>
#include <TVirtualPad.h>

#include <Buttons.h>
#include <Rtypes.h>
#include <RtypesCore.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// RooUnfold headers
#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ROOUNFOLD_ROOT / include)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"

// ---------------------------------
// Helper Functions
// ---------------------------------
template <typename T>
T* LoadHistogram(TFile* file, const TString& histName)
{
  T* h = (T*)file->Get(histName);
  if (!h) {
    std::cerr << "[ERROR] Histogram not found: " << histName << std::endl;
    return nullptr;
  }
  T* hClone = (T*)h->Clone(histName + "_clone");
  hClone->SetDirectory(nullptr);
  return hClone;
}

// Extract bin edges from a histogram into a vector
std::vector<double> GetBinEdges_SVD(const TH1D* h)
{
  int n = h->GetNbinsX();
  std::vector<double> edges(n + 1);
  for (int i = 1; i <= n + 1; ++i)
    edges[i - 1] = h->GetXaxis()->GetBinLowEdge(i);
  return edges;
}

// ---------------------------------
// Draw unfolding result for a given k
// (top: Truth / Measured / Unfolded,  bottom: Unfolded/Truth ratio)
// Truth and Unfolded are both in Gen (coarse) binning.
// Measured is shown in Reco (fine) binning on the same plot.
// ---------------------------------
void DrawSVDResult_RD(TVirtualPad* parentPad, int k, TH1D* hGen, TH1D* hData, TH1D* hUnfolded)
{
  parentPad->cd();

  // ---- Top Pad ----
  TPad* padTop = new TPad(Form("padTopSVD_RD_k%d", k), "padTop", 0.0, 0.3, 1.0, 1.0);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.12);
  padTop->SetLogy();
  padTop->Draw();
  padTop->cd();

  TH1D* hGenDraw = (TH1D*)hGen->Clone(Form("hGenDrawSVD_RD_k%d", k));
  hGenDraw->SetTitle(Form("SVD Unfolding (k = %d)", k));
  hGenDraw->SetLineColor(kBlack);
  hGenDraw->SetLineWidth(2);
  hGenDraw->GetXaxis()->SetLabelSize(0);
  hGenDraw->GetXaxis()->SetTitleSize(0);
  hGenDraw->GetYaxis()->SetTitle("Counts");
  hGenDraw->GetYaxis()->SetTitleSize(0.04);
  hGenDraw->GetYaxis()->SetLabelSize(0.04);

  // Rebin hData (2N Reco bins) into Gen (N) binning for display
  int nGenB = hGen->GetNbinsX();
  std::vector<double> genEdgesDraw(nGenB + 1);
  for (int i = 1; i <= nGenB + 1; ++i)
    genEdgesDraw[i - 1] = hGen->GetXaxis()->GetBinLowEdge(i);
  TH1D* hDataDraw = (TH1D*)hData->Rebin(nGenB, Form("hDataDrawSVD_RD_k%d", k), genEdgesDraw.data());

  double maxY = std::max({hGenDraw->GetMaximum(), hDataDraw->GetMaximum(), hUnfolded->GetMaximum()});
  double minY = hGenDraw->GetMinimum(0);
  if (minY <= 0)
    minY = 0.1;

  hGenDraw->SetMaximum(maxY * 10.0);
  hGenDraw->SetMinimum(minY * 0.1);
  hGenDraw->Draw("HIST");

  hDataDraw->SetLineColor(kBlue + 1);
  hDataDraw->SetMarkerColor(kBlue + 1);
  hDataDraw->SetMarkerStyle(24);
  hDataDraw->Draw("PE SAME");

  hUnfolded->SetLineColor(kRed + 1);
  hUnfolded->SetMarkerColor(kRed + 1);
  hUnfolded->SetMarkerStyle(20);
  hUnfolded->Draw("PE SAME");

  TLegend* leg = new TLegend(0.65, 0.65, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hGenDraw, "Truth (Test Gen)", "l");
  leg->AddEntry(hDataDraw, "Measured (Test Data)", "pe");
  leg->AddEntry(hUnfolded, Form("Unfolded (k=%d)", k), "pe");
  leg->Draw();

  // ---- Bottom Pad (ratio: Unfolded / Truth, both in Gen binning) ----
  parentPad->cd();
  TPad* padBot = new TPad(Form("padBotSVD_RD_k%d", k), "padBot", 0.0, 0.0, 1.0, 0.3);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.12);
  padBot->SetRightMargin(0.05);
  padBot->Draw();
  padBot->cd();

  TH1D* hRatio = (TH1D*)hUnfolded->Clone(Form("hRatioSVD_RD_k%d", k));
  hRatio->Divide(hGen);
  hRatio->SetTitle("");
  hRatio->GetYaxis()->SetTitle("Unfolded / Truth");
  hRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetLabelSize(0.11);
  hRatio->GetYaxis()->SetTitleSize(0.11);
  hRatio->GetYaxis()->SetTitleOffset(0.5);
  hRatio->GetXaxis()->SetTitle("p_{T}^{2} (GeV^{2}/c^{2})");
  hRatio->GetXaxis()->SetLabelSize(0.12);
  hRatio->GetXaxis()->SetTitleSize(0.13);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->Draw("PE");

  TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kBlack);
  line->Draw("SAME");
}

// ---------------------------------
// Draw d-vector diagnostic plot
// ---------------------------------
void DrawDVector_RD(TCanvas* c, TH1D* hDvec)
{
  c->cd();
  c->SetLogy();
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.14);
  c->SetRightMargin(0.06);
  c->SetTopMargin(0.08);

  TH1D* hD = (TH1D*)hDvec->Clone("hDvec_draw_rd");
  hD->SetTitle("");
  hD->SetLineColor(kBlue + 1);
  hD->SetMarkerColor(kBlue + 1);
  hD->SetMarkerStyle(20);
  hD->SetMarkerSize(1.3);
  hD->GetXaxis()->SetTitle("Component index  i");
  hD->GetXaxis()->SetTitleSize(0.05);
  hD->GetXaxis()->SetLabelSize(0.04);
  hD->GetYaxis()->SetTitle("|d_{i}|");
  hD->GetYaxis()->SetTitleSize(0.05);
  hD->GetYaxis()->SetLabelSize(0.04);
  hD->GetYaxis()->SetTitleOffset(1.4);
  hD->Draw("PE");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.16, 0.92, "#bf{SVD d-vector}  (choose k where |d_{i}| drops to noise level)");
}

// =====================
// Main
// =====================
void UnfoldingTestSVD_RecoDouble()
{
  // ===========================
  // Settings
  // ===========================
  // MC (Training for Response Matrix) — built with RebinFlatStats_RecoDouble
  const TString mcFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/15697/Coherent/flatstats/5bin/Step2_Rebinned_RecoDouble.root";
  const TString histNameGen = "hGenPt2_rebin";    // nGenBins bins (truth axis)
  const TString histNameReco = "hRecoPt2_rebin";  // 2*nGenBins bins (measured axis)
  const TString histNameMatrix = "hResponseMatrix";

  // Test Data (Closure Test) — fine-binned MC from same production
  const TString dataFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/15697/Coherent/Step1_merged.root";
  const TString histNameDataFine = "hRecoPt2_Test";
  const TString histNameGenFine = "hGenPt2_Test";

  // Output
  const TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/15697/Coherent/flatstats/5bin/";

  // Selected k values (k <= nGenBins; adjust after checking the d-vector plot)
  const std::vector<int> kSelected = {1, 2, 3, 4};

  // ===========================
  // Initialization
  // ===========================
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  // --- Load Training MC ---
  TFile* fMC = TFile::Open(mcFile, "READ");
  if (!fMC || fMC->IsZombie()) {
    std::cerr << "[ERROR] Cannot open MC file: " << mcFile << std::endl;
    return;
  }
  TH1D* hGenTrain = LoadHistogram<TH1D>(fMC, histNameGen);
  TH1D* hRecoTrain = LoadHistogram<TH1D>(fMC, histNameReco);
  TH2D* hMat = LoadHistogram<TH2D>(fMC, histNameMatrix);
  fMC->Close();

  if (!hGenTrain || !hRecoTrain || !hMat) {
    std::cerr << "[ERROR] Missing MC histograms. Aborting." << std::endl;
    return;
  }

  int nGenBins = hGenTrain->GetNbinsX();
  int nRecoBins = hRecoTrain->GetNbinsX();

  std::cout << "[INFO] Gen (Truth) bins : " << nGenBins << std::endl;
  std::cout << "[INFO] Reco (Measured) bins: " << nRecoBins << std::endl;
  if (nRecoBins != 2 * nGenBins) {
    std::cerr << "[WARN] Expected nRecoBins == 2 * nGenBins, but got "
              << nRecoBins << " vs 2*" << nGenBins << "=" << 2 * nGenBins << std::endl;
  }

  // --- Load Test Data ---
  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[ERROR] Cannot open Data file: " << dataFile << std::endl;
    return;
  }
  TH1D* hDataFine = LoadHistogram<TH1D>(fData, histNameDataFine);
  TH1D* hGenTestFine = LoadHistogram<TH1D>(fData, histNameGenFine);
  fData->Close();

  if (!hDataFine || !hGenTestFine) {
    std::cerr << "[ERROR] Missing Data histograms. Aborting." << std::endl;
    return;
  }

  // ===========================
  // Rebin Test Data (Reco bins) and Test Truth (Gen bins) separately
  // ===========================
  std::vector<double> genEdges = GetBinEdges_SVD(hGenTrain);
  std::vector<double> recoEdges = GetBinEdges_SVD(hRecoTrain);

  // Test truth → Gen (coarse) binning
  TH1D* hGenTest = (TH1D*)hGenTestFine->Rebin(nGenBins, "hGenTest_Rebinned", genEdges.data());

  // Test data → Reco (fine) binning
  TH1D* hDataTest = (TH1D*)hDataFine->Rebin(nRecoBins, "hDataTest_Rebinned", recoEdges.data());

  std::cout << "[INFO] Test Truth rebinned into " << nGenBins << " Gen bins." << std::endl;
  std::cout << "[INFO] Test Data  rebinned into " << nRecoBins << " Reco bins." << std::endl;
  std::cout << "[INFO] Valid k range: [1, " << nGenBins << "]" << std::endl;

  // ===========================
  // RooUnfolding: Build Response (non-square: Reco=2N, Gen=N)
  // ===========================
  std::cout << "[INFO] Constructing RooUnfoldResponse (non-square)" << std::endl;
  RooUnfoldResponse response(hRecoTrain, hGenTrain, hMat);

  // ===========================
  // Canvas 1: Selected k values
  // ===========================
  int nPads = (int)kSelected.size();
  int cols1 = (nPads > 2) ? 2 : nPads;
  int rows1 = (nPads > 2) ? 2 : 1;

  TCanvas* c1 = new TCanvas("c1_rd", "SVD Unfolding (RecoDouble) - Selected k", cols1 * 600, rows1 * 500);
  c1->Divide(cols1, rows1);

  for (size_t i = 0; i < kSelected.size(); ++i) {
    int k = kSelected[i];
    if (k < 1 || k > nGenBins) {
      std::cerr << "[WARN] k=" << k << " out of valid range [1," << nGenBins << "], skipping." << std::endl;
      continue;
    }
    std::cout << "[INFO] SVD Unfolding k=" << k << std::endl;

    RooUnfoldSvd unfoldSvd(&response, hDataTest, k);
    TH1D* hUnfolded = (TH1D*)unfoldSvd.Hreco()->Clone(Form("hUnfoldSVD_RD_k%d", k));

    // Unfolded result is in Gen binning — compare directly with hGenTest
    DrawSVDResult_RD(c1->cd(i + 1), k, hGenTest, hDataTest, hUnfolded);
  }

  c1->SaveAs(outDir + "Step3_2_SVD_RecoDouble_SelectedK.png");
  std::cout << "[Done] Saved: Step3_2_SVD_RecoDouble_SelectedK.png" << std::endl;

  // ===========================
  // Canvas 2: All k from 1 to nGenBins
  // ===========================
  int colsAll = 2;
  int rowsAll = (nGenBins + colsAll - 1) / colsAll;

  TCanvas* c2 = new TCanvas("c2_rd", "SVD Unfolding (RecoDouble) - All k", colsAll * 600, rowsAll * 500);
  c2->Divide(colsAll, rowsAll);

  for (int k = 1; k <= nGenBins; ++k) {
    std::cout << "[INFO] SVD Unfolding k=" << k << " (all-k canvas)" << std::endl;

    RooUnfoldSvd unfoldSvd_all(&response, hDataTest, k);
    TH1D* hUnfolded_all = (TH1D*)unfoldSvd_all.Hreco()->Clone(Form("hUnfoldSVD_RD_all_k%d", k));

    DrawSVDResult_RD(c2->cd(k), k, hGenTest, hDataTest, hUnfolded_all);
  }

  c2->SaveAs(outDir + "Step3_2_SVD_RecoDouble_AllK_1toN.png");
  std::cout << "[Done] Saved: Step3_2_SVD_RecoDouble_AllK_1toN.png" << std::endl;

  // ===========================
  // Canvas 3: d-vector diagnostic
  // Run with k=nGenBins to access all singular value components.
  // ===========================
  RooUnfoldSvd unfoldSvdFull(&response, hDataTest, nGenBins);
  unfoldSvdFull.Hreco();

  TH1D* hDvec = (TH1D*)unfoldSvdFull.Impl()->GetD()->Clone("hDvec_rd");
  hDvec->SetDirectory(nullptr);

  TCanvas* c3 = new TCanvas("c3_rd", "SVD d-vector diagnostic (RecoDouble)", 800, 600);
  DrawDVector_RD(c3, hDvec);
  c3->SaveAs(outDir + "Step3_2_SVD_RecoDouble_Dvector.png");
  std::cout << "[Done] Saved: Step3_2_SVD_RecoDouble_Dvector.png" << std::endl;
}
