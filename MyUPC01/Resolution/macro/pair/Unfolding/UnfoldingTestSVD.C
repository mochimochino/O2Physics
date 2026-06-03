// ===========================================================
// UnfoldingTestSVD.C
// SVD Unfolding (Hoecker & Kartvelishvili method)
//
// Outputs:
//   Step3_2_SVD_SelectedK.png  -- selected k values (2x2 canvas)
//   Step3_2_SVD_AllK_1toN.png  -- all k from 1 to nBins
//   Step3_2_SVD_Dvector.png    -- d-vector diagnostic for k selection
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

// ---------------------------------
// Draw unfolding result for a given k
// (top: Truth / Measured / Unfolded,  bottom: Unfolded/Truth ratio)
// ---------------------------------
void DrawSVDResult(TVirtualPad* parentPad, int k, TH1D* hGen, TH1D* hData, TH1D* hUnfolded)
{
  parentPad->cd();

  // ---- Top Pad ----
  TPad* padTop = new TPad(Form("padTopSVD_k%d", k), "padTop", 0.0, 0.3, 1.0, 1.0);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.12);
  padTop->SetLogy();
  padTop->Draw();
  padTop->cd();

  TH1D* hGenDraw = (TH1D*)hGen->Clone(Form("hGenDrawSVD_k%d", k));
  hGenDraw->SetTitle(Form("SVD Unfolding (k = %d)", k));
  hGenDraw->SetLineColor(kBlack);
  hGenDraw->SetLineWidth(2);
  hGenDraw->GetXaxis()->SetLabelSize(0);
  hGenDraw->GetXaxis()->SetTitleSize(0);
  hGenDraw->GetYaxis()->SetTitle("Counts");
  hGenDraw->GetYaxis()->SetTitleSize(0.04);
  hGenDraw->GetYaxis()->SetLabelSize(0.04);

  double maxY = std::max({hGenDraw->GetMaximum(), hData->GetMaximum(), hUnfolded->GetMaximum()});
  double minY = hGenDraw->GetMinimum(0); // positive minimum for log scale
  if (minY <= 0)
    minY = 0.1;

  hGenDraw->SetMaximum(maxY * 10.0);
  hGenDraw->SetMinimum(minY * 0.1);
  hGenDraw->Draw("HIST");

  TH1D* hDataDraw = (TH1D*)hData->Clone(Form("hDataDrawSVD_k%d", k));
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

  // ---- Bottom Pad (ratio) ----
  parentPad->cd();
  TPad* padBot = new TPad(Form("padBotSVD_k%d", k), "padBot", 0.0, 0.0, 1.0, 0.3);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.12);
  padBot->SetRightMargin(0.05);
  padBot->Draw();
  padBot->cd();

  TH1D* hRatio = (TH1D*)hUnfolded->Clone(Form("hRatioSVD_k%d", k));
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
//
// The d-vector d_i = (U^T b)_i is the measured data
// projected onto the singular vector basis.
// |d_i| is large for small i (signal) and falls to
// the noise floor for large i.
// Choose k at the index where the drop-off begins.
// ---------------------------------
void DrawDVector(TCanvas* c, TH1D* hDvec)
{
  c->cd();
  c->SetLogy();
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.14);
  c->SetRightMargin(0.06);
  c->SetTopMargin(0.08);

  TH1D* hD = (TH1D*)hDvec->Clone("hDvec_draw");
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

  lat.SetTextSize(0.033);
  lat.SetTextColor(kGray + 2);
  lat.DrawLatex(0.16, 0.85, "Small k  #rightarrow  over-smoothing (bias)");
  lat.DrawLatex(0.16, 0.80, "Large k  #rightarrow  under-smoothing (variance)");
}

// =====================
// Main
// =====================
void UnfoldingTestSVD()
{
  // ===========================
  // Settings
  // ===========================
  // MC (Training for Response Matrix)
  const TString mcFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/16140/Coherent/flatstats/5bin/Step2_Rebinned.root";
  const TString histNameGen = "hGenPt2_rebin";
  const TString histNameReco = "hRecoPt2_rebin";
  const TString histNameMatrix = "hResponseMatrix";

  // Test Data (Closure Test)
  const TString dataFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/16140/Coherent/Step1_merged.root";
  const TString histNameDataFine = "hRecoPt2_Test";
  const TString histNameGenFine = "hGenPt2_Test";

  // Output
  const TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/ResolutionAndEfficiency/16140/Coherent/flatstats/5bin/";

  // Selected k values for c1 (adjust after checking nBins and the d-vector plot)
  const std::vector<int> kSelected = {2, 3, 4, 5};

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
  // Data Rebin (Applying MC Bins to Test Data)
  // ===========================
  int nBins = hRecoTrain->GetNbinsX();
  double* binEdges = new double[nBins + 1];

  if (hRecoTrain->GetXaxis()->GetXbins()->GetSize() > 0) {
    const double* arr = hRecoTrain->GetXaxis()->GetXbins()->GetArray();
    std::copy(arr, arr + nBins + 1, binEdges);
  } else {
    for (int i = 1; i <= nBins + 1; ++i) {
      binEdges[i - 1] = hRecoTrain->GetXaxis()->GetBinLowEdge(i);
    }
  }

  TH1D* hDataTest = (TH1D*)hDataFine->Rebin(nBins, "hDataTest_Rebinned", binEdges);
  TH1D* hGenTest = (TH1D*)hGenTestFine->Rebin(nBins, "hGenTest_Rebinned", binEdges);
  delete[] binEdges;

  std::cout << "[INFO] Rebinned Test Data & Test Truth into " << nBins << " bins." << std::endl;
  std::cout << "[INFO] Valid k range: [1, " << nBins << "]" << std::endl;

  // ===========================
  // RooUnfolding: Build Response
  // ===========================
  std::cout << "[INFO] Constructing RooUnfoldResponse with Training MC" << std::endl;
  RooUnfoldResponse response(hRecoTrain, hGenTrain, hMat);

  // ===========================
  // Canvas 1: Selected k values
  // ===========================
  int nPads = (int)kSelected.size();
  int cols1 = (nPads > 2) ? 2 : nPads;
  int rows1 = (nPads > 2) ? 2 : 1;

  TCanvas* c1 = new TCanvas("c1", "SVD Unfolding - Selected k", cols1 * 600, rows1 * 500);
  c1->Divide(cols1, rows1);

  for (size_t i = 0; i < kSelected.size(); ++i) {
    int k = kSelected[i];
    if (k < 1 || k > nBins) {
      std::cerr << "[WARN] k=" << k << " out of valid range [1," << nBins << "], skipping." << std::endl;
      continue;
    }
    std::cout << "[INFO] SVD Unfolding k=" << k << std::endl;

    RooUnfoldSvd unfoldSvd(&response, hDataTest, k);
    TH1D* hUnfolded = (TH1D*)unfoldSvd.Hreco()->Clone(Form("hUnfoldSVD_k%d", k));

    DrawSVDResult(c1->cd(i + 1), k, hGenTest, hDataTest, hUnfolded);
  }

  c1->SaveAs(outDir + "Step3_2_SVD_SelectedK.png");
  std::cout << "[Done] Saved: Step3_2_SVD_SelectedK.png" << std::endl;

  // ===========================
  // Canvas 2: All k from 1 to nBins
  // ===========================
  int colsAll = 2;
  int rowsAll = (nBins + colsAll - 1) / colsAll; // ceil(nBins / 2)

  TCanvas* c2 = new TCanvas("c2", "SVD Unfolding - All k (1 to N)", colsAll * 600, rowsAll * 500);
  c2->Divide(colsAll, rowsAll);

  for (int k = 1; k <= nBins; ++k) {
    std::cout << "[INFO] SVD Unfolding k=" << k << " (all-k canvas)" << std::endl;

    RooUnfoldSvd unfoldSvd_all(&response, hDataTest, k);
    TH1D* hUnfolded_all = (TH1D*)unfoldSvd_all.Hreco()->Clone(Form("hUnfoldSVD_all_k%d", k));

    DrawSVDResult(c2->cd(k), k, hGenTest, hDataTest, hUnfolded_all);
  }

  c2->SaveAs(outDir + "Step3_2_SVD_AllK_1toN.png");
  std::cout << "[Done] Saved: Step3_2_SVD_AllK_1toN.png" << std::endl;

  // ===========================
  // Canvas 3: d-vector diagnostic
  //
  // Run SVD with k=nBins to access all singular value components.
  // The d-vector |d_i| = |(U^T b)_i| represents the measured data
  // in the singular vector basis.  Large values at small i indicate
  // signal; where |d_i| flattens to the noise floor, set k.
  //
  // API note: Impl() returns the underlying TSVDUnfold* object.
  // If compilation fails, check your RooUnfold version -- some versions
  // expose this via Svd() instead of Impl().
  // ===========================
  RooUnfoldSvd unfoldSvdFull(&response, hDataTest, nBins);
  unfoldSvdFull.Hreco(); // trigger internal SVD decomposition

  TH1D* hDvec = (TH1D*)unfoldSvdFull.Impl()->GetD()->Clone("hDvec");
  hDvec->SetDirectory(nullptr);

  TCanvas* c3 = new TCanvas("c3", "SVD d-vector diagnostic", 800, 600);
  DrawDVector(c3, hDvec);
  c3->SaveAs(outDir + "Step3_2_SVD_Dvector.png");
  std::cout << "[Done] Saved: Step3_2_SVD_Dvector.png" << std::endl;
}
