// ===========================================================
// UnfoldTestIterations.C
// Bayesian Unfolding
// ==========================================================
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
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TVirtualPad.h>

#include <Buttons.h>
#include <Rtypes.h>
#include <RtypesCore.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

// RooUnfold headers
#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ROOUNFOLD_ROOT / include)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

// ---------------------------------
// Helper Functions
// ---------------------------------
template <typename T>
T* LoadHistogram(TFile* file, const TString& histName)
{
  T* h = (T*)file->Get(histName);
  if (!h) {
    std::cerr << "[ERROR] Histogram not found:" << histName << std::endl;
    return nullptr;
  }
  T* hClone = (T*)h->Clone(histName + "_clone");
  hClone->SetDirectory(nullptr);
  return hClone;
}

void DrawUnfoldedResult(TVirtualPad* parentPad, int iter, TH1D* hGen, TH1D* hData, TH1D* hUnfolded)
{
  parentPad->cd();

  // ======================================
  // Top Pad
  // =====================================
  TPad* padTop = new TPad(Form("padTop_%d", iter), "padTop", 0.0, 0.3, 1.0, 1.0);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.12);
  padTop->SetLogy();
  padTop->Draw();
  padTop->cd();

  // Gen (Test Truth)
  TH1D* hGenDraw = (TH1D*)hGen->Clone(Form("hGenDraw_%d", iter));
  hGenDraw->SetTitle(Form("Bayes Unfolding (Iter = %d)", iter));
  hGenDraw->SetLineColor(kBlack);
  hGenDraw->SetLineWidth(2);
  hGenDraw->GetXaxis()->SetLabelSize(0);
  hGenDraw->GetXaxis()->SetTitleSize(0);
  hGenDraw->GetYaxis()->SetTitle("Counts");
  hGenDraw->GetYaxis()->SetTitleSize(0.04);
  hGenDraw->GetYaxis()->SetLabelSize(0.04);

  double maxY = std::max({hGenDraw->GetMaximum(), hData->GetMaximum(), hUnfolded->GetMaximum()});
  double minY = hGenDraw->GetMinimum();
  if (minY <= 0)
    minY = 0.1;

  hGenDraw->SetMaximum(maxY * 10.0);
  hGenDraw->SetMinimum(minY * 0.1);
  hGenDraw->Draw("HIST");

  // Draw DATA (Test Reco)
  TH1D* hDataDraw = (TH1D*)hData->Clone(Form("hDataDraw_%d", iter));
  hDataDraw->SetLineColor(kBlue + 1);
  hDataDraw->SetMarkerColor(kBlue + 1);
  hDataDraw->SetMarkerStyle(24);
  hDataDraw->Draw("PE SAME");

  // Draw UNFOLD
  hUnfolded->SetLineColor(kRed + 1);
  hUnfolded->SetMarkerColor(kRed + 1);
  hUnfolded->SetMarkerStyle(20);
  hUnfolded->Draw("PE SAME");

  // Legend
  TLegend* leg = new TLegend(0.65, 0.65, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hGenDraw, "Truth (Test Gen)", "l");
  leg->AddEntry(hDataDraw, "Measured (Test Data)", "pe");
  leg->AddEntry(hUnfolded, Form("Unfolded (iter=%d)", iter), "pe");
  leg->Draw();

  // =========================================
  // Bottom Pad
  // ========================================
  parentPad->cd();
  TPad* padBot = new TPad(Form("padBot_%d", iter), "padBot", 0.0, 0.0, 1.0, 0.3);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.12);
  padBot->SetRightMargin(0.05);
  padBot->Draw();
  padBot->cd();

  TH1D* hRatio = (TH1D*)hUnfolded->Clone(Form("hRatio_%d", iter));
  hRatio->Divide(hGen); // Unfolded / Truth (Test Gen)
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

// =====================
// Main
// =====================
void UnfoldingTestIteration()
{
  // ===========================
  // Setting
  // ===========================
  // MC (Training for Response Matrix)
  const TString mcFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0415/Coherent/constantdiff/5bin/Step2_Rebinned.root";
  const TString histNameGen = "hGenPt2_rebin";      // Needed to build Response object
  const TString histNameReco = "hRecoPt2_rebin";    // Needed to build Response object
  const TString histNameMatrix = "hResponseMatrix"; // Needed to build Response object

  // Test Data (Closure Test)
  const TString dataFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0415/Coherent/Step1_merged.root";
  const TString histNameDataFine = "hRecoPt2_Test"; // To be unfolded
  const TString histNameGenFine = "hGenPt2_Test";   // True distribution for evaluation

  // Output
  const TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0415/Coherent/constantdiff/5bin/";
  const TString outFileName = "Step3_1_4Pads_Comparison.png";

  const std::vector<int> iters = {1, 3, 5, 7};

  // ==========================
  // Initialization
  // =========================
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  // --- Load Training MC ---
  TFile* fMC = TFile::Open(mcFile, "READ");
  if (!fMC || fMC->IsZombie()) {
    std::cerr << "[ERROR] Cannot open MC file:" << mcFile << std::endl;
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
    std::cerr << "[ERROR] Cannot open Data file:" << dataFile << std::endl;
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

  // Get bin edges from Training MC histogram
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

  std::cout << "[INFO] Successfully rebinned Test Data & Test Truth into " << nBins << " bins." << std::endl;

  // ===========================
  // RooUnfolding Processing
  // ===========================
  std::cout << "[INFO] Constructing RooUnfoldResponse with Training Data" << std::endl;
  RooUnfoldResponse response(hRecoTrain, hGenTrain, hMat);

  // Canvas setting
  int nPads = iters.size();
  int cols = (nPads > 2) ? 2 : nPads;
  int rows = (nPads > 2) ? 2 : 1;

  TCanvas* c1 = new TCanvas("c1", "Bayes Iteration Comparison", cols * 600, rows * 500);
  c1->Divide(cols, rows);

  for (size_t i = 0; i < iters.size(); ++i) {
    int iter = iters[i];
    std::cout << "[INFO] Running Bayes Unfolding with Iteration = " << iter << std::endl;

    // Unfolding test
    RooUnfoldBayes unfoldBayes(&response, hDataTest, iter);
    TH1D* hUnfolded = (TH1D*)unfoldBayes.Hreco()->Clone(Form("hUnfold_iter%d", iter));

    DrawUnfoldedResult(c1->cd(i + 1), iter, hGenTest, hDataTest, hUnfolded);
  }

  // =============================
  // Save Output
  // =============================
  c1->SaveAs(outDir + outFileName);
  std::cout << "[Done] Plot saved to:" << outDir + outFileName << std::endl;
}
