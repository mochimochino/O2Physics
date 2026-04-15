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

  // Gen
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

  // Draw DATA
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
  leg->AddEntry(hGenDraw, "Truth (Gen)", "l");
  leg->AddEntry(hDataDraw, "Measured (Data)", "pe");
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
  hRatio->Divide(hGen); // Unfolded / Truth
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
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0414/Incoherent/constantdiff/15bin/";
  const TString inFile = inDir + "Step2_Rebinned.root";
  const TString outDir = inDir;
  const TString outFileName = "Step3_1_4Pads_Comparison.png";

  const TString histNameGen = "hGenPt2_rebin";
  const TString histNameReco = "hRecoPt2_rebin";
  const TString histNameMatrix = "hResponseMatrix";
  const TString histNameData = "hRecoPt2_rebin"; // Need to change to real data

  const std::vector<int> iters = {1, 3, 5, 7};

  // ==========================
  // Initialization
  // =========================
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file:" << inFile << std::endl;
    return;
  }

  // Read file with helper
  TH1D* hGen = LoadHistogram<TH1D>(fIn, histNameGen);
  TH1D* hReco = LoadHistogram<TH1D>(fIn, histNameReco);
  TH2D* hMat = LoadHistogram<TH2D>(fIn, histNameMatrix);
  TH1D* hData = LoadHistogram<TH1D>(fIn, histNameData);
  fIn->Close();

  if (!hGen || !hReco || !hMat || !hData) {
    std::cerr << "[ERROR] Missing histograms" << std::endl;
    return;
  }

  // ===========================
  // RooUnfolding Processing
  // ===========================
  std::cout << "[Info] Constructing RooUnfoldResponse" << std::endl;
  RooUnfoldResponse response(hReco, hGen, hMat);

  // Canvas setting
  int nPads = iters.size();
  int cols = (nPads > 2) ? 2 : nPads;
  int rows = (nPads > 2) ? 2 : 1;

  TCanvas* c1 = new TCanvas("c1", "Bayes Iteration Comparison", cols * 600, rows * 500);
  c1->Divide(cols, rows);

  for (Size_t i = 0; i < iters.size(); ++i) {
    int iter = iters[i];
    std::cout << "[INFO] Runnging Vayes Unfolding with Iteration =" << iter << std::endl;

    RooUnfoldBayes unfoldBayes(&response, hData, iter);
    TH1D* hUnfolded = (TH1D*)unfoldBayes.Hreco()->Clone(Form("hUnfold_iter%d", iter));

    DrawUnfoldedResult(c1->cd(i + 1), iter, hGen, hData, hUnfolded);
  }

  // =============================
  // Save Output
  // =============================
  c1->SaveAs(outDir + outFileName);
  std::cout << "[Done] Plot saved to:" << outDir + outFileName << std::endl;
}
