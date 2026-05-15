// ===========================================================
// UnfoldTestLcurve.C
// TUnfold with L-curve Optimization
// ==========================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TUnfoldDensity.h"
#include <TUnfold.h>

#include <algorithm>
#include <iostream>

// ヒストグラム読み込み用ヘルパー関数
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

void UnfoldTestLcurve()
{
  // ===========================
  // 1. Settings
  // ===========================
  const TString mcFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/TopDown/30/Step2_Rebinned.root";
  const TString histNameMatrix = "hResponseMatrix";

  const TString dataFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Coherent/Step1_merged.root";
  const TString histNameDataFine = "hRecoPt2_Test";
  const TString histNameGenFine = "hGenPt2_Test";

  gStyle->SetOptStat(0);

  // ===========================
  // 2. Load Data & MC
  // ===========================
  TFile* fMC = TFile::Open(mcFile, "READ");
  TH2D* hMat = LoadHistogram<TH2D>(fMC, histNameMatrix);
  fMC->Close();

  TFile* fData = TFile::Open(dataFile, "READ");
  TH1D* hDataFine = LoadHistogram<TH1D>(fData, histNameDataFine);
  TH1D* hGenTestFine = LoadHistogram<TH1D>(fData, histNameGenFine);
  fData->Close();

  if (!hMat || !hDataFine || !hGenTestFine)
    return;

  // ===========================
  // 3. Rebin Test Data based on Matrix
  // ===========================
  // ※ hMatのY軸(Reco)のビン数とエッジを取得すると仮定
  int nBinsReco = hMat->GetNbinsY();
  double* edgesReco = new double[nBinsReco + 1];
  for (int i = 1; i <= nBinsReco + 1; ++i) {
    edgesReco[i - 1] = hMat->GetYaxis()->GetBinLowEdge(i);
  }
  TH1D* hDataTest = (TH1D*)hDataFine->Rebin(nBinsReco, "hDataTest_Rebinned", edgesReco);

  // 真の分布(Gen)のリビニング -> 行列のX軸(5ビン)を使用
  int nBinsGen = hMat->GetNbinsX();
  double* edgesGen = new double[nBinsGen + 1];
  for (int i = 1; i <= nBinsGen + 1; ++i) {
    edgesGen[i - 1] = hMat->GetXaxis()->GetBinLowEdge(i);
  }
  TH1D* hGenTest = (TH1D*)hGenTestFine->Rebin(nBinsGen, "hGenTest_Rebinned", edgesGen);

  delete[] edgesReco;
  delete[] edgesGen;

  // ===========================
  // 4. Setup TUnfold
  // ===========================
  TUnfold::ERegMode regModeCur = TUnfold::kRegModeCurvature;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;

  TUnfoldDensity unfold(hMat, TUnfold::kHistMapOutputHoriz, regModeCur, constraintMode, densityFlags);

  if (unfold.SetInput(hDataTest) >= 10000) {
    std::cerr << "[WARNING] Unfolding result may be wrong. Check input bins." << std::endl;
  }

  // ===========================
  // 5. L-Curve Scan
  // ===========================
  Int_t nScan = 30; // スキャンするポイント数
  TSpline *logTauX, *logTauY;
  TGraph* lCurve;

  std::cout << "[INFO] Scanning L-curve..." << std::endl;
  Int_t iBest = unfold.ScanLcurve(nScan, 0.0, 0.0, &lCurve, &logTauX, &logTauY);

  std::cout << "Best tau = " << unfold.GetTau() << std::endl;
  std::cout << "chi**2 = " << unfold.GetChi2A() << " + " << unfold.GetChi2L() << " / " << unfold.GetNdf() << std::endl;

  // 最適点の座標取得
  Double_t t[1], x[1], y[1];
  logTauX->GetKnot(iBest, t[0], x[0]);
  logTauY->GetKnot(iBest, t[0], y[0]);
  TGraph* bestLcurve = new TGraph(1, x, y);
  bestLcurve->SetMarkerColor(kRed);
  bestLcurve->SetMarkerStyle(29);
  bestLcurve->SetMarkerSize(2.0);

  // ===========================
  // 6. Retrieve Results
  // ===========================
  TH1* hUnfolded = unfold.GetOutput("hUnfolded");

  // ===========================
  // 7. Plotting
  // ===========================
  TCanvas* c1 = new TCanvas("c1", "TUnfold L-curve Results", 1200, 800);
  c1->Divide(2, 2);

  // Pad 1: L-Curve
  c1->cd(1);
  lCurve->SetTitle("L-Curve;log(#chi^{2}_{L});log(#chi^{2}_{A})");
  lCurve->Draw("AL");
  bestLcurve->Draw("* SAME");

  // Pad 2: tau vs chi2
  c1->cd(2);
  logTauX->SetTitle("L-Curve X-coordinate vs log(#tau);log(#tau);log(#chi^{2}_{L})");
  logTauX->Draw();

  // Pad 3: Spectra Comparison
  c1->cd(3);
  gPad->SetLogy();
  hGenTest->SetLineColor(kBlack);
  hGenTest->SetLineWidth(2);
  hGenTest->SetTitle("Spectra Comparison;p_{T}^{2} (GeV^{2}/c^{2});Counts");
  hGenTest->Draw("HIST");

  hDataTest->SetLineColor(kBlue);
  hDataTest->SetMarkerStyle(24);
  hDataTest->SetMarkerColor(kBlue);
  hDataTest->Draw("PE SAME");

  hUnfolded->SetLineColor(kRed);
  hUnfolded->SetMarkerStyle(20);
  hUnfolded->SetMarkerColor(kRed);
  hUnfolded->Draw("PE SAME");

  TLegend* leg = new TLegend(0.5, 0.7, 0.88, 0.88);
  leg->AddEntry(hGenTest, "Truth (Gen)", "l");
  leg->AddEntry(hDataTest, "Measured (Reco)", "pe");
  leg->AddEntry(hUnfolded, "Unfolded", "pe");
  leg->Draw();

  // Pad 4: Ratio
  c1->cd(4);
  TH1D* hRatio = (TH1D*)hUnfolded->Clone("hRatio");
  hRatio->Divide(hGenTest);
  hRatio->SetTitle("Ratio (Unfolded / Truth);p_{T}^{2} (GeV^{2}/c^{2});Unfolded / Truth");
  hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hRatio->Draw("PE");
  TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw("SAME");

  c1->SaveAs("TUnfold_Lcurve_Result.png");
}
