// ===========================================================
// UnfoldTest.C
// TUnfold: Separated L-curve and Tau Scan Processes
// for Non-Square Matrix with Both Spectra Comparison
// Takuma Matsumoto
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

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

void UnfoldTest()
{
  // ===========================
  // Settings
  // ===========================
  const TString mcFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Incoherent/TopDown/40/Step2_Rebinned.root";
  const TString histNameMatrix = "hResponseMatrix";

  const TString dataFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/Incoherent/Step1_merged.root";
  const TString histNameDataFine = "hRecoPt2_Test";
  const TString histNameGenFine = "hGenPt2_Test";

  gStyle->SetOptStat(0);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetTitleOffset(1.6, "Y");
  gStyle->SetTitleOffset(1.2, "X");

  // ===========================
  // Load Data & MC
  // ===========================
  TFile* fMC = TFile::Open(mcFile, "READ");
  if (!fMC || fMC->IsZombie())
    return;
  TH2D* hMat = LoadHistogram<TH2D>(fMC, histNameMatrix);
  fMC->Close();

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie())
    return;
  TH1D* hDataFine = LoadHistogram<TH1D>(fData, histNameDataFine);
  TH1D* hGenTestFine = LoadHistogram<TH1D>(fData, histNameGenFine);
  fData->Close();

  if (!hMat || !hDataFine || !hGenTestFine)
    return;

  // ===========================
  // Extract Bin Edges from Matrix
  // ===========================
  int nBinsReco = hMat->GetNbinsY();
  std::vector<double> edgesReco(nBinsReco + 1);
  for (int i = 1; i <= nBinsReco; ++i) {
    edgesReco[i - 1] = hMat->GetYaxis()->GetBinLowEdge(i);
  }
  edgesReco[nBinsReco] = hMat->GetYaxis()->GetBinUpEdge(nBinsReco);

  int nBinsGen = hMat->GetNbinsX();
  std::vector<double> edgesGen(nBinsGen + 1);
  for (int i = 1; i <= nBinsGen; ++i) {
    edgesGen[i - 1] = hMat->GetXaxis()->GetBinLowEdge(i);
  }
  edgesGen[nBinsGen] = hMat->GetXaxis()->GetBinUpEdge(nBinsGen);

  // ===========================
  // Rebin Data
  // ===========================
  TH1D* hDataTest_RecoBins = (TH1D*)hDataFine->Rebin(nBinsReco, "hDataTest_RecoBins", edgesReco.data());
  TH1D* hDataTest_GenBins = (TH1D*)hDataFine->Rebin(nBinsGen, "hDataTest_GenBins", edgesGen.data());
  TH1D* hGenTest_GenBins = (TH1D*)hGenTestFine->Rebin(nBinsGen, "hGenTest_GenBins", edgesGen.data());

  // ===========================
  // Setup loop for 3 Regularization Modes
  // ===========================
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;

  const int nModes = 3;
  TUnfold::ERegMode modes[nModes] = {TUnfold::kRegModeSize, TUnfold::kRegModeDerivative, TUnfold::kRegModeCurvature};
  TString modeNames[nModes] = {"Size", "Derivative", "Curvature"};

  Int_t nScan = 40;

  for (int i = 0; i < nModes; ++i) {
    std::cout << "\n[INFO] =========================================" << std::endl;
    std::cout << "[INFO] Processing mode: " << modeNames[i] << std::endl;

    // 1. TUnfoldの初期化
    TUnfoldDensity unfold(hMat, TUnfold::kHistMapOutputHoriz, modes[i], constraintMode, densityFlags);
    if (unfold.SetInput(hDataTest_RecoBins) >= 10000) {
      std::cerr << "[WARNING] High condition number detected for mode " << modeNames[i] << std::endl;
    }

    // ===========================
    // 2. ScanLcurve
    // ===========================
    TSpline *logTauX = 0, *logTauY = 0;
    TGraph* lCurve = 0;
    unfold.ScanLcurve(nScan, 0.0, 0.0, &lCurve, &logTauX, &logTauY);
    Double_t tauLcurve = unfold.GetTau();
    TH1D* hUnfolded_L = (TH1D*)unfold.GetOutput(Form("hUnfolded_L_%s", modeNames[i].Data()));

    Double_t logTauBestL = std::log10(tauLcurve);
    Double_t xL_val = logTauX->Eval(logTauBestL);
    Double_t yL_val = logTauY->Eval(logTauBestL);
    TGraph* bestLcurve = new TGraph(1, &xL_val, &yL_val);
    bestLcurve->SetMarkerColor(kRed);
    bestLcurve->SetMarkerStyle(29);
    bestLcurve->SetMarkerSize(2.0);

    // ===========================
    // 3. ScanTau (RhoAvg)
    // ===========================
    TSpline* scanResult = 0;
    unfold.ScanTau(nScan, 1e-9, 1.0, &scanResult, TUnfoldDensity::kEScanTauRhoAvg, 0, 0, nullptr, nullptr, nullptr);
    Double_t tauRho = unfold.GetTau();
    TH1D* hUnfolded_Rho = (TH1D*)unfold.GetOutput(Form("hUnfolded_Rho_%s", modeNames[i].Data()));

    Double_t logTauBestRho = std::log10(tauRho);
    Double_t minRho = scanResult->Eval(logTauBestRho);
    TGraph* bestRho = new TGraph(1, &logTauBestRho, &minRho);
    bestRho->SetMarkerColor(kRed);
    bestRho->SetMarkerStyle(29);
    bestRho->SetMarkerSize(2.0);

    // ===========================
    // 4. Get Diagnostics (based on tauRho state)
    // ===========================
    TH2* hCorr = unfold.GetRhoIJtotal(Form("hCorr_%s", modeNames[i].Data()));
    TH1* hFolded = unfold.GetFoldedOutput(Form("hFolded_%s", modeNames[i].Data()));

    TH1D* hRhoI = new TH1D(Form("hRhoI_%s", modeNames[i].Data()), "Global Correlation #rho_{i};p_{T}^{2} (GeV^{2}/c^{2});#rho_{i}", nBinsGen, edgesGen.data());
    unfold.GetRhoI(hRhoI);

    std::cout << "  Best tau (L-curve) = " << tauLcurve << std::endl;
    std::cout << "  Best tau (RhoAvg)  = " << tauRho << std::endl;

    // ===========================
    // 5. Plotting (8 Pads per method)
    // ===========================
    TCanvas* c1 = new TCanvas(Form("c_%s", modeNames[i].Data()), Form("TUnfold Diagnostics - %s", modeNames[i].Data()), 2000, 1000);
    c1->Divide(4, 2);

    // --- Pad 1: L-Curve ---
    c1->cd(1);
    if (lCurve) {
      lCurve->SetTitle(Form("L-Curve (%s);log(#chi^{2}_{L});log(#chi^{2}_{A})", modeNames[i].Data()));
      lCurve->Draw("AL");
      bestLcurve->Draw("* SAME");
    }

    // --- Pad 2: tau vs RhoAvg ---
    c1->cd(2);
    if (scanResult) {
      scanResult->SetTitle(Form("Avg Global Correlation (%s);log(#tau);#bar{#rho}", modeNames[i].Data()));
      scanResult->SetLineColor(kBlue + 1);
      scanResult->SetLineWidth(2);
      scanResult->Draw();
      bestRho->Draw("* SAME");
    }

    // --- Pad 3: Spectra Comparison (Counts) ---
    c1->cd(3);
    // gPad->SetLogy();
    hGenTest_GenBins->SetLineColor(kBlack);
    hGenTest_GenBins->SetLineWidth(2);
    hGenTest_GenBins->SetTitle(Form("Spectra Comparison (%s);p_{T}^{2} (GeV^{2}/c^{2});Counts", modeNames[i].Data()));
    double yMax = hDataTest_GenBins->GetMaximum();
    double yMin = hDataTest_GenBins->GetMinimum();
    hGenTest_GenBins->GetYaxis()->SetRangeUser(yMin * 0.5, yMax * 1.5);
    // hGenTest_GenBins->GetYaxis()->SetMoreLogLabels(kTRUE);
    hGenTest_GenBins->Draw("HIST");

    hDataTest_GenBins->SetLineColor(kBlue);
    // hDataTest_GenBins->SetMarkerStyle(24);
    // hDataTest_GenBins->SetMarkerColor(kBlue);
    hDataTest_GenBins->SetLineWidth(2);
    hDataTest_GenBins->SetMarkerSize(0);
    hDataTest_GenBins->Draw("E SAME");

    hUnfolded_Rho->SetLineColor(kRed);
    // hUnfolded_Rho->SetMarkerStyle(24);
    // hUnfolded_Rho->SetMarkerColor(kRed);
    hUnfolded_Rho->SetLineWidth(2);
    hUnfolded_Rho->SetMarkerSize(0);
    hUnfolded_Rho->Draw("PE SAME");

    hUnfolded_L->SetLineColor(kGreen + 2);
    // hUnfolded_L->SetMarkerStyle(24);
    // hUnfolded_L->SetMarkerColor(kGreen + 2);
    hUnfolded_L->SetLineWidth(2);
    hUnfolded_L->SetMarkerSize(0);
    hUnfolded_L->Draw("PE SAME");

    TLegend* leg = new TLegend(0.58, 0.70, 0.90, 0.90);
    leg->AddEntry(hGenTest_GenBins, "Truth (Gen)", "l");
    leg->AddEntry(hDataTest_GenBins, "Measured", "l");
    leg->AddEntry(hUnfolded_Rho, "Unfolded (#tau_{#rho})", "l");
    leg->AddEntry(hUnfolded_L, "Unfolded (#tau_{L})", "l");
    leg->Draw();

    // --- Pad 4: Ratio (Unfolded / Gen) ---
    c1->cd(4);
    TH1D* hRatio_Rho = (TH1D*)hUnfolded_Rho->Clone(Form("hRatio_Rho_%s", modeNames[i].Data()));
    hRatio_Rho->Divide(hGenTest_GenBins);
    hRatio_Rho->SetTitle(Form("Ratio (Unfolded / Truth) - %s;p_{T}^{2} (GeV^{2}/c^{2});Unfolded / Truth", modeNames[i].Data()));
    hRatio_Rho->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio_Rho->Draw("PE");

    TH1D* hRatio_L = (TH1D*)hUnfolded_L->Clone(Form("hRatio_L_%s", modeNames[i].Data()));
    hRatio_L->Divide(hGenTest_GenBins);
    hRatio_L->Draw("PE SAME");

    TLine* lineRatio = new TLine(hRatio_Rho->GetXaxis()->GetXmin(), 1.0, hRatio_Rho->GetXaxis()->GetXmax(), 1.0);
    lineRatio->SetLineStyle(2);
    lineRatio->Draw("SAME");

    // --- Pad 5: Correlation Matrix ---
    c1->cd(5);
    gPad->SetRightMargin(0.15);
    hCorr->SetTitle(Form("Correlation Matrix (%s, #tau_{#rho});Bin i;Bin j", modeNames[i].Data()));
    hCorr->Draw("COLZ");

    // --- Pad 6: Data Consistency Check (Folded / Reco) ---
    c1->cd(6);
    TH1D* hRatioReco = (TH1D*)hFolded->Clone(Form("hRatioReco_%s", modeNames[i].Data()));
    hRatioReco->Divide(hDataTest_RecoBins);
    hRatioReco->SetTitle(Form("Data Consistency (Folded / Reco) - %s;p_{T}^{2} (GeV^{2}/c^{2});Folded / Measured", modeNames[i].Data()));
    hRatioReco->SetMarkerColor(kMagenta);
    hRatioReco->SetLineColor(kMagenta);
    hRatioReco->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatioReco->Draw("PE");
    TLine* lineReco = new TLine(hRatioReco->GetXaxis()->GetXmin(), 1.0, hRatioReco->GetXaxis()->GetXmax(), 1.0);
    lineReco->SetLineStyle(2);
    lineReco->Draw("SAME");

    // --- Pad 7: Response (Migration) Matrix ---
    c1->cd(7);
    gPad->SetRightMargin(0.15);

    gStyle->SetPaintTextFormat(".1f");

    TH2D* hMigration = (TH2D*)hMat->Clone(Form("hMigration_%s", modeNames[i].Data()));
    hMigration->SetTitle("Response (Migration) Matrix;Gen p_{T}^{2} (GeV^{2}/c^{2});Reco p_{T}^{2} (GeV^{2}/c^{2})");

    hMigration->Scale(0.01);

    hMigration->SetMinimum(1e-7);

    hMigration->SetMarkerSize(1.5);
    hMigration->SetMarkerColor(kBlack);

    hMigration->Draw("COLZ TEXT");

    /*
    // --- Pad 7: Spectra Comparison (Density) ---
    c1->cd(7);
    gPad->SetLogy();
    TH1D* hGenTest_Density = (TH1D*)hGenTest_GenBins->Clone(Form("hGenTest_Density_%s", modeNames[i].Data()));
    hGenTest_Density->Scale(1.0, "width");
    hGenTest_Density->SetTitle(Form("Spectra Comparison (Density) - %s;p_{T}^{2} (GeV^{2}/c^{2});dN / dp_{T}^{2}", modeNames[i].Data()));

    TH1D* hDataTest_Density = (TH1D*)hDataTest_GenBins->Clone(Form("hDataTest_Density_%s", modeNames[i].Data()));
    hDataTest_Density->Scale(1.0, "width");

    TH1D* hUnfolded_Rho_Density = (TH1D*)hUnfolded_Rho->Clone(Form("hUnfolded_Rho_Density_%s", modeNames[i].Data()));
    hUnfolded_Rho_Density->Scale(1.0, "width");

    TH1D* hUnfolded_L_Density = (TH1D*)hUnfolded_L->Clone(Form("hUnfolded_L_Density_%s", modeNames[i].Data()));
    hUnfolded_L_Density->Scale(1.0, "width");

    hGenTest_Density->Draw("HIST");
    hDataTest_Density->Draw("PE SAME");
    hUnfolded_Rho_Density->Draw("PE SAME");
    hUnfolded_L_Density->Draw("PE SAME");
    */

    // --- Pad 8: Global Correlation per Bin ---
    c1->cd(8);
    hRhoI->SetTitle(Form("Global Correlation per Bin (%s);p_{T}^{2} (GeV^{2}/c^{2});#rho_{i}", modeNames[i].Data()));
    hRhoI->SetLineColor(kMagenta + 2);
    hRhoI->SetLineWidth(2);
    hRhoI->GetYaxis()->SetRangeUser(0.0, 1.05);
    hRhoI->Draw("HIST");

    c1->Update();
    c1->SaveAs(Form("TUnfold_Diagnostics_Separated_%s.png", modeNames[i].Data()));

    delete c1;
  }
}
