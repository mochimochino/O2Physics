// ============================================================
//  Make2DHistFromBinnedHists.C
//  Helper task to make 2D histogram from binned 1D histograms
// ============================================================

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>
#include <vector>

// ----------------------------------------------------------
//  Get Bin edge from 1D histogram
// ----------------------------------------------------------
std::vector<double> GetBinEdges_Make2D(const TH1* h)
{
  std::vector<double> edges;
  if (!h)
    return edges;
  int n = h->GetNbinsX();
  edges.reserve(n + 1);
  for (int i = 1; i <= n + 1; ++i)
    edges.push_back(h->GetXaxis()->GetBinLowEdge(i));
  return edges;
}

// ----------------------------------------------------------
//  Main
// ----------------------------------------------------------
void Make2DHistFromBinnedHists(
  TString matrixFile = "", // User Settings
  TString matrixName = "hResponseMatrixPairPt2",
  TString binnedFile = "", // matrixFile
  TString xHistName = "hGenPt2",
  TString yHistName = "hRecoPt2",
  TString outFile = "Make2DHistFromBinnedHists_output.root")
{
  // ==========================================================
  // User Settings
  // ==========================================================
  if (matrixFile.IsNull() || matrixFile.IsWhitespace()) {
    matrixFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/CoherenrtJpsi/Step1_merged.root";
    matrixName = "hResponseMatrixPairPt2";
    binnedFile = ""; // matrixFile
    xHistName = "hGenPt2";
    yHistName = "hRecoPt2";
    outFile = "Make2DHistFromBinnedHists_output.root";
  }

  if (binnedFile.IsNull() || binnedFile.IsWhitespace())
    binnedFile = matrixFile;

  // ==========================================================
  // Read file and hist
  // ==========================================================
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  TFile* fMatrix = TFile::Open(matrixFile, "READ");
  if (!fMatrix || fMatrix->IsZombie()) {
    std::cerr << "[Error] Cannot open matrix file: " << matrixFile << std::endl;
    return;
  }
  TH2* hFineRaw = (TH2*)fMatrix->Get(matrixName);
  if (!hFineRaw) {
    std::cerr << "[Error] Histogram not found: " << matrixName << std::endl;
    fMatrix->Close();
    return;
  }
  TH2D* hFineMatrix = (TH2D*)hFineRaw->Clone("hFineMatrix_clone");
  hFineMatrix->SetDirectory(nullptr);
  fMatrix->Close();

  TFile* fBinned = TFile::Open(binnedFile, "READ");
  if (!fBinned || fBinned->IsZombie()) {
    std::cerr << "[Error] Cannot open binned file: " << binnedFile << std::endl;
    return;
  }
  TH1D* hGenX = (TH1D*)fBinned->Get(xHistName);
  TH1D* hGenY = (TH1D*)fBinned->Get(yHistName);
  if (!hGenX) {
    std::cerr << "[Error] X-axis histogram not found: " << xHistName << std::endl;
    fBinned->Close();
    return;
  }
  if (!hGenY) {
    std::cerr << "[Error] Y-axis histogram not found: " << yHistName << std::endl;
    fBinned->Close();
    return;
  }
  hGenX->SetDirectory(nullptr);
  hGenY->SetDirectory(nullptr);
  fBinned->Close();

  // ==========================================================
  // Get bin edge
  // ==========================================================
  std::vector<double> xEdges = GetBinEdges_Make2D(hGenX);
  std::vector<double> yEdges = GetBinEdges_Make2D(hGenY);
  const int nX = (int)xEdges.size() - 1;
  const int nY = (int)yEdges.size() - 1;

  if (nX <= 0 || nY <= 0) {
    std::cerr << "[Error] Invalid bin edges (nX=" << nX << ", nY=" << nY << ")." << std::endl;
    return;
  }

  std::cout << "\n[Info] Matrix file  : " << matrixFile << std::endl;
  std::cout << "[Info] Binned file  : " << binnedFile << std::endl;
  std::cout << "[Info] X-axis hist  : " << xHistName << "  (" << nX
            << " bins, " << xEdges.front() << " ~ " << xEdges.back() << ")" << std::endl;
  std::cout << "[Info] Y-axis hist  : " << yHistName << "  (" << nY
            << " bins, " << yEdges.front() << " ~ " << yEdges.back() << ")" << std::endl;

  // ==========================================================
  // 2D
  // ==========================================================
  TH2D* hMatrixRebinned = new TH2D(
    "hMatrixRebinned",
    "Response Matrix (Rebinned);Gen p_{T}^{2};Reco p_{T}^{2}",
    nX, xEdges.data(),
    nY, yEdges.data());

  for (int ix = 1; ix <= hFineMatrix->GetNbinsX(); ++ix) {
    double xCenter = hFineMatrix->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= hFineMatrix->GetNbinsY(); ++iy) {
      double yCenter = hFineMatrix->GetYaxis()->GetBinCenter(iy);
      double w = hFineMatrix->GetBinContent(ix, iy);
      if (w != 0.0)
        hMatrixRebinned->Fill(xCenter, yCenter, w);
    }
  }

  // ==========================================================
  // Purity & Stability
  // ==========================================================
  bool isSquare = (nX == nY);
  TH1D* hPurity = new TH1D("hPurity", "Purity;Gen bin;Purity", nX, xEdges.data());
  TH1D* hStability = new TH1D("hStability", "Stability;Reco bin;Stability", nY, yEdges.data());

  if (isSquare) {
    for (int i = 1; i <= nX; ++i) {
      double diag = hMatrixRebinned->GetBinContent(i, i);
      // double sumTrue = hMatrixRebinned->Integral(i, i, 1, nY); // Gen ビン i の行和
      // double sumReco = hMatrixRebinned->Integral(1, nX, i, i); // Reco ビン i の列和
      double sumTrue = hGenX->GetBinContent(i);
      double sumReco = hGenY->GetBinContent(i);
      hPurity->SetBinContent(i, (sumReco > 0) ? diag / sumReco : 0.0);
      hStability->SetBinContent(i, (sumTrue > 0) ? diag / sumTrue : 0.0);
    }
  }

  // ==========================================================
  // Cell / Gen [%] 2D ヒスト（全セル）
  //   Gen_ix     = hGenX (リビン済み Gen ヒスト)  の第 ix ビンの内容（X軸）
  //   Cell_ix_iy = 2D行列 (hMatrixRebinned) の第 ix, iy セルの内容
  // ==========================================================
  TH2D* hCellGenRatio = new TH2D(
    "hCellGenRatio",
    "(Reco + Gen) / Gen [%];Gen p_{T}^{2};Reco p_{T}^{2}",
    nX, xEdges.data(),
    nY, yEdges.data());

  for (int ix = 1; ix <= nX; ++ix) {
    double genVal = hGenX->GetBinContent(ix);

    double columnSum = 0.0;
    for (int iy = 1; iy <= nY; ++iy) {
      double cellVal = hMatrixRebinned->GetBinContent(ix, iy);
      columnSum += cellVal;
      double ratio = (genVal > 0) ? (cellVal / genVal) * 100.0 : 0.0;
      hCellGenRatio->SetBinContent(ix, iy, ratio);
    }
    // For check
    int fineBinMin = hFineMatrix->GetXaxis()->FindBin(xEdges[ix - 1] + 1e-6);
    int fineBinMax = hFineMatrix->GetXaxis()->FindBin(xEdges[ix] - 1e-6);

    double totalRecoIncludingOutliers = hFineMatrix->Integral(fineBinMin, fineBinMax, 0, hFineMatrix->GetNbinsY() + 1);

    double fullRatio = (genVal > 0) ? (totalRecoIncludingOutliers / genVal) * 100.0 : 0.0;

    std::cout << Form("Gen Bin %d: In-Range Sum = %.1f%% | Full Sum (incl. Overflow) = %.1f%%",
                      ix, (columnSum / genVal) * 100.0, fullRatio)
              << std::endl;
  }

  // ==========================================================
  // Output result on terminal
  // ==========================================================
  std::cout << "\n============================================================" << std::endl;
  std::cout << "  Rebinning Result Table" << std::endl;
  std::cout << "============================================================" << std::endl;

  if (isSquare) {
    std::cout << Form("%-5s | %-14s | %-14s | %-8s | %-8s | %-15s",
                      "Bin", "X Range (Gen)", "Y Range (Reco)",
                      "Purity%", "Stab%", "Diag/Gen (%)")
              << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    for (int i = 1; i <= nX; ++i) {
      double pur = hPurity->GetBinContent(i) * 100.0;
      double stab = hStability->GetBinContent(i) * 100.0;
      double ratio = hCellGenRatio->GetBinContent(i, i);
      std::cout << Form("%-5d | %6.4f - %6.4f | %6.4f - %6.4f | %-8.1f | %-8.1f | %-15.1f",
                        i,
                        xEdges[i - 1], xEdges[i],
                        yEdges[i - 1], yEdges[i],
                        pur, stab, ratio)
                << std::endl;
    }
  } else {
    std::cout << "Non-square matrix: printing Cell/Gen [%] for all cells." << std::endl;
    std::cout << Form("%-5s | %-5s | %-14s | %-14s | %-15s",
                      "iX", "iY", "X Range (Gen)", "Y Range (Reco)", "Cell/Gen (%)")
              << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    for (int ix = 1; ix <= nX; ++ix) {
      for (int iy = 1; iy <= nY; ++iy) {
        double ratio = hCellGenRatio->GetBinContent(ix, iy);
        std::cout << Form("%-5d | %-5d | %6.4f - %6.4f | %6.4f - %6.4f | %-15.1f",
                          ix, iy,
                          xEdges[ix - 1], xEdges[ix],
                          yEdges[iy - 1], yEdges[iy],
                          ratio)
                  << std::endl;
      }
    }
  }
  std::cout << "============================================================\n"
            << std::endl;

  // ==========================================================
  // Draw
  // ==========================================================
  gStyle->SetPaintTextFormat(".1f");
  TString outBase = outFile;
  outBase.ReplaceAll(".root", "");

  TCanvas* cMatrix = new TCanvas("cMatrix_2d", "Response Matrix (Rebinned)", 800, 700);
  cMatrix->SetRightMargin(0.15);
  hMatrixRebinned->Draw("COLZ");
  cMatrix->SetLogz();
  cMatrix->SaveAs(outBase + "_Matrix.png");

  TCanvas* cRatio = new TCanvas("cRatio_2d", "(Reco + Gen) / Gen [%]", 800, 700);
  cRatio->SetRightMargin(0.15);
  cRatio->SetLeftMargin(0.15);
  hCellGenRatio->Draw("COLZ TEXT");
  cRatio->SaveAs(outBase + "_CellGenRatio.png");

  if (isSquare) {
    TCanvas* cPS = new TCanvas("cPS_2d", "Purity & Stability", 900, 450);
    cPS->Divide(2, 1);
    cPS->cd(1);
    hPurity->GetYaxis()->SetRangeUser(0, 1);
    hPurity->Draw("HIST");
    cPS->cd(2);
    hStability->GetYaxis()->SetRangeUser(0, 1);
    hStability->Draw("HIST");
    cPS->SaveAs(outBase + "_PurityStability.png");
  }

  // ==========================================================
  // Save ROOT file
  // ==========================================================
  TFile* fOut = TFile::Open(outFile, "UPDATE");
  hMatrixRebinned->Write("hResponseMatrix");
  hCellGenRatio->Write();
  if (isSquare) {
    hPurity->Write();
    hStability->Write();
  }
  fOut->Close();

  std::cout << "[Done] Output saved to: " << outFile << std::endl;
}
