#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>

// -------- Make2DHistFromBinnedHists との連携 --------
#include "../Make2DHistFromBinnedHists.C"

/**
 * @brief Performs equal-width rebinning on Gen and Reco histograms and saves as PNG.
 */
void RebinEqualWidth(
  const TString inputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0415/Coherent/Step1_merged.root",
  const TString outputFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0415/Coherent/equalwidth/4bin/Step2_Rebinned.root",
  double binWidth = 0.016,
  double xMin = 0.0,
  double xMax = 0.065)
{
  // --- Style ---
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open input file: " << inputFile << std::endl;
    return;
  }

  TH1D* hGen = (TH1D*)fIn->Get("hGenPt2");
  TH1D* hReco = (TH1D*)fIn->Get("hRecoPt2");

  if (!hGen || !hReco) {
    std::cerr << "[Error] hGenPt2 or hRecoPt2 not found in " << inputFile << std::endl;
    fIn->Close();
    return;
  }

  // Calculate number of bins
  int nBins = (int)((xMax - xMin) / binWidth);
  double xMaxAdjusted = xMin + nBins * binWidth;

  std::cout << "[Info] Rebinning to " << nBins << " bins with width " << binWidth
            << " from " << xMin << " to " << xMaxAdjusted << std::endl;

  // Create new histograms
  TH1D* hGenRebin = new TH1D("hGenPt2_rebin", "Gen p_{T}^{2} (Rebinned)", nBins, xMin, xMaxAdjusted);
  TH1D* hRecoRebin = new TH1D("hRecoPt2_rebin", "Reco p_{T}^{2} (Rebinned)", nBins, xMin, xMaxAdjusted);

  hGenRebin->GetXaxis()->SetTitle(hGen->GetXaxis()->GetTitle());
  hGenRebin->GetYaxis()->SetTitle(hGen->GetYaxis()->GetTitle());
  hRecoRebin->GetXaxis()->SetTitle(hReco->GetXaxis()->GetTitle());
  hRecoRebin->GetYaxis()->SetTitle(hReco->GetYaxis()->GetTitle());

  // Fill loop (Summing contents and errors)
  auto fillRebin = [&](TH1D* hOld, TH1D* hNew) {
    for (int i = 1; i <= hOld->GetNbinsX(); ++i) {
      double x = hOld->GetBinCenter(i);
      double val = hOld->GetBinContent(i);
      double err = hOld->GetBinError(i);
      if (x >= xMin && x < xMaxAdjusted) {
        int targetBin = hNew->FindBin(x);
        hNew->SetBinContent(targetBin, hNew->GetBinContent(targetBin) + val);
        hNew->SetBinError(targetBin, std::sqrt(std::pow(hNew->GetBinError(targetBin), 2) + std::pow(err, 2)));
      }
    }
  };

  fillRebin(hGen, hGenRebin);
  fillRebin(hReco, hRecoRebin);

  // --- Save to ROOT ---
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  hGenRebin->Write();
  hRecoRebin->Write();
  fOut->Close();
  std::cout << "[Info] Saved ROOT: " << outputFile << std::endl;

  // --- Draw and Save PNG ---
  TString outBase = outputFile;
  outBase.ReplaceAll(".root", "");

  auto drawAndSave = [&](TH1D* h, const TString& title, const TString& suffix) {
    TCanvas* c = new TCanvas("c_" + suffix, title, 800, 600);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetGrid();

    h->SetLineColor(kBlue + 1);
    h->SetLineWidth(2);
    h->Draw("HIST E");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.04);
    // tex.DrawLatex(0.65, 0.84, Form("#bf{%s}", title.Data()));
    tex.DrawLatex(0.65, 0.79, Form("Bin Width: %.3f", binWidth));

    c->SaveAs(outBase + "_" + suffix + ".png");
    delete c;
  };

  drawAndSave(hGenRebin, "Gen p_{T}^{2} (Rebinned)", "Gen");
  drawAndSave(hRecoRebin, "Reco p_{T}^{2} (Rebinned)", "Reco");

  std::cout << "[Info] Saved PNGs: " << outBase << "_Gen/Reco.png" << std::endl;

  fIn->Close();

  // -------- 2D ヒストグラム作成 --------
  // リビン済みヒストのビン境界を使って 2D レスポンス行列をリビンし、
  // (Diag + Gen) / Gen [%] を出力する。
  //
  // matrixFile: 高解像度 2D 行列の ROOT ファイル (inputFile と同じでよい場合が多い)
  // 必要に応じて下記を書き換えてください。
  const TString matrixFile = inputFile;                // 高解像度行列のファイル
  const TString matrixName = "hResponseMatrixPairPt2"; // 行列ヒスト名
  const TString outFile2D = outputFile;                // 2D 結果の出力先 (Rebin出力と同じでも可)

  Make2DHistFromBinnedHists(
    matrixFile,       // 行列ファイル
    matrixName,       // 行列ヒスト名
    outputFile,       // リビン済み 1D ヒストが入ったファイル
    "hGenPt2_rebin",  // X軸（Gen）用ヒスト名
    "hRecoPt2_rebin", // Y軸（Reco）用ヒスト名
    outFile2D         // 2D 出力 ROOT ファイル
  );
}
