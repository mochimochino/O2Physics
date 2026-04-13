// ============================================================
//  Step3_2_TUnfold.C (TUnfold Optimization and Validation)
// ============================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"
#include "TSpline.h"
#include "TStyle.h"

// TUnfold headers
#include "TUnfoldDensity.h"

#include <iostream>

void Step3_2_TUnfold()
{
  // --- Settings ---
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0409/IncoherentJpsi/50/";
  const TString inFile = inDir + "Step2_Response_for_Unfolding.root";
  const TString outDir = inDir;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open " << inFile << std::endl;
    return;
  }

  TH1D* hGen = (TH1D*)fIn->Get("hMCGen_Rebinned");     // Truth (MC Gen)
  TH2D* hMat = (TH2D*)fIn->Get("hResponseMatrix");     // Response Matrix
  TH1D* hData = (TH1D*)fIn->Get("hDataReco_Rebinned"); // Measured (Data)

  std::cout << "[Info] Constructing TUnfoldDensity..." << std::endl;

  // TUnfoldの設定
  // kHistMapOutputVert: Y軸がTruth（Unfolded）, X軸がMeasured（Reco）であることを前提とします
  // 論文に合わせて正則化条件は2階微分（曲率: kRegModeCurvature）を指定します
  TUnfoldDensity unfold(hMat, TUnfold::kHistMapOutputVert,
                        TUnfold::kRegModeCurvature,
                        TUnfold::kEConstraintNone,
                        TUnfoldDensity::kDensityModeNone);

  // 測定データのセット
  unfold.SetInput(hData);

  // ==========================================================
  // Tau (τ) の最適化: 平均グローバル相関係数の最小化
  // ==========================================================
  std::cout << "[Info] Scanning tau to minimize average global correlation..." << std::endl;
  TSpline* scanResult = 0;

  // ScanTau(探索点数, 最小τ, 最大τ, スプラインへのポインタ, スキャンモード)
  // kEScanTauRhoAvg は平均グローバル相関係数を最小化するモードです
  Int_t iBest = unfold.ScanTau(100, 0.0001, 1.0, &scanResult, TUnfoldDensity::kEScanTauRhoAvg);
  Double_t bestTau = unfold.GetTau();

  std::cout << "[Info] Best tau found: " << bestTau << std::endl;

  // 最適化された展開結果の取得
  TH1D* hUnfolded = (TH1D*)unfold.GetOutput("hUnfolded_TUnfold");
  hUnfolded->SetLineColor(kRed + 1);
  hUnfolded->SetMarkerColor(kRed + 1);
  hUnfolded->SetMarkerStyle(20);

  // ==========================================================
  // 描画セットアップ
  // ==========================================================
  TCanvas* c1 = new TCanvas("c1", "TUnfold Results", 1200, 500);
  c1->Divide(2, 1);

  // ----------------------------------------------------------
  // Pad 1: Tau スキャンの可視化 (ρ vs log(τ))
  // ----------------------------------------------------------
  c1->cd(1);
  TPad* padScan = (TPad*)gPad;
  padScan->SetLeftMargin(0.15);
  padScan->SetRightMargin(0.05);
  padScan->SetGrid();

  if (scanResult) {
    scanResult->SetTitle("Global Correlation Scan;#log_{10}(#tau);Average Global Correlation #langle#rho#rangle");
    scanResult->SetLineWidth(2);
    scanResult->SetLineColor(kBlue + 1);
    scanResult->Draw("C"); // 曲線として描画

    // 最小値（最適τ）にマーカーを打つ
    Double_t bestTauLog = std::log10(bestTau);
    Double_t minRho = scanResult->Eval(bestTauLog);

    TMarker* mBest = new TMarker(bestTauLog, minRho, 29); // Star marker
    mBest->SetMarkerColor(kRed);
    mBest->SetMarkerSize(2.5);
    mBest->Draw("SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.8, Form("Best #tau = %.4f", bestTau));
  }

  // ----------------------------------------------------------
  // Pad 2: 展開結果と比の描画 (上段・下段分割)
  // ----------------------------------------------------------
  c1->cd(2);
  TPad* padRight = (TPad*)gPad;
  padRight->Divide(1, 2);

  // Pad 2-1: Spectra
  TPad* padTop = (TPad*)padRight->GetPad(1);
  padTop->SetPad(0.0, 0.3, 1.0, 1.0);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.12);
  padTop->SetRightMargin(0.05);
  padTop->SetLogy();
  padTop->Draw();
  padTop->cd();

  TH1D* hGenDraw = (TH1D*)hGen->Clone("hGenDraw");
  hGenDraw->SetTitle("TUnfold Result (Optimized #tau)");
  hGenDraw->SetLineColor(kBlack);
  hGenDraw->SetLineWidth(2);
  hGenDraw->GetXaxis()->SetLabelSize(0);
  hGenDraw->GetXaxis()->SetTitleSize(0);
  hGenDraw->GetYaxis()->SetTitle("Counts");
  hGenDraw->GetYaxis()->SetTitleSize(0.05);
  hGenDraw->GetYaxis()->SetLabelSize(0.045);

  double maxY = std::max({hGenDraw->GetMaximum(), hData->GetMaximum(), hUnfolded->GetMaximum()});
  double minY = hGenDraw->GetMinimum();
  hGenDraw->SetMaximum(maxY * 10.0);
  hGenDraw->SetMinimum(minY * 0.1);
  hGenDraw->Draw("HIST");

  TH1D* hDataDraw = (TH1D*)hData->Clone("hDataDraw");
  hDataDraw->SetLineColor(kBlue + 1);
  hDataDraw->SetMarkerColor(kBlue + 1);
  hDataDraw->SetMarkerStyle(24);
  hDataDraw->Draw("PE SAME");
  hUnfolded->Draw("PE SAME");

  TLegend* leg = new TLegend(0.60, 0.65, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hGenDraw, "Truth", "l");
  leg->AddEntry(hDataDraw, "Measured", "pe");
  leg->AddEntry(hUnfolded, "TUnfold", "pe");
  leg->Draw();

  // Pad 2-2: Ratio
  padRight->cd();
  TPad* padBot = (TPad*)padRight->GetPad(2);
  padBot->SetPad(0.0, 0.0, 1.0, 0.3);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.12);
  padBot->SetRightMargin(0.05);
  padBot->Draw();
  padBot->cd();

  TH1D* hRatio = (TH1D*)hUnfolded->Clone("hRatio");
  hRatio->Divide(hGen);
  hRatio->SetTitle("");
  hRatio->GetYaxis()->SetTitle("Unfolded/Truth");
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
  line->Draw("SAME");

  c1->SaveAs(outDir + "Step3_2_TUnfold_Optimization.png");
  std::cout << "[Done] Step3-2: TUnfold Tau optimization plot created." << std::endl;
}
