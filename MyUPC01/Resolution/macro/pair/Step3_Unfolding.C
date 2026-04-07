// ============================================================
//  Step3_Unfolding.C (Closure Test with Bayes & SVD)
// ============================================================
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TStyle.h"

#include <iostream>

// RooUnfold headers
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include <TSystem.h>

void Step3_Unfolding()
{
  // ----------------------------------------------------------
  // 0. RooUnfoldライブラリの動的ロード
  // ----------------------------------------------------------
  if (gSystem->Load("libRooUnfold") < 0) {
    std::cerr << "[Error] Failed to load libRooUnfold.so! Please check your environment." << std::endl;
    return;
  }
  std::cout << "[Info] Successfully loaded libRooUnfold." << std::endl;

  // ----------------------------------------------------------
  // Settings
  // ----------------------------------------------------------
  const TString inDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0403/without/";
  const TString inFile = inDir + "Step2_Response_for_Unfolding.root";
  const TString outDir = inDir;

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);

  TFile* fIn = TFile::Open(inFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[Error] Cannot open " << inFile << std::endl;
    return;
  }

  // ----------------------------------------------------------
  // 1. ヒストグラムの読み込み
  // ----------------------------------------------------------
  TH1D* hGen = (TH1D*)fIn->Get("hMCGen_Rebinned");     // True distribution
  TH1D* hReco = (TH1D*)fIn->Get("hMCReco_Rebinned");   // Measured distribution (MC)
  TH2D* hMat = (TH2D*)fIn->Get("hResponseMatrix");     // Response Matrix
  TH1D* hData = (TH1D*)fIn->Get("hDataReco_Rebinned"); // Data to unfold (currently = hReco)

  if (!hGen || !hReco || !hMat || !hData) {
    std::cerr << "[Error] Required histograms not found in " << inFile << std::endl;
    return;
  }

  // ----------------------------------------------------------
  // 2. RooUnfoldResponse の構築
  // ----------------------------------------------------------
  std::cout << "[Info] Constructing RooUnfoldResponse..." << std::endl;
  // 引数: (Measured 1D, True 1D, Response Matrix 2D)
  RooUnfoldResponse response(hReco, hGen, hMat);

  // ----------------------------------------------------------
  // 3. Method A: Iterative Bayesian (D'Agostini)
  // ----------------------------------------------------------
  int nIterations = 4; // ベイズ展開の反復回数（通常 3〜5 あたりで最適化します）
  std::cout << "[Info] Running Bayesian Unfolding (iterations = " << nIterations << ")..." << std::endl;
  RooUnfoldBayes unfoldBayes(&response, hData, nIterations);

  TH1D* hUnfoldedBayes = (TH1D*)unfoldBayes.Hreco();
  hUnfoldedBayes->SetName("hUnfoldedBayes");
  hUnfoldedBayes->SetLineColor(kRed);
  hUnfoldedBayes->SetMarkerColor(kRed);
  hUnfoldedBayes->SetMarkerStyle(20);

  // ----------------------------------------------------------
  // 4. Method B: SVD (Singular Value Decomposition)
  // ----------------------------------------------------------
  int kTerm = hGen->GetNbinsX() / 2; // SVDの正則化パラメータ（ビンの半数程度が初期値の目安）
  std::cout << "[Info] Running SVD Unfolding (k-term = " << kTerm << ")..." << std::endl;
  RooUnfoldSvd unfoldSvd(&response, hData, kTerm);

  TH1D* hUnfoldedSvd = (TH1D*)unfoldSvd.Hreco();
  hUnfoldedSvd->SetName("hUnfoldedSvd");
  hUnfoldedSvd->SetLineColor(kBlue);
  hUnfoldedSvd->SetMarkerColor(kBlue);
  hUnfoldedSvd->SetMarkerStyle(21);

  // ----------------------------------------------------------
  // 5. 描画: Closure Test (Unfolded vs Gen)
  // ----------------------------------------------------------
  hGen->SetLineColor(kBlack);
  hGen->SetLineWidth(2);

  TCanvas* c1 = new TCanvas("c1", "Closure Test", 800, 800);
  c1->Divide(1, 2);

  // 上段: 物理分布の比較
  TPad* pad1 = (TPad*)c1->GetPad(1);
  pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  hGen->SetTitle("Closure Test: Unfolded MC vs Gen MC");
  hGen->GetXaxis()->SetLabelSize(0);
  hGen->GetXaxis()->SetTitleSize(0);
  hGen->Draw("HIST");
  hUnfoldedBayes->Draw("PE SAME");
  hUnfoldedSvd->Draw("PE SAME");

  TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->AddEntry(hGen, "True (MC Gen)", "l");
  leg->AddEntry(hUnfoldedBayes, Form("Bayes (iter=%d)", nIterations), "pe");
  leg->AddEntry(hUnfoldedSvd, Form("SVD (k=%d)", kTerm), "pe");
  leg->Draw();

  c1->cd();

  // 下段: Ratio (Unfolded / Gen)
  TPad* pad2 = (TPad*)c1->GetPad(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1D* hRatioBayes = (TH1D*)hUnfoldedBayes->Clone("hRatioBayes");
  hRatioBayes->Divide(hGen);
  hRatioBayes->SetTitle("");
  hRatioBayes->GetYaxis()->SetTitle("Unfolded / True");
  hRatioBayes->GetYaxis()->SetRangeUser(0.1, 1.9); // 比が1.0に収まっているかを確認
  hRatioBayes->GetYaxis()->SetNdivisions(505);
  hRatioBayes->GetYaxis()->SetLabelSize(0.1);
  hRatioBayes->GetYaxis()->SetTitleSize(0.12);
  hRatioBayes->GetYaxis()->SetTitleOffset(0.4);

  hRatioBayes->GetXaxis()->SetTitle("p_{T}^{2} (GeV^{2}/c^{2})");
  hRatioBayes->GetXaxis()->SetLabelSize(0.1);
  hRatioBayes->GetXaxis()->SetTitleSize(0.12);

  hRatioBayes->Draw("PE");

  TH1D* hRatioSvd = (TH1D*)hUnfoldedSvd->Clone("hRatioSvd");
  hRatioSvd->Divide(hGen);
  hRatioSvd->Draw("PE SAME");

  TLine* line = new TLine(hGen->GetXaxis()->GetXmin(), 1.0, hGen->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw("SAME");

  c1->SaveAs(outDir + "Step3_ClosureTest.png");

  std::cout << "[Done] Closure test complete. Check Step3_ClosureTest.png" << std::endl;
}
