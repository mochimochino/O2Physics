#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>

// ==========================================
// Main macro
// ==========================================
void pT_plot()
{
  // --- 初期設定 ---
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);

  // --- パラメータ設定 ---
  TString filePath = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/test/pT/AnalysisResults.root";
  TString dirName = "my-upc-mass-jpsi/";
  int rebinFactor = 3;
  double drawMin = 0.0;
  double drawMax = 2.5;

  // 1. ファイルとヒストグラムの読み込み
  TFile* f = TFile::Open(filePath, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open file: " << filePath << std::endl;
    return;
  }

  TH1* hUnlike = dynamic_cast<TH1*>(f->Get(dirName + "PtMuonUnlike"));
  TH1* hLike = dynamic_cast<TH1*>(f->Get(dirName + "PtMuonLike"));
  if (!hUnlike || !hLike) {
    std::cerr << "ERROR: Cannot get histograms from: " << dirName << std::endl;
    return;
  }

  // 2. Rebin
  hUnlike->Rebin(rebinFactor);
  hLike->Rebin(rebinFactor);

  // 3. Signal Extraction (Unlike - Like)
  TH1* hSignal = (TH1*)hUnlike->Clone("hSignal");
  // hSignal->Add(hLike, -1.0);

  double binWidth = hSignal->GetXaxis()->GetBinWidth(1);
  hSignal->SetTitle(Form("Signal p_{T} (Unlike); p_{T} [GeV/c]; Counts / (%.3f GeV/c)", binWidth));
  hSignal->SetLineColor(kBlack);
  hSignal->SetMarkerColor(kBlack);
  hSignal->SetLineWidth(2);
  hSignal->SetMarkerStyle(21);
  hSignal->GetXaxis()->SetRangeUser(drawMin, drawMax);

  // 4. 描画設定 (単一のキャンバス)
  TCanvas* c = new TCanvas("c", "pT Distribution Canvas", 800, 600);
  c->SetGrid();
  c->SetLogy(); // 対数スケールで見たい場合はコメントアウトを外す

  hSignal->Draw("E");

  TLegend* leg = new TLegend(0.55, 0.75, 0.90, 0.90);
  leg->SetBorderSize(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(hSignal, "Data (Unlike - Like)", "lep");
  leg->Draw();

  c->SaveAs("pT_distribution.png");
  f->Close();
}
