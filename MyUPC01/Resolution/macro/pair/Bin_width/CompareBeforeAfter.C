#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include <iostream>

/**
 * @brief リビン前後の1Dヒストグラムを比較するマクロ
 */
void CompareBeforeAfter(
  const TString originalFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0414/Coherent/Step1_merged.root",
  const TString rebinnedFile = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0414/test/Step2_Rebinned_EqualWidth.root",
  const TString histName = "hGenPt2" // "hGenPt2" or "hRecoPt2"
) {
  gStyle->SetOptStat(0);
  
  TFile* fOrig = TFile::Open(originalFile, "READ");
  TFile* fRebin = TFile::Open(rebinnedFile, "READ");
  
  if (!fOrig || !fRebin) return;

  TH1D* hOrig = (TH1D*)fOrig->Get(histName);
  // リビン後のファイルでは名前が "_rebin" などになっている可能性があるため調整
  TString rebinName = histName;
  if (!fRebin->Get(rebinName)) {
      if (fRebin->Get(histName + "_rebin")) rebinName = histName + "_rebin";
      else if (fRebin->Get(histName + "_flat")) rebinName = histName + "_flat";
      else if (fRebin->Get(histName + "_cdiff")) rebinName = histName + "_cdiff";
  }
  TH1D* hRebin = (TH1D*)fRebin->Get(rebinName);

  if (!hOrig || !hRebin) {
      std::cerr << "Histograms not found: " << histName << " or " << rebinName << std::endl;
      return;
  }

  TCanvas* c = new TCanvas("c_comp", "Before vs After Rebinning", 800, 600);
  c->SetLogy();
  c->SetGrid();

  // オリジナル（細いビン）
  hOrig->SetLineColor(kGray+1);
  hOrig->SetFillColorAlpha(kGray, 0.3);
  hOrig->SetTitle("Comparison: Before vs After Rebinning");
  hOrig->Draw("HIST");

  // リビン後（太いビン）
  // 面積を合わせる（積分値を等しくする）ために、正規化はせずそのまま描画
  // ただし、ビン幅が違うので、見た目上のカウントは増えることに注意
  hRebin->SetLineColor(kRed+1);
  hRebin->SetLineWidth(2);
  hRebin->SetMarkerStyle(20);
  hRebin->SetMarkerSize(0.8);
  hRebin->SetMarkerColor(kRed+1);
  hRebin->Draw("SAME P E");

  TLegend* leg = new TLegend(0.5, 0.7, 0.88, 0.88);
  leg->AddEntry(hOrig, "Original (Fine bins)", "f");
  leg->AddEntry(hRebin, "Rebinned", "lp");
  leg->Draw();

  TString outName = rebinnedFile;
  outName.ReplaceAll(".root", "_Comparison_" + histName + ".png");
  c->SaveAs(outName);
}
