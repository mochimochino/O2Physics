// ======================================================================
// FitMass_MC_ExtractDSCB.C
//
// MCデータからJ/psiのDSCBパラメータを抽出し、
// 必要なパラメータとデータ情報のみをプロット上に描画するマクロ
// ======================================================================

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>

// DSCB関数の定義
double DSCB(double* x, double* p)
{
  double m = x[0];
  double N = p[0];
  double m0 = p[1];
  double sigma = p[2];
  double alpha1 = p[3];
  double n1 = p[4];
  double alpha2 = p[5];
  double n2 = p[6];

  double alpha = (m - m0) / sigma;

  if (alpha < -alpha1) {
    double A = TMath::Power(n1 / TMath::Abs(alpha1), n1) * TMath::Exp(-0.5 * alpha1 * alpha1);
    double B = n1 / TMath::Abs(alpha1) - TMath::Abs(alpha1);
    return N * A * TMath::Power(B - alpha, -n1);
  } else if (alpha > alpha2) {
    double C = TMath::Power(n2 / TMath::Abs(alpha2), n2) * TMath::Exp(-0.5 * alpha2 * alpha2);
    double D = n2 / TMath::Abs(alpha2) - TMath::Abs(alpha2);
    return N * C * TMath::Power(D + alpha, -n2);
  } else {
    return N * TMath::Exp(-0.5 * alpha * alpha);
  }
}

void FitMass_MC_ExtractDSCB(
  // const char* filename = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0508/jpsi-coh.root",
  // const char* histoPath = "my-upc-muon-pair-resolution-eta/registry/hPairMassReco_PostCut",
  const char* filename = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/jpsi-incoh.root",
  const char* histoPath = "my-upc-muon-pair-resolution/registry/hPairMassReco_PostCut",
  int rebin = 2,
  const char* dataLabel = "J/#psi incoherent")
{
  // デフォルトのStatsボックスとFitボックスを非表示にする
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetErrorX(0.5);

  int fontCode = 42;
  gStyle->SetLabelFont(fontCode, "XYZ");
  gStyle->SetTitleFont(fontCode, "XYZ");

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Cannot open file." << std::endl;
    return;
  }

  TH1D* hMC = (TH1D*)f->Get(histoPath);
  if (!hMC) {
    std::cerr << "Error: Histogram not found." << std::endl;
    f->Close();
    return;
  }
  hMC->SetDirectory(0);
  f->Close();

  // Rebin処理
  if (rebin > 1) {
    hMC->Rebin(rebin);
  }

  hMC->SetTitle(""); // ヒストグラムタイトルを非表示

  // 軸タイトルの設定
  double binWidthMeV = hMC->GetBinWidth(1) * 1000.0;
  hMC->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/#it{c}^{2})");
  hMC->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/#it{c}^{2})", binWidthMeV));
  hMC->GetYaxis()->SetTitleOffset(1.3);

  // フィット範囲の設定
  double fitMin = 2.0;
  double fitMax = 4.0;
  hMC->GetXaxis()->SetRangeUser(2.0, 4.0);
  hMC->GetYaxis()->SetRangeUser(0, hMC->GetMaximum() * 1.5);

  // TF1のセットアップとパラメータの制限
  TF1* fJpsi = new TF1("fJpsi", DSCB, fitMin, fitMax, 7);
  fJpsi->SetNpx(1000);
  fJpsi->SetLineColor(kRed);
  fJpsi->SetLineWidth(2);

  // 初期値と制限 (Limits) の設定
  fJpsi->SetParameter(0, hMC->GetMaximum());
  fJpsi->SetParameter(1, 3.096);
  fJpsi->SetParLimits(1, 3.0, 3.2); // m0
  fJpsi->SetParameter(2, 0.07);
  fJpsi->SetParLimits(2, 0.03, 0.15); // sigma
  fJpsi->SetParameter(3, 1.0);
  fJpsi->SetParLimits(3, 0.5, 5.0); // alpha1
  fJpsi->SetParameter(4, 10.0);
  fJpsi->SetParLimits(4, 1.0, 50.0); // n1
  fJpsi->SetParameter(5, 2.0);
  fJpsi->SetParLimits(5, 0.5, 5.0); // alpha2
  fJpsi->SetParameter(6, 10.0);
  fJpsi->SetParLimits(6, 1.0, 50.0); // n2

  // フィット実行
  hMC->Fit(fJpsi, "L R S Q"); // Qを追加してターミナル出力を少し静かにする

  // ==========================================
  // 描画
  // ==========================================
  TCanvas* c1 = new TCanvas("c1", "MC J/psi DSCB Fit", 900, 650);
  c1->SetLeftMargin(0.13);
  c1->SetBottomMargin(0.13);
  c1->SetTopMargin(0.06);
  c1->SetRightMargin(0.05);
  c1->SetGrid();

  hMC->SetMarkerStyle(20);
  hMC->SetMarkerSize(0.8);
  hMC->SetMarkerColor(kBlack);
  hMC->SetLineColor(kBlack);
  hMC->Draw("E");
  fJpsi->Draw("SAME");

  // ==========================================
  // カスタムテキストの描画 (データ情報とFit結果)
  // ==========================================
  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(fontCode);

  tex.SetTextSize(0.036);
  tex.DrawLatex(0.16, 0.88, "#bf{Invariant Mass (Unlike)}");
  tex.DrawLatex(0.16, 0.83, "#bf{LHC26b8} #font[52]{(MC, UPC #mu pair)}");
  tex.DrawLatex(0.16, 0.78, Form("Data: %s", dataLabel));
  // tex.DrawLatex(0.16, 0.73, "StandaloneMuon");
  tex.DrawLatex(0.16, 0.73, "GlobalMuon");

  double chi2ndf = fJpsi->GetChisquare() / fJpsi->GetNDF();
  tex.SetTextSize(0.032);
  tex.DrawLatex(0.73, 0.55, Form("#chi^{2}/ndf = %.2f", chi2ndf));
  tex.DrawLatex(0.73, 0.90, Form("m_{0} = %.4f GeV/#it{c}^{2}", fJpsi->GetParameter(1)));
  tex.DrawLatex(0.73, 0.85, Form("#sigma = %.4f GeV/#it{c}^{2}", fJpsi->GetParameter(2)));
  tex.DrawLatex(0.73, 0.80, Form("#alpha_{1} = %.3f", fJpsi->GetParameter(3)));
  tex.DrawLatex(0.73, 0.75, Form("n_{1} = %.3f", fJpsi->GetParameter(4)));
  tex.DrawLatex(0.73, 0.70, Form("#alpha_{2} = %.3f", fJpsi->GetParameter(5)));
  tex.DrawLatex(0.73, 0.65, Form("n_{2} = %.3f", fJpsi->GetParameter(6)));

  // ==========================================
  // ターミナル出力
  // ==========================================
  std::cout << "\n==================================================" << std::endl;
  std::cout << " [MC Parameters Extracted for Data Fit] " << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << Form("double mc_m0     = %.5f;", fJpsi->GetParameter(1)) << std::endl;
  std::cout << Form("double mc_sigma  = %.5f;", fJpsi->GetParameter(2)) << std::endl;
  std::cout << Form("double mc_alpha1 = %.5f;", fJpsi->GetParameter(3)) << std::endl;
  std::cout << Form("double mc_n1     = %.5f;", fJpsi->GetParameter(4)) << std::endl;
  std::cout << Form("double mc_alpha2 = %.5f;", fJpsi->GetParameter(5)) << std::endl;
  std::cout << Form("double mc_n2     = %.5f;", fJpsi->GetParameter(6)) << std::endl;
  std::cout << "==================================================\n"
            << std::endl;

  c1->SaveAs("MC_DSCB_Extraction_Custom.png");
}
