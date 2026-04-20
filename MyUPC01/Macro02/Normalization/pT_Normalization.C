#include "Math/MinimizerOptions.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>
#include <vector>

void pT_Normalization()
{
  // Setting parameter
  TString dataDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0408/Resolution/";
  TString histTag = "my-upc-muon-pair-resolution/registry/";
  TString histName = histTag + "hPairPtReco_PreCut";
  TString outDir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0416/CoherentFeedDown/";

  struct Sample {
    TString file;
    TString label;
  };

  std::vector<Sample> samples = {
    //{dataDir + "jpsi-incoh.root", "J/#psi incoh."},
    // {dataDir + "jpsi-coh.root", "J/#psi coh."},
    // {dataDir + "psi2s-incoh-fd.root", "#psi(2S) incoh. fd"},
    {dataDir + "psi2s-coh-fd.root", "#psi(2S) coh. fd"},
  };

  const int rebinFactor = 5;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  for (const auto& sample : samples) {
    TString inputFile = sample.file;
    std::cout << "\n======================================" << std::endl;
    std::cout << "[Processing] " << inputFile << std::endl;

    // File check
    TFile* fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
      std::cerr << "ERROR: Cannot open file: " << inputFile << std::endl;
      continue;
    }

    // Histograms
    TH1D* hOrigin = (TH1D*)fIn->Get(histName);
    if (!hOrigin) {
      std::cerr << "ERROR: Cannot find histogram: " << histName << " in " << inputFile << std::endl;
      fIn->Close();
      continue;
    }

    TString baseName = inputFile(inputFile.Last('/') + 1, inputFile.Length());
    baseName.ReplaceAll(".root", "");

    TString outHistName = "hTemplate_" + baseName;
    TH1D* hTemplate = (TH1D*)hOrigin->Clone(outHistName);

    if (!hTemplate->GetSumw2N()) {
      hTemplate->Sumw2();
    }
    hTemplate->SetDirectory(0);
    fIn->Close();

    // Normalization
    double integral = hTemplate->Integral();

    if (integral > 0) {
      hTemplate->Scale(1.0 / integral);
    } else {
      std::cerr << "ERROR: Integral is zero or negative for " << inputFile << std::endl;
      delete hTemplate;
      continue;
    }

    hTemplate->Rebin(rebinFactor);

    hTemplate->SetLineColor(kBlue + 1);
    hTemplate->SetLineWidth(2);
    hTemplate->SetMarkerStyle(20);
    hTemplate->SetMarkerColor(kBlue + 1);
    hTemplate->SetMarkerSize(0.8);

    hTemplate->GetXaxis()->SetTitle("p_{T}^{#mu#mu} (GeV/c)");
    hTemplate->GetYaxis()->SetTitle("Probability Density");
    hTemplate->GetXaxis()->SetRangeUser(0.0, 2.7);
    hTemplate->GetYaxis()->SetRangeUser(10e-5, 10e-1);
    hTemplate->GetXaxis()->SetTitleSize(0.045);
    hTemplate->GetYaxis()->SetTitleSize(0.045);
    hTemplate->GetXaxis()->SetLabelSize(0.04);
    hTemplate->GetYaxis()->SetLabelSize(0.04);
    hTemplate->GetYaxis()->SetTitleOffset(1.2);

    // Save histogram
    TString outFileName = outDir + "Template_" + baseName + ".root";
    TFile* fOut = TFile::Open(outFileName, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
      std::cerr << "ERROR: Cannot open output file: " << outFileName << std::endl;
      delete hTemplate;
      continue;
    }
    hTemplate->Write();
    fOut->Close();
    std::cout << "[Info] Saved ROOT : " << outFileName << std::endl;

    // Draw & Save PNG
    TCanvas* c = new TCanvas("c_" + baseName, "Normalization", 800, 600);
    c->SetLogy();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);

    c->cd();
    hTemplate->Draw("HIST E");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.65, 0.88, "MC, LHC26b8");
    latex.DrawLatex(0.65, 0.83, sample.label + " template");

    TString outImgName = outDir + "Norm_" + baseName + ".png";
    c->SaveAs(outImgName);
    std::cout << "[Info] Saved Image: " << outImgName << std::endl;

    delete c;
    delete hTemplate;
  }

  std::cout << "\n[Success] All processes completed!" << std::endl;
}
