#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include <TString.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

void PlotLumiCrossCheck(const char* fileList = "/home/takuma/alice/O2Physics/MyUPC01/Luminosity/lumi.txt",
                        const char* histoPathTCE = "eventselection-run3/luminosity/hLumiTCE",
                        const char* histoPathZNC = "eventselection-run3/luminosity/hLumiZNC")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::ifstream infile(fileList);
  if (!infile.is_open()) {
    std::cerr << "Error: Cannot open file list " << fileList << std::endl;
    return;
  }

  std::map<int, double> runLumiMapTCE;
  std::map<int, double> runLumiMapZNC;
  std::map<int, std::string> runToPeriodMap;

  double totalLumiTCE = 0.0;
  double totalLumiZNC = 0.0;
  std::string rootFileName;

  while (infile >> rootFileName) {
    TFile* file = TFile::Open(rootFileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
      if (file)
        file->Close();
      continue;
    }

    std::string periodName = "";
    size_t pos = rootFileName.find("LHC23");
    if (pos != std::string::npos) {
      size_t endPos = rootFileName.find("_", pos);
      if (endPos == std::string::npos)
        endPos = rootFileName.find("/", pos);
      if (endPos != std::string::npos)
        periodName = rootFileName.substr(pos + 5, endPos - pos - 5);
    }

    TH1* hOrigTCE = (TH1*)file->Get(histoPathTCE);
    TH1* hOrigZNC = (TH1*)file->Get(histoPathZNC);

    if (hOrigTCE && hOrigZNC) {
      int nBinsOrig = hOrigTCE->GetNbinsX();
      for (int i = 1; i <= nBinsOrig; ++i) {
        double lumiTCE = hOrigTCE->GetBinContent(i);
        double lumiZNC = hOrigZNC->GetBinContent(i);

        if (lumiTCE > 0) {
          int runNumber = 0;
          const char* label = hOrigTCE->GetXaxis()->GetBinLabel(i);
          if (label && label[0] != '\0')
            runNumber = std::stoi(label);
          else
            runNumber = (int)hOrigTCE->GetXaxis()->GetBinCenter(i);

          runLumiMapTCE[runNumber] += lumiTCE;
          runLumiMapZNC[runNumber] += lumiZNC;
          runToPeriodMap[runNumber] = periodName;

          totalLumiTCE += lumiTCE;
          totalLumiZNC += lumiZNC;
        }
      }
    } else {
      std::cerr << "Warning: Cannot find histograms in " << rootFileName << std::endl;
    }
    file->Close();
  }

  int nValidRuns = runLumiMapTCE.size();
  if (nValidRuns == 0)
    return;

  TH1D* hTCE = new TH1D("hTCE", "Luminosity Cross Check: TCE vs ZNC", nValidRuns, 0, nValidRuns);
  TH1D* hZNC = new TH1D("hZNC", "Luminosity Cross Check: TCE vs ZNC", nValidRuns, 0, nValidRuns);
  TH1D* hRatio = new TH1D("hRatio", "", nValidRuns, 0, nValidRuns);

  int bin = 1;
  for (auto const& [run, lumiTCE] : runLumiMapTCE) {
    double lumiZNC = runLumiMapZNC[run];

    hTCE->SetBinContent(bin, lumiTCE);
    hZNC->SetBinContent(bin, lumiZNC);

    if (lumiTCE > 0)
      hRatio->SetBinContent(bin, lumiZNC / lumiTCE);
    else
      hRatio->SetBinContent(bin, 0);

    std::string runStr = std::to_string(run);
    hTCE->GetXaxis()->SetBinLabel(bin, runStr.c_str());
    hZNC->GetXaxis()->SetBinLabel(bin, runStr.c_str());
    hRatio->GetXaxis()->SetBinLabel(bin, runStr.c_str());

    bin++;
  }

  TCanvas* c1 = new TCanvas("c1", "Luminosity Cross Check", 1800, 800);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetRightMargin(0.05);
  pad1->Draw();
  pad1->cd();

  hTCE->SetLineColor(kBlue);
  hTCE->SetLineWidth(2);
  hTCE->GetYaxis()->SetTitle("Integrated Luminosity [1/#mub]");
  hTCE->GetXaxis()->SetLabelSize(0);
  hTCE->SetMaximum(hTCE->GetMaximum() * 1.3);
  hTCE->GetYaxis()->SetTitleSize(20);
  hTCE->GetYaxis()->SetTitleFont(43);
  hTCE->GetYaxis()->SetTitleOffset(1.2);
  hTCE->GetYaxis()->SetLabelFont(43);
  hTCE->GetYaxis()->SetLabelSize(15);

  hZNC->SetLineColor(kRed);
  hZNC->SetLineWidth(2);

  hTCE->Draw("HIST");
  hZNC->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.80, 0.75, 0.93, 0.88);
  leg->AddEntry(hTCE, Form("TCE (Total: %.1f)", totalLumiTCE), "l");
  leg->AddEntry(hZNC, Form("ZNC (Total: %.1f)", totalLumiZNC), "l");
  leg->Draw();

  c1->Update();
  double ymax = gPad->GetUymax();
  TLine* line = new TLine();
  line->SetLineColor(kMagenta);
  line->SetLineWidth(2);
  TLatex* texPeriod = new TLatex();
  texPeriod->SetTextSize(0.04);
  texPeriod->SetTextColor(kMagenta);

  std::string currentPeriod = "";
  int currentBin = 1;
  for (auto const& [run, lumi] : runLumiMapTCE) {
    std::string period = runToPeriodMap[run];
    if (period != currentPeriod && !period.empty()) {
      double x_line = currentBin - 1;
      line->DrawLine(x_line, 0, x_line, ymax);
      texPeriod->DrawLatex(x_line + 0.3, ymax * 0.95, period.c_str());
      currentPeriod = period;
    }
    currentBin++;
  }

  c1->cd();
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin(0.05);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.8);

  hRatio->GetYaxis()->SetTitle("ZNC / TCE");
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(20);
  hRatio->GetYaxis()->SetTitleFont(43);
  hRatio->GetYaxis()->SetTitleOffset(1.2);
  hRatio->GetYaxis()->SetLabelFont(43);
  hRatio->GetYaxis()->SetLabelSize(15);
  hRatio->SetMinimum(0.8);
  hRatio->SetMaximum(1.2);

  hRatio->GetXaxis()->SetTitle("Run Number");
  hRatio->GetXaxis()->LabelsOption("v");
  hRatio->GetXaxis()->SetLabelSize(0.08);
  hRatio->GetXaxis()->SetTitleSize(20);
  hRatio->GetXaxis()->SetTitleFont(43);
  hRatio->GetXaxis()->SetTitleOffset(2.0);

  hRatio->Draw("P");

  TLine* lineRatio = new TLine(0, 1.0, nValidRuns, 1.0);
  lineRatio->SetLineColor(kBlue);
  lineRatio->SetLineStyle(2);
  lineRatio->Draw();

  c1->SaveAs("LumiCrossCheck_TCE_vs_ZNC.png");
}
