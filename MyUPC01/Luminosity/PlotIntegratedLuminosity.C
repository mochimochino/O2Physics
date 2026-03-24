#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include <TString.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

void PlotIntegratedLuminosity(const char* fileList = "/home/takuma/alice/O2Physics/MyUPC01/Luminosity/lumi.txt",
                              const char* histoPath = "eventselection-run3/luminosity/hLumiTCE")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::ifstream infile(fileList);
  if (!infile.is_open()) {
    std::cerr << "Error: Cannot open file list " << fileList << std::endl;
    return;
  }

  std::map<int, double> runLumiMap;
  std::map<int, std::string> runToPeriodMap;
  double totalLumi = 0.0;
  std::string rootFileName;
  while (infile >> rootFileName) {
    TFile* file = TFile::Open(rootFileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: Cannot open file " << rootFileName << std::endl;
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
      if (endPos != std::string::npos) {
        periodName = rootFileName.substr(pos + 5, endPos - pos - 5);
      }
    }

    TH1* hOriginal = (TH1*)file->Get(histoPath);
    if (!hOriginal) {
      std::cerr << "Warning: Cannot find histogram " << histoPath << " in " << rootFileName << std::endl;
      file->Close();
      continue;
    }

    int nBinsOrig = hOriginal->GetNbinsX();
    for (int i = 1; i <= nBinsOrig; ++i) {
      double lumi = hOriginal->GetBinContent(i);
      if (lumi > 0) {
        int runNumber = 0;
        const char* label = hOriginal->GetXaxis()->GetBinLabel(i);
        if (label && label[0] != '\0') {
          runNumber = std::stoi(label);
        } else {
          runNumber = (int)hOriginal->GetXaxis()->GetBinCenter(i);
        }
        runLumiMap[runNumber] += lumi;
        runToPeriodMap[runNumber] = periodName;
        totalLumi += lumi;
      }
    }
    file->Close();
  }

  int nValidRuns = runLumiMap.size();
  if (nValidRuns == 0) {
    std::cerr << "Error: No valid data points found in any files." << std::endl;
    return;
  }

  TH1D* hLumiHist = new TH1D("hLumiHist",
                             "Integrated luminosity as function of run number - apass5 (LumiTCE)",
                             nValidRuns, 0, nValidRuns);

  int bin = 1;
  for (auto const& [run, lumi] : runLumiMap) {
    hLumiHist->SetBinContent(bin, lumi);
    hLumiHist->GetXaxis()->SetBinLabel(bin, std::to_string(run).c_str());
    std::cout << "Mapped Run: " << run << " | Total L_int: " << lumi << std::endl;
    bin++;
  }

  std::cout << ">>> Total Integrated Luminosity: " << totalLumi << " <<< " << std::endl;

  TCanvas* c1 = new TCanvas("c1", "Luminosity Canvas", 1800, 600);

  c1->SetBottomMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.05);

  hLumiHist->GetXaxis()->SetTitle("Run Number");
  hLumiHist->GetYaxis()->SetTitle("Integrated Luminosity [1/#mub]");

  hLumiHist->GetXaxis()->SetTitleOffset(2.0);
  hLumiHist->GetYaxis()->SetTitleOffset(0.5);
  hLumiHist->GetXaxis()->LabelsOption("v");
  hLumiHist->GetXaxis()->SetLabelSize(0.035);

  hLumiHist->SetLineColor(kBlue);
  hLumiHist->SetLineWidth(1);

  hLumiHist->SetMaximum(hLumiHist->GetMaximum() * 1.2);

  hLumiHist->Draw("HIST");
  c1->Update();

  double ymax = gPad->GetUymax();

  TLine* line = new TLine();
  line->SetLineColor(kMagenta);
  line->SetLineWidth(2);

  TLatex* texPeriod = new TLatex();
  texPeriod->SetTextSize(0.04);
  texPeriod->SetTextColor(kMagenta);
  texPeriod->SetTextAlign(13);

  std::string currentPeriod = "";
  int currentBin = 1;

  for (auto const& [run, lumi] : runLumiMap) {
    std::string period = runToPeriodMap[run];
    if (period != currentPeriod && !period.empty()) {
      double x_line = currentBin - 1;

      line->DrawLine(x_line, 0, x_line, ymax);
      texPeriod->DrawLatex(x_line + 0.3, ymax * 0.96, period.c_str());

      currentPeriod = period;
    }
    currentBin++;
  }

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->SetTextAlign(31);

  latex->DrawLatex(0.90, 0.92, Form("Total L_{int} = %.2f /#mub", totalLumi));

  c1->SaveAs("hLumiTCE.png");
}
