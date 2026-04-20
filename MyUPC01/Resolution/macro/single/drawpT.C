#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>

#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

void drawpT(const char* filename = "AO2D_merged.root")
{
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return;
  }

  std::string treename = "";
  for (auto key : *f->GetListOfKeys()) {
    std::string name = key->GetName();
    if (name.find("DF_") == 0) {
      treename = name + "/O2udmcparticle";
      break;
    }
  }

  f->Close();

  if (treename.empty()) {
    std::cerr << "Error: Could not find DF_ directory in " << filename << std::endl;
    return;
  }

  std::cout << "Auto-detected tree path: " << treename << std::endl;

  ROOT::RDataFrame df(treename.c_str(), filename);

  auto df_pt = df.Define("pT", "sqrt(fPx*fPx + fPy*fPy)");

  auto h_pt = df_pt.Histo1D({"h_pT", "Muon Candidate Single p_{T} Distribution;p_{T} (GeV/c);Entries", 100, 0.3, 5.0}, "pT");

  auto c1 = new TCanvas("c1", "pT Canvas", 800, 600);
  c1->SetLeftMargin(0.13);

  h_pt->SetLineColor(kBlue);
  h_pt->SetLineWidth(2);
  h_pt->Draw();

  c1->SaveAs("pT_distribution.png");

  std::cout << "Successfully generated pT_distribution.png" << std::endl;
}
