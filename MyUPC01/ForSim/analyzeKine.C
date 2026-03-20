#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>

int analyzeKine()
{
  std::string path{"o2sim_Kine.root"};
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open ROOT file " << path << "\n";
    return 1;
  }
  auto tree = (TTree*)file.Get("o2sim");
  if (!tree) {
    std::cerr << "Cannot find tree o2sim in file " << path << "\n";
    return 1;
  }

  std::vector<o2::MCTrack>* tracks{};
  tree->SetBranchAddress("MCTrack", &tracks);

  auto nEvents = tree->GetEntries();

  for (int i = 0; i < nEvents; i++) {

    cout << "Event " << i << endl;
    auto check = tree->GetEntry(i);

    for (int idxMCTrack = 0; idxMCTrack < tracks->size() - 1; ++idxMCTrack) {
      auto track = tracks->at(idxMCTrack);
      cout << track.GetPdgCode() << " " << track.GetPt() << " " << track.GetEta() << " " << track.getMotherTrackId() << " " << track.isTransported() << " "
           << track.getFirstDaughterTrackId() << endl;
    }
  }
  return 0;
}
