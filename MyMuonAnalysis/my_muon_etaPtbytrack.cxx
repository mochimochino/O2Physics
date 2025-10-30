// FWD Muon Pt Histogram

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include <array>
#include <string>

using namespace o2;
using namespace o2::framework;

struct my_muon_etaPtbytrack {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 30, "N bins in Eta histo"};
  Configurable<float> minEta{"minEta", -5.0, "min Eta in Eta histo"};
  Configurable<float> maxEta{"maxEta", -2.0, "max Eta in Eta histo"};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0, "min pT in pT histo"};
  Configurable<float> maxPt{"maxPt", 10, "max pT in pT histo"};


  const std::array<int, 4> trackTypes = {0, 2, 3, 4};


  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
    const AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T}"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta}, true);
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt}, true);
    histos.add("etaPtHistogram", "#eta vs pT", kTH2F, {axisPt, axisEta});
  
    // each track types

    for (int type : trackTypes) {
        std::string nameEta = "etaHistogram_Type" + std::to_string(type);
        std::string titleEta = "Eta (Type" + std::to_string(type) +  ")";

        std::string namePt = "ptHistogram_Type" + std::to_string(type);
        std::string titlePt = "pT (Type" + std::to_string(type) + ")";

        std::string nameEtaPt = "EtaPtHistogram_Type" + std::to_string(type);
        std::string titleEtaPt = "Eta vs pT (Type" + std::to_string(type) + ")";


        // create histogram
        histos.add(nameEta.c_str(), titleEta.c_str(), kTH1F, {axisEta}, true);
        histos.add(namePt.c_str(), titlePt.c_str(), kTH1F, {axisPt}, true);
        histos.add(nameEtaPt.c_str(), titleEtaPt.c_str(), kTH2F, {axisPt, axisEta});
    }
  
  }



// have to change to fwd tracks
// Tutorial TracksIU -> FwdTrack
  void process(aod::FwdTracks const& fwdtracks)
  {
    for (auto& fwdtrack : fwdtracks) {
      histos.fill(HIST("etaHistogram"), fwdtrack.eta());
      histos.fill(HIST("ptHistogram"), fwdtrack.pt());
      histos.fill(HIST("etaPtHistogram"), fwdtrack.pt(), fwdtrack.eta());
      
      // type
      int type = fwdtrack.trackType();

      if (type == 0 || type == 2 || type == 3 || type == 4) {

          std::string nameEta   = "etaHistogram_Type" + std::to_string(type);
          std::string namePt    = "ptHistogram_Type" + std::to_string(type);
          std::string nameEtaPt = "etaPtHistogram_Type" + std::to_string(type);

          histos.fill(HIST("nameEta"), fwdtrack.eta());
          histos.fill(HIST("namePt"), fwdtrack.pt());
          histos.fill(HIST("nameEtaPt"), fwdtrack.pt(), fwdtrack.eta());
        }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<my_muon_etaPtbytrack>(cfgc)};
}
