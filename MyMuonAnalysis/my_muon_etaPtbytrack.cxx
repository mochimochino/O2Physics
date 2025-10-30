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

        std::string nameEtaPt = "etaPtHistogram_Type" + std::to_string(type);
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
      
      const float eta = fwdtrack.eta();
      const float pt = fwdtrack.pt();
      const int type = fwdtrack.trackType();

      histos.fill(HIST("etaHistogram"), eta);
      histos.fill(HIST("ptHistogram"), pt);
      histos.fill(HIST("etaPtHistogram"), pt, eta);
      
      switch (type) {
        case 0:
            histos.fill(HIST("etaHistogram_Type0"), eta);
            histos.fill(HIST("ptHistogram_Type0"), pt);
            histos.fill(HIST("etaPtHistogram_Type0"), pt, eta);
            break;
        case 2:
            histos.fill(HIST("etaHistogram_Type2"), eta);
            histos.fill(HIST("ptHistogram_Type2"), pt);
            histos.fill(HIST("etaPtHistogram_Type2"), pt, eta);
            break;
        case 3:
            histos.fill(HIST("etaHistogram_Type3"), eta);
            histos.fill(HIST("ptHistogram_Type3"), pt);
            histos.fill(HIST("etaPtHistogram_Type3"), pt, eta);
            break;
        case 4:
            histos.fill(HIST("etaHistogram_Type4"), eta);
            histos.fill(HIST("ptHistogram_Type4"), pt);
            histos.fill(HIST("etaPtHistogram_Type4"), pt, eta);
            break;
        default:
            
            break;
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<my_muon_etaPtbytrack>(cfgc)};
}