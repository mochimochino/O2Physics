// FWD Muon Pt Histogram
// Using MymuonCut
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

//#include "PWGDQ/Core/CutsLibrary.h"

#include "MymuonCut.h"

using namespace o2;
using namespace o2::framework;

struct my_muonPt {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 30, "N bins in Eta histo"};
  Configurable<float> minEta{"minEta", -5.0, "min Eta in Eta histo"};
  Configurable<float> maxEta{"maxEta", -2.0, "max Eta in Eta histo"};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0, "min pT in pT histo"};
  Configurable<float> maxPt{"maxPt", 10, "max pT in pT histo"};

  Configurable<int> nBinsPDca{"nBinsPDca", 200, "N bins in PDca histo"};
  Configurable<float> minPDca{"minPDca", 0.0, "min PDca in PDca histo"};
  Configurable<float> maxPDca{"maxPDca", 1000.0, "max PDca in PDca histo"};

// add PWGDQ info
  Configurable<int> nBinsRAtAbsorberEnd{"nBinsRAtAbsorberEnd", 200, "N bins in RAtAbsorberEnd histo"};
  Configurable<float> minRAtAbsorberEnd{"minRAtAbsorberEnd", 0.0, "min RAtAbsorberEnd"};
  Configurable<float> maxRAtAbsorberEnd{"maxRAtAbsorberEnd", 100.0, "max RAtAbsorberEnd"};


  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
    const AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T}"};
    const AxisSpec axisPDca{nBinsPDca, minPDca, maxPDca, "PDca muon"};
    const AxisSpec axisRAtAbsorberEnd{nBinsRAtAbsorberEnd, minRAtAbsorberEnd, maxRAtAbsorberEnd, "R at Absorber End"};


    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta}, true);
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt}, true);
    histos.add("etaPtHistogram", "#eta vs pT", kTH2F, {axisPt, axisEta});
    histos.add("PDcamuonHistogram", "PDca muon Histogram", kTH1F, {axisPDca}, true);
    histos.add("RAtAbsorberEndHistogram", "R at Absorber End", kTH1F, {axisRAtAbsorberEnd}, true);
  }

// have to change to fwd tracks
// Tutorial TracksIU -> FwdTrack
  void process(aod::FwdTracks const& fwdtracks)
  {
    for (auto& fwdtrack : fwdtracks) {
      if (!isGoodMuonTrack(fwdtrack)) {
        continue;
      }
      histos.fill(HIST("etaHistogram"), fwdtrack.eta());
      histos.fill(HIST("ptHistogram"), fwdtrack.pt());
      histos.fill(HIST("etaPtHistogram"), fwdtrack.pt(), fwdtrack.eta());
      histos.fill(HIST("PDcamuonHistogram"), fwdtrack.pDca());
      histos.fill(HIST("RAtAbsorberEndHistogram"), fwdtrack.rAtAbsorberEnd());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<my_muonPt>(cfgc)};
}
