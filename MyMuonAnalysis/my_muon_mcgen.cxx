#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

struct my_muon_mcgen {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 30, "N bins in Eta histo"};
  Configurable<float> minEta{"minEta", -5.0, "min Eta in Eta histo"};
  Configurable<float> maxEta{"maxEta", -2.0, "max Eta in Eta histo"};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0, "min pT in pT histo"};
  Configurable<float> maxPt{"maxPt", 6, "max pT in pT histo"};

  //Configurable<int> nBinsPhi{"nBinsPhi", 100, "N bins in Phi histo"};
  //Configurable<float> minPhi{"minPhi", -M_PI, "min pT in Phi histo"};
  //Configurable<float> maxPhi{"maxPhi", M_PI, "max pT in Phi histo"};

void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
    const AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T}"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta}, true);
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt}, true);
    histos.add("etaPtHistogram", "#eta vs pT", kTH2F, {axisPt, axisEta});
  }

// have to change to fwd tracks
// Tutorial TracksIU -> FwdTrack
void process(aod::McParticles const& mcParticles)
  {
  for (auto& mc : mcParticles) {
    if (abs(mc.pdgCode()) != 13) continue; // muon
    
    const float eta = mc.eta();
    const float pt = mc.pt();

    if (eta < minEta || eta > maxEta || pt < minPt || pt > maxPt) continue;

    histos.fill(HIST("etaHistogram"), eta);
    histos.fill(HIST("ptHistogram"), pt);
    histos.fill(HIST("etaPtHistogram"), pt, eta);
   }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) {
  return WorkflowSpec{
    adaptAnalysisTask<my_muon_mcgen>(cfg)
  };
}
