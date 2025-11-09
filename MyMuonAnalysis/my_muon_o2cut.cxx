// FWD Muon Pt Histogram
// Using MymuonCut
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/AnalysisCut.h"

//#include "MymuonCut.h"

using namespace o2;
using namespace o2::framework;

using o2::aod::dqcuts::GetAnalysisCut;

struct my_muon_o2cut {
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

  AnalysisCut* mMuonCut = nullptr;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
    const AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T}"};
    const AxisSpec axisPDca{nBinsPDca, minPDca, maxPDca, "PDca muon"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta}, true);
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt}, true);
    histos.add("etaPtHistogram", "#eta vs pT", kTH2F, {axisPt, axisEta});
    histos.add("PDcamuonHistogram", "PDca muon Histogram", kTH1F, {axisPDca}, true);

    mMuonCut = o2::aod::dqcuts::GetAnalysisCut("matchedQualityCuts"); 
    
  }

// have to change to fwd tracks
// Tutorial TracksIU -> FwdTrack
  void process(aod::FwdTracks const& fwdtracks)
  {
    for (auto& fwdtrack : fwdtracks) {
      VarManager::FillTrack<VarManager::Muon>(fwdtrack, VarManager::fgValues); 
      
      // 2. mMuonCut->IsSelected に VarManager を渡してカット判定
      if (!mMuonCut->IsSelected(VarManager::fgValues)) {
          continue; 
      }
      histos.fill(HIST("etaHistogram"), fwdtrack.eta());
      histos.fill(HIST("ptHistogram"), fwdtrack.pt());
      histos.fill(HIST("etaPtHistogram"), fwdtrack.pt(), fwdtrack.eta());
      histos.fill(HIST("PDcamuonHistogram"), fwdtrack.pDca());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<my_muon_o2cut>(cfgc)};
}
