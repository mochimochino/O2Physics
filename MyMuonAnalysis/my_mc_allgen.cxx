#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"


using namespace o2;
using namespace o2::framework;

struct mc_mc_allgen {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 30, "N bins in Eta histo"};
  Configurable<float> minEta{"minEta", -5.0, "min Eta in Eta histo"};
  Configurable<float> maxEta{"maxEta", -2.0, "max Eta in Eta histo"};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0, "min pT in pT histo"};
  Configurable<float> maxPt{"maxPt", 6, "max pT in pT histo"};

  void init(InitContext const&)
  {
    AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
    AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T} [GeV/c]"};

    // all muon
    histos.add("pt_mc_all", "MC muon pT (all)", kTH1F, {axisPt});
    histos.add("eta_mc_all", "MC muon eta (all)", kTH1F, {axisEta});
    // primary, decay-from-pion/kaon, other-secondaries
    histos.add("pt_mc_primary", "MC muon pT (primary)", kTH1F, {axisPt});
    histos.add("pt_mc_fromPiK", "MC muon pT (from pi/K)", kTH1F, {axisPt});
    histos.add("pt_mc_otherSec", "MC muon pT (other secondary)", kTH1F, {axisPt});

    histos.add("etaPt_mc_all", "MC muon eta vs pT (all)", kTH2F, {axisPt, axisEta});
  }

  void process(aod::McParticles const& mcParticles)
  {
    for (auto& mc : mcParticles) {
      if (std::abs(mc.pdgCode()) != 13) continue; // muon only

      float pt = mc.pt();
      float eta = mc.eta();

      histos.fill(HIST("pt_mc_all"), pt);
      histos.fill(HIST("eta_mc_all"), eta);
      histos.fill(HIST("etaPt_mc_all"), pt, eta);

      if (mc.isPhysicalPrimary() && mc.producedByGenerator()) {
        histos.fill(HIST("pt_mc_primary"), pt);
        continue;
      }

      // from π/K
      bool fromPiK = false;
      if (mc.has_mothers()) {
        for (auto mother : mc.mothers_as<aod::McParticles>()) {
          int mpdg = mother.pdgCode();
          if (std::abs(mpdg) == 211 || std::abs(mpdg) == 321) {
            fromPiK = true;
            break;
          }
        }
      }

      if (fromPiK) {
        histos.fill(HIST("pt_mc_fromPiK"), pt);
      } else {
        // pother secondary
        if (!mc.producedByGenerator()) {
          histos.fill(HIST("pt_mc_otherSec"), pt);
        } else {

        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<mc_mc_allgen>(cfg)
  };
}
