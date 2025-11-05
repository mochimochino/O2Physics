#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <set>

using namespace o2;
using namespace o2::framework;

struct my_muon_tracking_efficiency {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable parameters pT Eta Phi
  Configurable<int> nBinsPt{"nBinsPt", 60, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0.0, "min pT"};
  Configurable<float> maxPt{"maxPt", 6.0, "max pT"};
  Configurable<int> nBinsEta{"nBinsEta", 50, "N bins in Eta histo"};
  Configurable<float> minEta{"minEta", -4.0, "min Eta"};
  Configurable<float> maxEta{"maxEta", -2.5, "max Eta"};
  Configurable<int> nBinsPhi{"nBinsPhi", 64, "N bins in Phi histo"};
  Configurable<float> minPhi{"minPhi", - TMath::Pi(), "min Phi"};
  Configurable<float> maxPhi{"maxPhi", TMath::Pi(), "max Phi"};

  const std::array<int, 4> trackTypes = {0, 2, 3, 4};

  void init(InitContext const&)
    {
        AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T} [GeV/c]"}; 
        AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
        AxisSpec axisPhi{nBinsPhi, minPhi, maxPhi, "#phi [rad]"};

        histos.add("hPt_Gen", "Generated pT (Denominator)", kTH1F, {axisPt});
        histos.add("hEta_Gen", "Generated eta (Denominator)", kTH1F, {axisEta});
        histos.add("hPhi_Gen", "Generated phi (Denominator)", kTH1F, {axisPhi});

        histos.add("hPt_Reco_Matched_All", "Matched Reco pT (Numerator, All Types)", kTH1F, {axisPt});
        histos.add("hEta_Reco_Matched_All", "Matched Reco eta (Numerator, All Types)", kTH1F, {axisEta});
        histos.add("hPhi_Reco_Matched_All", "Matched Reco phi (Numerator, All Types)", kTH1F, {axisPhi});

        for (int type : trackTypes) {
            std::string t = std::to_string(type);
            histos.add(("hPt_Reco_Matched_Type" + t).c_str(), ("Matched Reco pT (Numerator, Type " + t + ")").c_str(), kTH1F, {axisPt});
            histos.add(("hEta_Reco_Matched_Type" + t).c_str(), ("Matched Reco eta (Numerator, Type " + t + ")").c_str(), kTH1F, {axisEta});
            histos.add(("hPhi_Reco_Matched_Type" + t).c_str(), ("Matched Reco phi (Numerator, Type " + t + ")").c_str(), kTH1F, {axisPhi});
        }
    } 

  void process(soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks, aod::McParticles const& mcParticles)
   {
    auto normPhi = [](double a) {
      while (a > TMath::Pi()) a -= 2.0 * TMath::Pi();
      while (a <= -TMath::Pi()) a += 2.0 * TMath::Pi();
      return a;
    };


    for (size_t i = 0; i < mcParticles.size(); ++i) {
        auto mc = mcParticles.iteratorAt(i);

        if (std::abs(mc.pdgCode()) != 13) continue; 
        
        float ptGen = mc.pt();
        float etaGen = mc.eta();
        float phiGen = normPhi(mc.phi());

        if (ptGen < minPt || ptGen > maxPt) continue;
        if (etaGen < minEta || etaGen > maxEta) continue;
        // if (phiGen < minPhi || phiGen > maxPhi) continue;


        histos.fill(HIST("hPt_Gen"), ptGen);
        histos.fill(HIST("hEta_Gen"), etaGen);
        histos.fill(HIST("hPhi_Gen"), phiGen);
    }

    std::set<int> countedMcIds;
    for (auto& trk : tracks) {
        int mcId = trk.mcParticleId();
        if (mcId < 0 || mcId >= static_cast<int>(mcParticles.size())) continue;

        auto mc = mcParticles.iteratorAt(mcId);

        if (std::abs(mc.pdgCode()) != 13) continue;

        // 2. Acceptance cuts
        float ptGen = mc.pt();
        float etaGen = mc.eta();
        float phiGen = normPhi(mc.phi());

        if (ptGen < minPt || ptGen > maxPt) continue;
        if (etaGen < minEta || etaGen > maxEta) continue;
        // if (phiGen < minPhi || phiGen > maxPhi) continue;


        histos.fill(HIST("hPt_Reco_Matched_All"), ptGen);
        histos.fill(HIST("hEta_Reco_Matched_All"), etaGen);
        histos.fill(HIST("hPhi_Reco_Matched_All"), phiGen);

        // Fill per-type histograms
        int type = trk.trackType();
        switch (type) {
        case 0:
            histos.fill(HIST("hPt_Reco_Matched_Type0"), ptGen);
            histos.fill(HIST("hEta_Reco_Matched_Type0"), etaGen);
            histos.fill(HIST("hPhi_Reco_Matched_Type0"), phiGen);
            break;
        case 2:
            histos.fill(HIST("hPt_Reco_Matched_Type2"), ptGen);
            histos.fill(HIST("hEta_Reco_Matched_Type2"), etaGen);
            histos.fill(HIST("hPhi_Reco_Matched_Type2"), phiGen);
            break;
        case 3:
            histos.fill(HIST("hPt_Reco_Matched_Type3"), ptGen);
            histos.fill(HIST("hEta_Reco_Matched_Type3"), etaGen);
            histos.fill(HIST("hPhi_Reco_Matched_Type3"), phiGen);
            break;
        case 4:
            histos.fill(HIST("hPt_Reco_Matched_Type4"), ptGen);
            histos.fill(HIST("hEta_Reco_Matched_Type4"), etaGen);
            histos.fill(HIST("hPhi_Reco_Matched_Type4"), phiGen);
            break;
        default:
            break;
        }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<my_muon_tracking_efficiency>(cfg)};
}