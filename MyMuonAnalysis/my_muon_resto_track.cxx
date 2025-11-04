#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

struct mc_muon_resto_track {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable parameters
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0.0, "min pT"};
  Configurable<float> maxPt{"maxPt", 6.0, "max pT"};
  Configurable<int> nBinsRes{"nBinsRes", 200, "N bins in resolution histo"};
  Configurable<float> minRes{"minRes", -0.5, "min resolution"};
  Configurable<float> maxRes{"maxRes", 0.5, "max resolution"};

  const std::array<int, 4> trackTypes = {0, 2, 3, 4};

  void init(InitContext const&)
  {
    AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T}^{gen} [GeV/c]"};
    AxisSpec axisRes{nBinsRes, minRes, maxRes, "(p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen}"};
    AxisSpec axisEtaRes{200, -0.25, 0.25, "#eta_{reco} - #eta_{gen}"};
    AxisSpec axisPhiRes{200, -0.25, 0.25, "#phi_{reco} - #phi_{gen}"};

    histos.add("pt_reco_vs_gen", "pT_reco vs pT_gen", kTH2F, {axisPt, axisPt});
    histos.add("pt_resolution", "(pT_reco - pT_gen)/pT_gen", kTH1F, {axisRes});
    histos.add("eta_resolution", "eta_reco - eta_gen", kTH1F, {axisEtaRes});
    histos.add("phi_resolution", "phi_reco - phi_gen", kTH1F, {axisPhiRes});
    histos.add("ptRes_vs_pt", "Resolution vs pT_gen", kTH2F, {axisPt, axisRes});
    histos.add("etaRes_vs_pt", "Eta Resolution vs pT_gen", kTH2F, {axisPt, axisEtaRes});
    histos.add("phiRes_vs_pt", "Phi Resolution vs pT_gen", kTH2F, {axisPt, axisPhiRes});

    for (int type : trackTypes) {
      std::string t = std::to_string(type);
      histos.add(("pt_resolution_Type" + t).c_str(), ("pT resolution (Type " + t + ")").c_str(), kTH1F, {axisRes});
      histos.add(("eta_resolution_Type" + t).c_str(), ("eta resolution (Type " + t + ")").c_str(), kTH1F, {axisEtaRes});
      histos.add(("phi_resolution_Type" + t).c_str(), ("phi resolution (Type " + t + ")").c_str(), kTH1F, {axisPhiRes});
      histos.add(("ptRes_vs_pt_Type" + t).c_str(), ("pT res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisRes});
      histos.add(("etaRes_vs_pt_Type" + t).c_str(), ("eta res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisEtaRes});
      histos.add(("phiRes_vs_pt_Type" + t).c_str(), ("phi res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisPhiRes});
    }
  }

  void process(soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks,
               aod::McParticles const& mcParticles)
  {
    for (auto& trk : tracks) {
      int mcId = trk.mcParticleId();
      if (mcId < 0 || mcId >= static_cast<int>(mcParticles.size())) continue;

      auto mc = mcParticles.iteratorAt(mcId);
      if (std::abs(mc.pdgCode()) != 13) continue; // muon only

      float ptReco = trk.pt();
      float ptGen = mc.pt();
      if (ptGen <= 0) continue;

      float etaReco = trk.eta();
      float etaGen = mc.eta();
      float phiReco = trk.phi();
      float phiGen = mc.phi();

      float resPt = (ptReco - ptGen) / ptGen;
      float resEta = etaReco - etaGen;
      float resPhi = phiReco - phiGen;

      int type = trk.trackType();

      histos.fill(HIST("pt_reco_vs_gen"), ptReco, ptGen);
      histos.fill(HIST("pt_resolution"), resPt);
      histos.fill(HIST("eta_resolution"), resEta);
      histos.fill(HIST("phi_resolution"), resPhi);
      histos.fill(HIST("ptRes_vs_pt"), ptGen, resPt);
      histos.fill(HIST("etaRes_vs_pt"), ptGen, resEta);
      histos.fill(HIST("phiRes_vs_pt"), ptGen, resPhi);

      switch (type) {
        case 0: //global MUON (MFT+MCH+MID)
          histos.fill(HIST("pt_resolution_Type0"), resPt);
          histos.fill(HIST("eta_resolution_Type0"), resEta);
          histos.fill(HIST("phi_resolution_Type0"), resPhi);
          histos.fill(HIST("ptRes_vs_pt_Type0"), ptGen, resPt);
          histos.fill(HIST("etaRes_vs_pt_Type0"), ptGen, resEta);
          histos.fill(HIST("phiRes_vs_pt_Type0"), ptGen, resPhi);
          break;

        case 2: // MFT+MCH track
          histos.fill(HIST("pt_resolution_Type2"), resPt);
          histos.fill(HIST("eta_resolution_Type2"), resEta);
          histos.fill(HIST("phi_resolution_Type2"), resPhi);
          histos.fill(HIST("ptRes_vs_pt_Type2"), ptGen, resPt);
          histos.fill(HIST("etaRes_vs_pt_Type2"), ptGen, resEta);
          histos.fill(HIST("phiRes_vs_pt_Type2"), ptGen, resPhi);
          break;

        case 3: //standalone MUON (MCH+MID)
          histos.fill(HIST("pt_resolution_Type3"), resPt);
          histos.fill(HIST("eta_resolution_Type3"), resEta);
          histos.fill(HIST("phi_resolution_Type3"), resPhi);
          histos.fill(HIST("ptRes_vs_pt_Type3"), ptGen, resPt);
          histos.fill(HIST("etaRes_vs_pt_Type3"), ptGen, resEta);
          histos.fill(HIST("phiRes_vs_pt_Type3"), ptGen, resPhi);
          break;

        case 4: //standalone MCH track
          histos.fill(HIST("pt_resolution_Type4"), resPt);
          histos.fill(HIST("eta_resolution_Type4"), resEta);
          histos.fill(HIST("phi_resolution_Type4"), resPhi);
          histos.fill(HIST("ptRes_vs_pt_Type4"), ptGen, resPt);
          histos.fill(HIST("etaRes_vs_pt_Type4"), ptGen, resEta);
          histos.fill(HIST("phiRes_vs_pt_Type4"), ptGen, resPhi);
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
    adaptAnalysisTask<mc_muon_resto_track>(cfg)};
}
