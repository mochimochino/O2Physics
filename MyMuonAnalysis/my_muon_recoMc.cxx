#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"


using namespace o2;
using namespace o2::framework;

struct mc_muon_recoMc {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    AxisSpec axisPt{200, 0.0, 6.0, "p_{T}^{gen} [GeV/c]"};
    AxisSpec axisRes{200, -0.5, 0.5, "(pT_{reco} - pT_{gen}) / pT_{gen}"};
    AxisSpec axisEtaRes{200, -0.25, 0.25, "#eta_{reco} - #eta_{gen}"};
    AxisSpec axisPhiRes{200, -0.25, 0.25, "#phi_{reco} - #phi_{gen}"};


    histos.add("pt_reco_vs_gen", "pT_reco vs pT_gen", kTH2F, {axisPt, axisPt});
    histos.add("pt_resolution", "(pT_reco - pT_gen)/pT_gen", kTH1F, {axisRes});
    histos.add("eta_resolution", "eta_reco - eta_gen", kTH1F, {axisEtaRes});
    histos.add("phi_resolution", "phi_reco - phi_gen", kTH1F, {axisPhiRes});

    // resolution vs pT
    histos.add("ptRes_vs_pt", "Resolution vs pT_gen", kTH2F, {axisPt, axisRes});
    histos.add("etaRes_vs_pt", "Eta Resolution vs pT_gen", kTH2F, {axisPt, axisEtaRes});
    histos.add("phiRes_vs_pt", "Phi Resolution vs pT_gen", kTH2F, {axisPt, axisPhiRes});
  }

  // Join FwdTracks with McFwdTrackLabels
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

      float resPt = (ptGen - ptReco) / ptGen;
      float resEta = etaGen - etaReco;
      float resPhi = phiGen - phiReco;


      histos.fill(HIST("pt_reco_vs_gen"), ptReco, ptGen);
      histos.fill(HIST("pt_resolution"), resPt);
      histos.fill(HIST("eta_resolution"), resEta);
      histos.fill(HIST("ptRes_vs_pt"), ptGen, resPt);
      // misstta
      histos.fill(HIST("etaRes_vs_pt"), ptGen, resEta);
      histos.fill(HIST("phiRes_vs_pt"), ptGen, resPhi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<mc_muon_recoMc>(cfg)
  };
}
