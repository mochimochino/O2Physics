#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <cmath>
#include <TMath.h>
#include <TVector2.h>

using namespace o2;
using namespace o2::framework;

struct MchMftResiduals {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> maxphires{"maxphires", 1.0, "Max phi residual [rad]"};
  Configurable<float> minphires{"minphires", -1.0, "Min phi residual [rad]"};
  Configurable<float> maxetares{"maxetares", 1.0, "Max eta residual"};
  Configurable<float> minetares{"minetares", -1.0, "Min eta residual"};

  AxisSpec axisPt{100, 0.0, 10.0, "p_{T} (MCH) [GeV/c]"};
  AxisSpec axisResEta{200, minetares, maxetares, "#Delta#eta (#eta_{MCH} - #eta_{MFT})"};
  AxisSpec axisResPhi{200, minphires, maxphires, "#Delta#phi (#phi_{MCH} - #phi_{MFT}) [rad]"};

  void init(InitContext const&)
  {
    histos.add("True/hResEta", "Residual #eta (True);#Delta#eta;Counts", kTH1F, {axisResEta});
    histos.add("True/hResPhi", "Residual #phi (True);#Delta#phi;Counts", kTH1F, {axisResPhi});
    histos.add("True/hResEtaPhi", "Residual 2D (True);#Delta#eta;#Delta#phi", kTH2F, {axisResEta, axisResPhi});
    
    histos.add("Fake/hResEta", "Residual #eta (Fake);#Delta#eta;Counts", kTH1F, {axisResEta});
    histos.add("Fake/hResPhi", "Residual #phi (Fake);#Delta#phi;Counts", kTH1F, {axisResPhi});
    histos.add("Fake/hResEtaPhi", "Residual 2D (Fake);#Delta#eta;#Delta#phi", kTH2F, {axisResEta, axisResPhi});
    histos.add("Fake/hPt_vs_ResPhi", "p_{T} vs #Delta#phi (Fake);p_{T} [GeV/c];#Delta#phi", kTH2F, {axisPt, axisResPhi});
  }

void process(aod::FwdTracks const& mchTracks,
               aod::MFTTracks const& mftTracks,
               aod::McFwdTrackLabels const& mchLabels,
               aod::McMFTTrackLabels const& mftLabels)
  {
    for (auto const& mchTr : mchTracks) {
      
      int iMFT = mchTr.matchMFTTrackId();

      if (iMFT < 0) {
        continue;
      }

      auto const& mftTr = mftTracks.iteratorAt(iMFT);


      int mcIdMCH = -1;
      int mcIdMFT = -1;

      if (mchTr.globalIndex() < mchLabels.size()) {
          auto const& lbl = mchLabels.iteratorAt(mchTr.globalIndex());
          mcIdMCH = lbl.mcParticleId();
      }

      if (iMFT < mftLabels.size()) {
          auto const& lbl = mftLabels.iteratorAt(iMFT);
          mcIdMFT = lbl.mcParticleId();
      }

      bool isTrue = (mcIdMCH >= 0) && (mcIdMFT >= 0) && (mcIdMCH == mcIdMFT);


      float etaMCH = mchTr.eta();
      float phiMCH = mchTr.phi();
      float ptMCH  = mchTr.pt();

      float etaMFT = mftTr.eta();
      float phiMFT = mftTr.phi();

      float dEta = etaMCH - etaMFT;
      float dPhi = TVector2::Phi_mpi_pi(phiMCH - phiMFT);


      if (isTrue) {
          histos.fill(HIST("True/hResEta"), dEta);
          histos.fill(HIST("True/hResPhi"), dPhi);
          histos.fill(HIST("True/hResEtaPhi"), dEta, dPhi);
      } else {
          histos.fill(HIST("Fake/hResEta"), dEta);
          histos.fill(HIST("Fake/hResPhi"), dPhi);
          histos.fill(HIST("Fake/hResEtaPhi"), dEta, dPhi);
          histos.fill(HIST("Fake/hPt_vs_ResPhi"), ptMCH, dPhi);
      }
    }
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<MchMftResiduals>(cfg)
  };
}