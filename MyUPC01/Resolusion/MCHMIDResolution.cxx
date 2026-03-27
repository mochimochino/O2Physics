// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file \"COPYING\".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \brief  pT / eta / phi resolution of MCH-MID muon tracks using MC truth.
///         Step 1 of the binning-strategy study:
///           - Fill (reco - MC) residuals as a function of MC-truth pT, eta, phi
///           - Fill 2D (pT_reco vs pT_MC) for migration-matrix inspection
///         Only "true" matched tracks (mcMask == 0) with standard acceptance cuts
///         are used so that fake matches do not bias the resolution estimate.
/// \author Takuma
/// \date   2026
///

#include <cmath>
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "PWGDQ/Core/VarManager.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// -----------------------------------------------------------------------
// Table joins
// MCH-MID tracks (type 3) with covariance and MC labels
using MCHMIDTracks = soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov, o2::aod::McFwdTrackLabels>;
using MyEvents     = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using ExtBCs       = soa::Join<aod::BCs, aod::Timestamps>;
// -----------------------------------------------------------------------

struct MCHMIDResolutionTask {

  // -----------------------------------------------------------------------
  // CCDB
  // -----------------------------------------------------------------------
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::parameters::GRPMagField*        grpmag     = nullptr;
  int                                 mCurrentRun = -1;

  Configurable<std::string> ccdburl    {"ccdb-url",    "http://alice-ccdb.cern.ch",    "URL of the CCDB repository"};
  Configurable<std::string> grpmagPath {"grpmagPath",  "GLO/Config/GRPMagField",       "CCDB path of GRPMagField"};
  Configurable<std::string> geoPath    {"geoPath",
                                         "GLO/Config/GeometryAligned",
                                         "Path of the aligned geometry"};

  // -----------------------------------------------------------------------
  // Acceptance configurable cuts  (same defaults as the reference task)
  // -----------------------------------------------------------------------
  Configurable<float> cfgEtaMin  {"cfgEtaMin",  -3.6f, "Minimum eta for MCH-MID tracks"};
  Configurable<float> cfgEtaMax  {"cfgEtaMax",  -2.5f, "Maximum eta for MCH-MID tracks"};
  Configurable<float> cfgRabsMin {"cfgRabsMin",  17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> cfgRabsMax {"cfgRabsMax",  89.5f, "Maximum R at absorber end [cm]"};
  // pDCA cut thresholds
  Configurable<float> cfgPDCAcut1z {"cfgPDCAcut1z", 26.5f, "R boundary for the two pDCA slopes [cm]"};
  Configurable<float> cfgPDCAcut1  {"cfgPDCAcut1",  594.0f, "pDCA cut for R < cfgPDCAcut1z"};
  Configurable<float> cfgPDCAcut2  {"cfgPDCAcut2",  324.0f, "pDCA cut for R >= cfgPDCAcut1z"};

  // -----------------------------------------------------------------------
  // Histogram registry
  // -----------------------------------------------------------------------
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ---- axis definitions -------------------------------------------------
  // MC-truth kinematics
  AxisSpec axisPtMC  {120,  0.0,  12.0, "p_{T}^{MC} [GeV/c]"};
  AxisSpec axisEtaMC {100, -4.5,  -2.0, "#eta^{MC}"};
  AxisSpec axisPhiMC { 64, -M_PI, M_PI, "#phi^{MC} [rad]"};

  // Reco kinematics (for 2D migration matrices)
  AxisSpec axisPtReco  {120,  0.0,  12.0, "p_{T}^{reco} [GeV/c]"};
  AxisSpec axisEtaReco {100, -4.5,  -2.0, "#eta^{reco}"};
  AxisSpec axisPhiReco { 64, -M_PI, M_PI, "#phi^{reco} [rad]"};

  // Resolution axes  (reco - MC) / MC  or  (reco - MC)
  // pT relative resolution: (pT_reco - pT_MC) / pT_MC
  AxisSpec axisDeltaPtRel {200, -0.5,   0.5,  "(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}"};
  // pT absolute difference [GeV/c]
  AxisSpec axisDeltaPtAbs {200, -2.0,   2.0,  "p_{T}^{reco} - p_{T}^{MC} [GeV/c]"};
  // eta difference
  AxisSpec axisDeltaEta   {200, -0.2,   0.2,  "#eta^{reco} - #eta^{MC}"};
  // phi difference
  AxisSpec axisDeltaPhi   {200, -0.2,   0.2,  "#phi^{reco} - #phi^{MC} [rad]"};
  // 1/pT relative
  AxisSpec axisDeltaInvPt {200, -0.5,   0.5,  "(1/p_{T}^{reco} - 1/p_{T}^{MC}) / (1/p_{T}^{MC})"};

  // charge
  AxisSpec axisCharge {3, -1.5, 1.5, "charge q"};

  // -----------------------------------------------------------------------
  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    VarManager::SetDefaultVarNames();

    // ---- 1D resolution histograms (projected onto MC pT slices) ----------
    //   DeltaPtRel  vs  pT_MC
    histos.add("Resolution/hDeltaPtRel_vs_PtMC",
               "#Deltap_{T}/p_{T}^{MC} vs p_{T}^{MC}; p_{T}^{MC} [GeV/c]; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisPtMC, axisDeltaPtRel});

    //   DeltaPtAbs  vs  pT_MC
    histos.add("Resolution/hDeltaPtAbs_vs_PtMC",
               "#Deltap_{T} (abs) vs p_{T}^{MC}; p_{T}^{MC} [GeV/c]; p_{T}^{reco}-p_{T}^{MC} [GeV/c]",
               kTH2F, {axisPtMC, axisDeltaPtAbs});

    //   Delta(1/pT) vs pT_MC  (most Gaussian due to tracking)
    histos.add("Resolution/hDeltaInvPtRel_vs_PtMC",
               "#Delta(1/p_{T})/(1/p_{T}^{MC}) vs p_{T}^{MC}; p_{T}^{MC} [GeV/c]; #Delta(1/p_{T})/(1/p_{T}^{MC})",
               kTH2F, {axisPtMC, axisDeltaInvPt});

    //   Delta eta vs eta_MC
    histos.add("Resolution/hDeltaEta_vs_EtaMC",
               "#Delta#eta vs #eta^{MC}; #eta^{MC}; #eta^{reco}-#eta^{MC}",
               kTH2F, {axisEtaMC, axisDeltaEta});

    //   Delta phi vs phi_MC
    histos.add("Resolution/hDeltaPhi_vs_PhiMC",
               "#Delta#phi vs #phi^{MC}; #phi^{MC} [rad]; #phi^{reco}-#phi^{MC} [rad]",
               kTH2F, {axisPhiMC, axisDeltaPhi});

    // ---- Resolution vs pT_MC (also useful: eta slice) --------------------
    histos.add("Resolution/hDeltaPtRel_vs_EtaMC",
               "#Deltap_{T}/p_{T}^{MC} vs #eta^{MC}; #eta^{MC}; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisEtaMC, axisDeltaPtRel});

    // ---- 2D migration matrix (pT_reco vs pT_MC) --------------------------
    histos.add("Migration/hPtReco_vs_PtMC",
               "p_{T} migration matrix; p_{T}^{MC} [GeV/c]; p_{T}^{reco} [GeV/c]",
               kTH2F, {axisPtMC, axisPtReco});

    histos.add("Migration/hEtaReco_vs_EtaMC",
               "#eta migration matrix; #eta^{MC}; #eta^{reco}",
               kTH2F, {axisEtaMC, axisEtaReco});

    histos.add("Migration/hPhiReco_vs_PhiMC",
               "#phi migration matrix; #phi^{MC} [rad]; #phi^{reco} [rad]",
               kTH2F, {axisPhiMC, axisPhiReco});

    // ---- 1D projections of reco and MC kinematics (sanity check) ---------
    histos.add("Kinematics/hPtMC",   "MC truth p_{T}; p_{T}^{MC} [GeV/c]; Counts", kTH1F, {axisPtMC});
    histos.add("Kinematics/hPtReco", "Reco p_{T};     p_{T}^{reco} [GeV/c]; Counts", kTH1F, {axisPtReco});
    histos.add("Kinematics/hEtaMC",  "MC truth #eta;  #eta^{MC}; Counts",   kTH1F, {axisEtaMC});
    histos.add("Kinematics/hEtaReco","Reco #eta;      #eta^{reco}; Counts",  kTH1F, {axisEtaReco});
    histos.add("Kinematics/hPhiMC",  "MC truth #phi;  #phi^{MC} [rad]; Counts",  kTH1F, {axisPhiMC});
    histos.add("Kinematics/hPhiReco","Reco #phi;      #phi^{reco} [rad]; Counts", kTH1F, {axisPhiReco});

    // ---- 3D: DeltaPtRel vs (pT_MC, eta_MC)  -----------------------------
    // Useful to see acceptance-dependent resolution
    histos.add("Resolution/hDeltaPtRel_PtMC_EtaMC",
               "#Deltap_{T}/p_{T}^{MC} vs (p_{T}^{MC}, #eta^{MC}); p_{T}^{MC} [GeV/c]; #eta^{MC}; #Deltap_{T}/p_{T}^{MC}",
               kTH3F, {axisPtMC, axisEtaMC, axisDeltaPtRel});

    // ---- charge split (mu+ / mu-)  ---------------------------------------
    histos.add("Resolution/hDeltaPtRel_vs_PtMC_Pos",
               "#Deltap_{T}/p_{T}^{MC} vs p_{T}^{MC} (#mu^{+}); p_{T}^{MC} [GeV/c]; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisPtMC, axisDeltaPtRel});
    histos.add("Resolution/hDeltaPtRel_vs_PtMC_Neg",
               "#Deltap_{T}/p_{T}^{MC} vs p_{T}^{MC} (#mu^{-}); p_{T}^{MC} [GeV/c]; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisPtMC, axisDeltaPtRel});

    // ---- QA counter -------------------------------------------------------
    histos.add("QA/hTrackTypeAll",  "Track type (before cuts); Track type; Counts", kTH1I,
               {{5, -0.5, 4.5, "track type"}});
    histos.add("QA/hTrackTypePass", "Track type (after cuts);  Track type; Counts", kTH1I,
               {{5, -0.5, 4.5, "track type"}});
    histos.add("QA/hMcMaskPass",    "mcMask of selected tracks; mcMask; Counts",    kTH1I,
               {{20, -0.5, 19.5, "mcMask"}});
  }

  // -----------------------------------------------------------------------
  // Helper: standard acceptance cuts for MCH-MID tracks (type == 3)
  // Returns true if the track passes all cuts.
  // -----------------------------------------------------------------------
  template <typename TTrack>
  bool acceptMCHMID(TTrack const& tr) const
  {
    if (tr.trackType() != 3) return false;                // MCH-MID only
    if (tr.eta() < cfgEtaMin || tr.eta() > cfgEtaMax) return false;
    const float rAbs = tr.rAtAbsorberEnd();
    if (rAbs < cfgRabsMin || rAbs > cfgRabsMax) return false;
    const float pDca = tr.pDca();
    if (pDca < 0.f) return false;
    if (rAbs < cfgPDCAcut1z) {
      if (pDca > cfgPDCAcut1) return false;
    } else {
      if (pDca > cfgPDCAcut2) return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------
  // Main process function
  // -----------------------------------------------------------------------
  void process(MCHMIDTracks const& tracks,
               MyEvents       const& collisions,
               ExtBCs         const& bcs,
               aod::McParticles const& mcParticles)
  {
    // ---- CCDB / magnetic field update (once per run) --------------------
    if (bcs.size() > 0) {
      const int runNumber = bcs.begin().runNumber();
      if (runNumber != mCurrentRun) {
        const long ts = bcs.begin().timestamp();
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, ts);
        if (grpmag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmag);
          VarManager::SetMagneticField(grpmag->getNominalL3Field());
          VarManager::SetMatchingPlane(-77.5f);
        }
        VarManager::SetupMuonMagField();
        mCurrentRun = runNumber;
        LOG(info) << "Run " << mCurrentRun << ": magnetic field updated.";
      }
    }

    // ---- Loop over tracks -----------------------------------------------
    for (auto const& tr : tracks) {

      histos.fill(HIST("QA/hTrackTypeAll"), tr.trackType());

      // --- acceptance selection ---
      if (!acceptMCHMID(tr)) continue;

      histos.fill(HIST("QA/hTrackTypePass"), tr.trackType());
      histos.fill(HIST("QA/hMcMaskPass"),    tr.mcMask());

      // --- require true match (no mismatched hits) -----------------------
      if (tr.mcMask() != 0) continue;

      // --- retrieve MC particle ------------------------------------------
      const int mcId = tr.mcParticleId();
      if (mcId < 0) continue;

      auto const& mc = mcParticles.iteratorAt(mcId);

      // MC kinematics
      const double ptMC  = mc.pt();
      const double etaMC = mc.eta();
      const double phiMC = mc.phi();

      // Reco kinematics (raw values stored in FwdTracks)
      // For MCH-MID (type 3) no MFT propagation is applied here;
      // we compare what the tracker reconstructed versus the generator.
      const double ptReco  = tr.pt();
      const double etaReco = tr.eta();
      const double phiReco = tr.phi();

      // Relative pT resolution
      const double deltaPtRel   = (ptMC > 0.) ? (ptReco - ptMC) / ptMC   : -999.;
      const double deltaPtAbs   = ptReco - ptMC;
      const double deltaEta     = etaReco - etaMC;

      // phi difference mapped to [-pi, pi]
      double deltaPhi = phiReco - phiMC;
      while (deltaPhi >  M_PI) deltaPhi -= 2. * M_PI;
      while (deltaPhi < -M_PI) deltaPhi += 2. * M_PI;

      // Delta(1/pT) relative
      const double deltaInvPt = (ptMC > 0.) ? (1./ptReco - 1./ptMC) / (1./ptMC) : -999.;

      // Charge from reconstructed signed 1/pT
      const double invQPt  = tr.invQPt(); // signed q/pT
      const int    charge  = (invQPt > 0.) ? 1 : -1;

      // ---- Fill resolution histograms ----------------------------------
      histos.fill(HIST("Resolution/hDeltaPtRel_vs_PtMC"),     ptMC, deltaPtRel);
      histos.fill(HIST("Resolution/hDeltaPtAbs_vs_PtMC"),     ptMC, deltaPtAbs);
      histos.fill(HIST("Resolution/hDeltaInvPtRel_vs_PtMC"),  ptMC, deltaInvPt);
      histos.fill(HIST("Resolution/hDeltaEta_vs_EtaMC"),      etaMC, deltaEta);
      histos.fill(HIST("Resolution/hDeltaPhi_vs_PhiMC"),      phiMC, deltaPhi);
      histos.fill(HIST("Resolution/hDeltaPtRel_vs_EtaMC"),    etaMC, deltaPtRel);
      histos.fill(HIST("Resolution/hDeltaPtRel_PtMC_EtaMC"),  ptMC, etaMC, deltaPtRel);

      // charge-split pT resolution
      if (charge > 0) {
        histos.fill(HIST("Resolution/hDeltaPtRel_vs_PtMC_Pos"), ptMC, deltaPtRel);
      } else {
        histos.fill(HIST("Resolution/hDeltaPtRel_vs_PtMC_Neg"), ptMC, deltaPtRel);
      }

      // ---- Fill migration matrices -------------------------------------
      histos.fill(HIST("Migration/hPtReco_vs_PtMC"),   ptMC, ptReco);
      histos.fill(HIST("Migration/hEtaReco_vs_EtaMC"), etaMC, etaReco);
      histos.fill(HIST("Migration/hPhiReco_vs_PhiMC"), phiMC, phiReco);

      // ---- Sanity-check 1D distributions --------------------------------
      histos.fill(HIST("Kinematics/hPtMC"),    ptMC);
      histos.fill(HIST("Kinematics/hPtReco"),  ptReco);
      histos.fill(HIST("Kinematics/hEtaMC"),   etaMC);
      histos.fill(HIST("Kinematics/hEtaReco"), etaReco);
      histos.fill(HIST("Kinematics/hPhiMC"),   phiMC);
      histos.fill(HIST("Kinematics/hPhiReco"), phiReco);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<MCHMIDResolutionTask>(cfg)
  };
}
