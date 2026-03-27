// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
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

  Configurable<std::string> ccdburl    {"ccdb-url",   "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<std::string> grpmagPath {"grpmagPath", "GLO/Config/GRPMagField",    "CCDB path of GRPMagField"};
  Configurable<std::string> geoPath    {"geoPath",    "GLO/Config/GeometryAligned","Path of the aligned geometry"};

  // -----------------------------------------------------------------------
  // Acceptance cuts
  // -----------------------------------------------------------------------
  Configurable<float> cfgEtaMin    {"cfgEtaMin",    -3.6f,  "Minimum eta for MCH-MID tracks"};
  Configurable<float> cfgEtaMax    {"cfgEtaMax",    -2.5f,  "Maximum eta for MCH-MID tracks"};
  Configurable<float> cfgRabsMin   {"cfgRabsMin",   17.6f,  "Minimum R at absorber end [cm]"};
  Configurable<float> cfgRabsMax   {"cfgRabsMax",   89.5f,  "Maximum R at absorber end [cm]"};
  Configurable<float> cfgPDCAcut1z {"cfgPDCAcut1z", 26.5f,  "R boundary for the two pDCA slopes [cm]"};
  Configurable<float> cfgPDCAcut1  {"cfgPDCAcut1",  594.0f, "pDCA cut for R < cfgPDCAcut1z"};
  Configurable<float> cfgPDCAcut2  {"cfgPDCAcut2",  324.0f, "pDCA cut for R >= cfgPDCAcut1z"};

  // -----------------------------------------------------------------------
  // Histogram binning Configurables
  // (defined here so AxisSpec can be built from them inside init())
  // -----------------------------------------------------------------------
  // MC-truth pT axis
  Configurable<int>   nBinsPtMC  {"nBinsPtMC",   120,  "N bins for MC-truth pT axis"};
  Configurable<float> minPtMC    {"minPtMC",      0.0f, "Min of MC-truth pT axis [GeV/c]"};
  Configurable<float> maxPtMC    {"maxPtMC",     12.0f, "Max of MC-truth pT axis [GeV/c]"};
  // eta axis (shared for MC and reco)
  Configurable<int>   nBinsEta   {"nBinsEta",    100,  "N bins for eta axis"};
  Configurable<float> minEta     {"minEta",      -4.5f, "Min of eta axis"};
  Configurable<float> maxEta     {"maxEta",      -2.0f, "Max of eta axis"};
  // phi axis (shared)
  Configurable<int>   nBinsPhi   {"nBinsPhi",     64,  "N bins for phi axis"};
  // Resolution axes
  Configurable<int>   nBinsResRel{"nBinsResRel", 200,  "N bins for relative resolution axis"};
  Configurable<float> minResRel  {"minResRel",   -0.5f, "Min of relative resolution axis"};
  Configurable<float> maxResRel  {"maxResRel",    0.5f, "Max of relative resolution axis"};
  Configurable<int>   nBinsResAbs{"nBinsResAbs", 200,  "N bins for absolute pT resolution axis"};
  Configurable<float> minResAbs  {"minResAbs",   -2.0f, "Min of absolute pT resolution axis [GeV/c]"};
  Configurable<float> maxResAbs  {"maxResAbs",    2.0f, "Max of absolute pT resolution axis [GeV/c]"};
  Configurable<int>   nBinsResAngle{"nBinsResAngle", 200, "N bins for eta/phi residual axis"};
  Configurable<float> minResAngle{"minResAngle", -0.2f, "Min of eta/phi residual axis"};
  Configurable<float> maxResAngle{"maxResAngle",  0.2f, "Max of eta/phi residual axis"};

  // -----------------------------------------------------------------------
  // Histogram registry
  // -----------------------------------------------------------------------
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // -----------------------------------------------------------------------
  void init(InitContext const&)
  {
    // ---- CCDB setup -------------------------------------------------------
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    VarManager::SetDefaultVarNames();

    // ---- Build axes from Configurables (same style as MyUPCMass02) --------
    const AxisSpec axisPtMC    {nBinsPtMC,    minPtMC,    maxPtMC,    "p_{T}^{MC} [GeV/c]"};
    const AxisSpec axisPtReco  {nBinsPtMC,    minPtMC,    maxPtMC,    "p_{T}^{reco} [GeV/c]"};
    const AxisSpec axisEtaMC   {nBinsEta,     minEta,     maxEta,     "#eta^{MC}"};
    const AxisSpec axisEtaReco {nBinsEta,     minEta,     maxEta,     "#eta^{reco}"};
    const AxisSpec axisPhiMC   {nBinsPhi,    -M_PI,       M_PI,       "#phi^{MC} [rad]"};
    const AxisSpec axisPhiReco {nBinsPhi,    -M_PI,       M_PI,       "#phi^{reco} [rad]"};

    const AxisSpec axisDeltaPtRel  {nBinsResRel,   minResRel,   maxResRel,   "(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}"};
    const AxisSpec axisDeltaPtAbs  {nBinsResAbs,   minResAbs,   maxResAbs,   "p_{T}^{reco} - p_{T}^{MC} [GeV/c]"};
    const AxisSpec axisDeltaInvPt  {nBinsResRel,   minResRel,   maxResRel,   "(1/p_{T}^{reco} - 1/p_{T}^{MC}) / (1/p_{T}^{MC})"};
    const AxisSpec axisDeltaEta    {nBinsResAngle, minResAngle, maxResAngle, "#eta^{reco} - #eta^{MC}"};
    const AxisSpec axisDeltaPhi    {nBinsResAngle, minResAngle, maxResAngle, "#phi^{reco} - #phi^{MC} [rad]"};

    // ---- Resolution 2D histograms -----------------------------------------
    histos.add("Resolution/hDeltaPtRel_vs_PtMC",
               "#Deltap_{T}/p_{T}^{MC} vs p_{T}^{MC}; p_{T}^{MC} [GeV/c]; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisPtMC, axisDeltaPtRel});

    histos.add("Resolution/hDeltaPtAbs_vs_PtMC",
               "#Deltap_{T} (abs) vs p_{T}^{MC}; p_{T}^{MC} [GeV/c]; p_{T}^{reco}-p_{T}^{MC} [GeV/c]",
               kTH2F, {axisPtMC, axisDeltaPtAbs});

    histos.add("Resolution/hDeltaInvPtRel_vs_PtMC",
               "#Delta(1/p_{T})/(1/p_{T}^{MC}) vs p_{T}^{MC}; p_{T}^{MC} [GeV/c]; #Delta(1/p_{T})/(1/p_{T}^{MC})",
               kTH2F, {axisPtMC, axisDeltaInvPt});

    histos.add("Resolution/hDeltaEta_vs_EtaMC",
               "#Delta#eta vs #eta^{MC}; #eta^{MC}; #eta^{reco}-#eta^{MC}",
               kTH2F, {axisEtaMC, axisDeltaEta});

    histos.add("Resolution/hDeltaPhi_vs_PhiMC",
               "#Delta#phi vs #phi^{MC}; #phi^{MC} [rad]; #phi^{reco}-#phi^{MC} [rad]",
               kTH2F, {axisPhiMC, axisDeltaPhi});

    histos.add("Resolution/hDeltaPtRel_vs_EtaMC",
               "#Deltap_{T}/p_{T}^{MC} vs #eta^{MC}; #eta^{MC}; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisEtaMC, axisDeltaPtRel});

    // ---- 3D: resolution vs (pT_MC, eta_MC) --------------------------------
    histos.add("Resolution/hDeltaPtRel_PtMC_EtaMC",
               "#Deltap_{T}/p_{T}^{MC} vs (p_{T}^{MC}, #eta^{MC}); p_{T}^{MC} [GeV/c]; #eta^{MC}; #Deltap_{T}/p_{T}^{MC}",
               kTH3F, {axisPtMC, axisEtaMC, axisDeltaPtRel});

    // ---- Charge-split pT resolution ---------------------------------------
    histos.add("Resolution/hDeltaPtRel_vs_PtMC_Pos",
               "#Deltap_{T}/p_{T}^{MC} vs p_{T}^{MC} (#mu^{+}); p_{T}^{MC} [GeV/c]; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisPtMC, axisDeltaPtRel});
    histos.add("Resolution/hDeltaPtRel_vs_PtMC_Neg",
               "#Deltap_{T}/p_{T}^{MC} vs p_{T}^{MC} (#mu^{-}); p_{T}^{MC} [GeV/c]; (p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}",
               kTH2F, {axisPtMC, axisDeltaPtRel});

    // ---- Migration matrices -----------------------------------------------
    histos.add("Migration/hPtReco_vs_PtMC",
               "p_{T} migration matrix; p_{T}^{MC} [GeV/c]; p_{T}^{reco} [GeV/c]",
               kTH2F, {axisPtMC, axisPtReco});
    histos.add("Migration/hEtaReco_vs_EtaMC",
               "#eta migration matrix; #eta^{MC}; #eta^{reco}",
               kTH2F, {axisEtaMC, axisEtaReco});
    histos.add("Migration/hPhiReco_vs_PhiMC",
               "#phi migration matrix; #phi^{MC} [rad]; #phi^{reco} [rad]",
               kTH2F, {axisPhiMC, axisPhiReco});

    // ---- Kinematics sanity check ------------------------------------------
    histos.add("Kinematics/hPtMC",    "MC truth p_{T}; p_{T}^{MC} [GeV/c]; Counts",    kTH1F, {axisPtMC});
    histos.add("Kinematics/hPtReco",  "Reco p_{T};     p_{T}^{reco} [GeV/c]; Counts",  kTH1F, {axisPtReco});
    histos.add("Kinematics/hEtaMC",   "MC truth #eta;  #eta^{MC}; Counts",              kTH1F, {axisEtaMC});
    histos.add("Kinematics/hEtaReco", "Reco #eta;      #eta^{reco}; Counts",            kTH1F, {axisEtaReco});
    histos.add("Kinematics/hPhiMC",   "MC truth #phi;  #phi^{MC} [rad]; Counts",        kTH1F, {axisPhiMC});
    histos.add("Kinematics/hPhiReco", "Reco #phi;      #phi^{reco} [rad]; Counts",      kTH1F, {axisPhiReco});

    // ---- QA ---------------------------------------------------------------
    histos.add("QA/hTrackTypeAll",  "Track type (before cuts); Track type; Counts", kTH1I,
               {{5, -0.5, 4.5, "track type"}});
    histos.add("QA/hTrackTypePass", "Track type (after cuts);  Track type; Counts", kTH1I,
               {{5, -0.5, 4.5, "track type"}});
    histos.add("QA/hMcMaskPass",    "mcMask of selected tracks; mcMask; Counts",    kTH1I,
               {{20, -0.5, 19.5, "mcMask"}});
  }

  // -----------------------------------------------------------------------
  // Acceptance helper: MCH-MID standard cuts (type == 3)
  // -----------------------------------------------------------------------
  template <typename TTrack>
  bool acceptMCHMID(TTrack const& tr) const
  {
    if (tr.trackType() != 3) return false;
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
  // Main process
  // -----------------------------------------------------------------------
  void process(MCHMIDTracks    const& tracks,
               MyEvents        const& /*collisions*/,
               ExtBCs          const& bcs,
               aod::McParticles const& mcParticles)
  {
    // ---- Magnetic field update (once per run) ----------------------------
    if (bcs.size() > 0) {
      const int runNumber = bcs.begin().runNumber();
      if (runNumber != mCurrentRun) {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(
                   grpmagPath, bcs.begin().timestamp());
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

    // ---- Track loop -------------------------------------------------------
    for (auto const& tr : tracks) {

      histos.fill(HIST("QA/hTrackTypeAll"), tr.trackType());

      if (!acceptMCHMID(tr)) continue;

      histos.fill(HIST("QA/hTrackTypePass"), tr.trackType());
      histos.fill(HIST("QA/hMcMaskPass"),    tr.mcMask());

      // Require true match only
      if (tr.mcMask() != 0) continue;

      const int mcId = tr.mcParticleId();
      if (mcId < 0) continue;

      auto const& mc = mcParticles.iteratorAt(mcId);

      // ---- Kinematics ---------------------------------------------------
      const double ptMC    = mc.pt();
      const double etaMC   = mc.eta();
      const double phiMC   = mc.phi();

      const double ptReco  = tr.pt();
      const double etaReco = tr.eta();
      const double phiReco = tr.phi();

      // ---- Residuals ----------------------------------------------------
      const double deltaPtRel = (ptMC > 0.) ? (ptReco - ptMC) / ptMC : -999.;
      const double deltaPtAbs = ptReco - ptMC;
      const double deltaEta   = etaReco - etaMC;

      double deltaPhi = phiReco - phiMC;
      while (deltaPhi >  M_PI) deltaPhi -= 2. * M_PI;
      while (deltaPhi < -M_PI) deltaPhi += 2. * M_PI;

      const double deltaInvPt = (ptMC > 0.) ? (1./ptReco - 1./ptMC) / (1./ptMC) : -999.;

      const int charge = (tr.invQPt() > 0.) ? 1 : -1;

      // ---- Fill resolution histos ---------------------------------------
      histos.fill(HIST("Resolution/hDeltaPtRel_vs_PtMC"),    ptMC, deltaPtRel);
      histos.fill(HIST("Resolution/hDeltaPtAbs_vs_PtMC"),    ptMC, deltaPtAbs);
      histos.fill(HIST("Resolution/hDeltaInvPtRel_vs_PtMC"), ptMC, deltaInvPt);
      histos.fill(HIST("Resolution/hDeltaEta_vs_EtaMC"),     etaMC, deltaEta);
      histos.fill(HIST("Resolution/hDeltaPhi_vs_PhiMC"),     phiMC, deltaPhi);
      histos.fill(HIST("Resolution/hDeltaPtRel_vs_EtaMC"),   etaMC, deltaPtRel);
      histos.fill(HIST("Resolution/hDeltaPtRel_PtMC_EtaMC"), ptMC, etaMC, deltaPtRel);

      if (charge > 0) {
        histos.fill(HIST("Resolution/hDeltaPtRel_vs_PtMC_Pos"), ptMC, deltaPtRel);
      } else {
        histos.fill(HIST("Resolution/hDeltaPtRel_vs_PtMC_Neg"), ptMC, deltaPtRel);
      }

      // ---- Fill migration matrices --------------------------------------
      histos.fill(HIST("Migration/hPtReco_vs_PtMC"),   ptMC,  ptReco);
      histos.fill(HIST("Migration/hEtaReco_vs_EtaMC"), etaMC, etaReco);
      histos.fill(HIST("Migration/hPhiReco_vs_PhiMC"), phiMC, phiReco);

      // ---- Sanity checks -----------------------------------------------
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
