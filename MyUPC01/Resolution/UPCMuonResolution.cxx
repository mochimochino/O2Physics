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

/// \file   UPCMuonResolution.cxx
/// \brief  pT resolution analysis for forward muons in UPC events.
///         Reads MC reco information (udFwdTracks + udMcParticles) and fills
///         two 2D histograms:
///           1) hPtResoVsPtTrue  : (pTreco - pTtrue) / pTtrue  vs  pTtrue
///           2) hResponseMatrix  : pTreco  vs  pTtrue  (for unfolding)
/// \author Takuma Matsumoto

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include <unordered_map>
#include <vector>
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ============================================================================

struct UPCMuonResolution {

  Service<o2::framework::O2DatabasePDG> pdg;

  // Track type aliases  -------------------------------------------------------
  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  // Histogram registries  -----------------------------------------------------
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurables  ------------------------------------------------------------
  Configurable<float> etaMin{"etaMin", -4.0f, "Minimum pseudorapidity for muon tracks"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum pseudorapidity for muon tracks"};
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};
  Configurable<int>   nBinsPtTrue{"nBinsPtTrue", 100, "Number of bins on the pT_true axis"};
  Configurable<float> ptTrueMax{"ptTrueMax", 10.0f, "Maximum pT_true [GeV/c]"};
  Configurable<int>   nBinsRelRes{"nBinsRelRes", 200, "Number of bins on the relative-resolution axis"};
  Configurable<float> relResMax{"relResMax", 0.5f, "Half-range of relative-resolution axis"};
  Configurable<int>   nBinsPtReco{"nBinsPtReco", 100, "Number of bins on the pT_reco axis"};
  Configurable<float> ptRecoMax{"ptRecoMax", 10.0f, "Maximum pT_reco [GeV/c]"};

  static constexpr int kMuonPDG = 13;

  // ---------------------------------------------------------------------------
  void init(InitContext&)
  {
    // --- Axis definitions ---
    // X axis: pT_true
    const AxisSpec axisPtTrue{nBinsPtTrue, 0.f, ptTrueMax, "#it{p}_{T}^{true} (GeV/#it{c})"};
    // Y axis for resolution histogram: (pTreco - pTtrue) / pTtrue
    const AxisSpec axisRelRes{nBinsRelRes, -relResMax, relResMax, "(#it{p}_{T}^{reco} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}"};
    // Y axis for response matrix: pT_reco
    const AxisSpec axisPtReco{nBinsPtReco, 0.f, ptRecoMax, "#it{p}_{T}^{reco} (GeV/#it{c})"};

    // Event counter
    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("eventCounter", "Processed Events (MC Reco)", kTH1F, {axisCounter});

    // --- 2D Histograms ---
    // (1) Resolution: Y = (pTreco - pTtrue) / pTtrue  vs  X = pTtrue
    registry.add("hPtResoVsPtTrue",
                 "pT Resolution;#it{p}_{T}^{true} (GeV/#it{c});(#it{p}_{T}^{reco} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}",
                 kTH2F, {axisPtTrue, axisRelRes});

    // (2) Response matrix: Y = pTreco  vs  X = pTtrue
    registry.add("hResponseMatrix",
                 "Response Matrix;#it{p}_{T}^{true} (GeV/#it{c});#it{p}_{T}^{reco} (GeV/#it{c})",
                 kTH2F, {axisPtTrue, axisPtReco});

    // 1D projections for quick inspection
    registry.add("hPtTrue",  "pT True;#it{p}_{T}^{true} (GeV/#it{c});Counts",  kTH1F, {axisPtTrue});
    registry.add("hPtReco",  "pT Reco;#it{p}_{T}^{reco} (GeV/#it{c});Counts",  kTH1F, {axisPtReco});
  }

  // ---------------------------------------------------------------------------
  /// Collect track global-indices per UDCollision ID
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId >= 0 && tr.has_udMcParticle()) {
        tracksPerCand[candId].push_back(tr.globalIndex());
      }
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply quality cuts on a single track
  template <typename TTrack>
  bool passTrackCuts(const TTrack& tr)
  {
    // R at absorber end
    float rAbs = tr.rAtAbsorberEnd();
    if (rAbs < rAbsMin || rAbs > rAbsMax)
      return false;

    // pDCA cut
    float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
    if (tr.pDca() > pDcaMax)
      return false;

    // MCH-MID match required
    if (tr.chi2MatchMCHMID() <= 0)
      return false;

    return true;
  }

  // ---------------------------------------------------------------------------
  /// Fill resolution histograms for a single (reco, mc) track pair
  template <typename TTrack, typename TMcPart>
  void fillResolution(const TTrack& tr, const TMcPart& mc)
  {
    // Verify the MC particle is a muon
    if (std::abs(mc.pdgCode()) != kMuonPDG)
      return;

    // Build 4-vectors
    const float mMu = 0.10566f; // muon mass in GeV
    TLorentzVector recoVec, trueVec;
    recoVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    trueVec.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);

    float etaReco = recoVec.Eta();
    if (etaReco <= etaMin || etaReco >= etaMax)
      return;

    float ptTrue = trueVec.Pt();
    float ptReco = recoVec.Pt();

    if (ptTrue <= 0.f)
      return;

    float relRes = (ptReco - ptTrue) / ptTrue;

    // Fill histograms
    registry.fill(HIST("hPtResoVsPtTrue"), ptTrue, relRes);
    registry.fill(HIST("hResponseMatrix"), ptTrue, ptReco);
    registry.fill(HIST("hPtTrue"), ptTrue);
    registry.fill(HIST("hPtReco"), ptReco);
  }

  // ---------------------------------------------------------------------------
  void processMcReco(CandidatesFwd const& eventCandidates,
                     CompleteFwdTracks const& fwdTracks,
                     o2::aod::UDMcCollisions const&,
                     o2::aod::UDMcParticles const& /*McParts*/)
  {
    registry.fill(HIST("eventCounter"), 0.5);

    // Collect tracks (with MC truth link) per candidate
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    for (const auto& item : tracksPerCand) {
      int32_t candId = item.first;
      const auto& trkIds = item.second;
      if (trkIds.size() < 1)
        continue;

      // Count MCH-MID matched tracks (quality requirement for the event)
      int nMchMid = 0;
      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        if (tr.chi2MatchMCHMID() > 0)
          nMchMid++;
      }
      // Require exactly 2 MCH-MID tracks for dimuon candidates
      if (nMchMid < 2)
        continue;

      // Loop over individual tracks
      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        if (!passTrackCuts(tr))
          continue;
        const auto& mc = tr.udMcParticle();
        fillResolution(tr, mc);
      }

      (void)candId; // suppress unused-variable warning
      (void)eventCandidates;
    }
  }
  PROCESS_SWITCH(UPCMuonResolution, processMcReco, "Fill resolution histograms from MC reco data", true);
};

// ============================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCMuonResolution>(cfgc, TaskName{"my-upc-muon-resolution"}),
  };
}
