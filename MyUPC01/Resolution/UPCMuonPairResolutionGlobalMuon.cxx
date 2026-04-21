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

/// \file   UPCMuonPairResolution.cxx
/// \brief  Resolution analysis for the dimuon pair pT and pT^2 in UPC J/psi photoproduction.
///         Since |t| ≈ pT^2 for coherent UPC events, pT^2 resolution directly reflects
///         the momentum-transfer resolution critical for nuclear structure studies.
///         Produces: 1D Pair pT (Reco/MC), 1D Pair pT^2 (Reco/MC), and 2D response
///         matrix (pT^2_reco vs pT^2_MC) as well as a cut-flow histogram.
/// \author Takuma Matsumoto

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ============================================================================

struct UPCMuonPairResolution {

  Service<o2::framework::O2DatabasePDG> pdg;

  // Track type aliases  -------------------------------------------------------
  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  // Histogram registry  -------------------------------------------------------
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ===========================================================================
  // Configurables
  // ===========================================================================
  // Event Level
  Configurable<int> reqMatchMFT{"reqMatchMFT", 2, "Required number of MCH-MFT matched tracks"};

  // Track Type (0: GlobalMuon, 3: Standalone)
  Configurable<int> reqTrackType{"reqTrackType", 0, "Track type required for the analysis"};

  // Single Track Level
  Configurable<float> etaMin{"etaMin", -4.0f, "Minimum pseudorapidity for muon tracks"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum pseudorapidity for muon tracks"};
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};

  // Dimuon Pair Level
  Configurable<float> pairPtMax{"pairPtMax", 10.0f, "Maximum dimuon pair pT [GeV/c]"}; // 10.0
  Configurable<float> pairRapidityMin{"pairRapidityMin", -4.0f, "Minimum dimuon pair rapidity"};
  Configurable<float> pairRapidityMax{"pairRapidityMax", -2.5f, "Maximum dimuon pair rapidity"};
  Configurable<float> pairMassMin{"pairMassMin", 1.0f, "Minimum dimuon pair mass [GeV/c^2]"};
  Configurable<float> pairMassMax{"pairMassMax", 10.0f, "Maximum dimuon pair mass [GeV/c^2]"};

  // Histogram Binning
  Configurable<int> nBinsPt{"nBinsPt", 1000, "Number of bins on pT axis"};
  Configurable<float> ptMax{"ptMax", 5.0f, "Upper edge of pT axis [GeV/c]"}; // 5
  Configurable<int> nBinsPt2{"nBinsPt2", 5000, "Number of bins on pT^2 axis"};
  Configurable<float> pt2Max{"pt2Max", 2.5f, "Upper edge of pT^2 axis [GeV^2/c^2]"}; // 2.5

  static constexpr int kMuonPDG = 13;
  float mMu = 0.0f;

  // ---------------------------------------------------------------------------
  void init(InitContext&)
  {
    auto* particle = TDatabasePDG::Instance()->GetParticle(kMuonPDG);
    mMu = particle ? particle->Mass() : 0.105658f;

    // --- Cut Flow Histogram ---
    auto hCutFlow = registry.add<TH1>("hCutFlow", "Selection Cut Flow;;Counts", HistType::kTH1I, {{16, 0., 16.}});
    TString CutNames[16] = {
      "0: All Cand",            // 0
      "1: Pass Exact MatchMFT", // 1
      "2: Track All",           // 2
      "3: Track Pass Type",     // 3  <- トラックタイプの選択を最初に配置
      "4: Track Pass rAbs",     // 4
      "5: Track Pass pDCA",     // 5
      "6: Track Pass MatchMFT", // 6  <- MCH-MIDからMCH-MFTのチェックに変更
      "7: Track Pass Eta",      // 7
      "8: Pair All",            // 8
      "9: Pair Unlike-sign",    // 9
      "10: Pair Both MC Muon",  // 10
      "11: Pair Pass Pt",       // 11
      "12: Pair Pass Rapidity", // 12
      "13: Pair Pass Mass",     // 13
      "14: Filler",             // 14
      "15: Filler"              // 15
    };
    for (int i = 0; i < 16; i++) {
      hCutFlow->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
    }

    // --- Axis definitions ---
    // Pair pT  (wide range for pre-cut; tight range re-used for post-cut)
    const AxisSpec axisPairPtMC{nBinsPt, 0.f, ptMax, "#it{p}_{T,#mu#mu}^{MC} (GeV/#it{c})"};
    const AxisSpec axisPairPtReco{nBinsPt, 0.f, ptMax, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};

    // Pair pT^2
    const AxisSpec axisPairPt2MC{nBinsPt2, 0.f, pt2Max, "#it{p}_{T,#mu#mu}^{2,MC} (GeV^{2}/#it{c}^{2})"};
    const AxisSpec axisPairPt2Reco{nBinsPt2, 0.f, pt2Max, "#it{p}_{T,#mu#mu}^{2,reco} (GeV^{2}/#it{c}^{2})"};

    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});

    // =========================================================================
    // 1D Histograms: Pair pT  -- PreCut (after unlike-sign+MC-muon, before kin. cuts)
    // =========================================================================
    registry.add("hPairPtMC_PreCut", "Dimuon Pair #it{p}_{T} MC truth (Pre Kin. Cuts)", kTH1F, {axisPairPtMC});
    registry.add("hPairPtReco_PreCut", "Dimuon Pair #it{p}_{T} Reco (Pre Kin. Cuts)", kTH1F, {axisPairPtReco});

    // 1D Histograms: Pair pT  -- PostCut (all cuts passed)
    registry.add("hPairPtMC_PostCut", "Dimuon Pair #it{p}_{T} MC truth (Post All Cuts)", kTH1F, {axisPairPtMC});
    registry.add("hPairPtReco_PostCut", "Dimuon Pair #it{p}_{T} Reco (Post All Cuts)", kTH1F, {axisPairPtReco});

    // =========================================================================
    // 1D Histograms: Pair pT^2  (|t| ≈ pT^2)  -- PreCut
    // =========================================================================
    registry.add("hPairPt2MC_PreCut", "Dimuon Pair #it{p}_{T}^{2} MC truth (Pre Kin. Cuts)", kTH1F, {axisPairPt2MC});
    registry.add("hPairPt2Reco_PreCut", "Dimuon Pair #it{p}_{T}^{2} Reco (Pre Kin. Cuts)", kTH1F, {axisPairPt2Reco});

    // 1D Histograms: Pair pT^2  -- PostCut
    registry.add("hPairPt2MC_PostCut", "Dimuon Pair #it{p}_{T}^{2} MC truth (Post All Cuts)", kTH1F, {axisPairPt2MC});
    registry.add("hPairPt2Reco_PostCut", "Dimuon Pair #it{p}_{T}^{2} Reco (Post All Cuts)", kTH1F, {axisPairPt2Reco});

    // =========================================================================
    // 2D Response Matrices (PostCut only)
    // =========================================================================
    // pT response matrix
    registry.add("hResponseMatrixPairPt",
                 "Pair #it{p}_{T} Response Matrix (#it{p}_{T}^{reco} vs #it{p}_{T}^{MC})",
                 kTH2F, {axisPairPtMC, axisPairPtReco});

    // pT^2 response matrix
    registry.add("hResponseMatrixPairPt2",
                 "Pair #it{p}_{T}^{2} Response Matrix (#it{p}_{T}^{2,reco} vs #it{p}_{T}^{2,MC})",
                 kTH2F, {axisPairPt2MC, axisPairPt2Reco});

    // =========================================================================
    // 2D Relative Residuals (PostCut only)
    // =========================================================================
    const AxisSpec axisPtRelRes{200, -1.0f, 1.0f,
                                "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};
    registry.add("hPairPtResoVsPtMC",
                 "Pair #it{p}_{T} Resolution vs #it{p}_{T}^{MC}",
                 kTH2F, {axisPairPtMC, axisPtRelRes});

    const AxisSpec axisPt2RelRes{200, -5.0f, 30.0f,
                                 "(#it{p}_{T}^{2,reco} - #it{p}_{T}^{2,MC}) / #it{p}_{T}^{2,MC}"};
    registry.add("hPairPt2ResoVsPt2MC",
                 "Pair #it{p}_{T}^{2} Resolution vs #it{p}_{T}^{2,MC}",
                 kTH2F, {axisPairPt2MC, axisPt2RelRes});
  }

  // ---------------------------------------------------------------------------
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand,
                      TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId >= 0 && tr.has_udMcParticle()) {
        tracksPerCand[candId].push_back(tr.globalIndex());
      }
    }
  }

  // ---------------------------------------------------------------------------
  template <typename TTrack>
  bool passTrackCuts(const TTrack& tr)
  {
    registry.fill(HIST("hCutFlow"), 2); // 2: Track All

    // 1. トラックタイプの一致を確認
    if (tr.trackType() != reqTrackType)
      return false;
    registry.fill(HIST("hCutFlow"), 3); // 3: Track Pass Type

    // 2. rAbs のカット
    float rAbs = tr.rAtAbsorberEnd();
    if (rAbs < rAbsMin || rAbs > rAbsMax)
      return false;
    registry.fill(HIST("hCutFlow"), 4); // 4: Track Pass rAbs

    // 3. pDCA のカット
    float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
    if (tr.pDca() > pDcaMax)
      return false;
    registry.fill(HIST("hCutFlow"), 5); // 5: Track Pass pDCA

    // 4. MCH-MFT マッチングのチェック
    if (tr.chi2MatchMCHMFT() <= 0)
      return false;
    registry.fill(HIST("hCutFlow"), 6); // 6: Track Pass MatchMFT

    // 5. Eta のカット
    TLorentzVector recoVec;
    recoVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    if (recoVec.Eta() <= etaMin || recoVec.Eta() >= etaMax)
      return false;
    registry.fill(HIST("hCutFlow"), 7); // 7: Track Pass Eta

    return true;
  }

  // ---------------------------------------------------------------------------
  void processMcReco(CandidatesFwd const& eventCandidates,
                     CompleteFwdTracks const& fwdTracks,
                     o2::aod::UDMcCollisions const&,
                     o2::aod::UDMcParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    for (const auto& item : tracksPerCand) {
      int32_t candId = item.first;
      const auto& trkIds = item.second;
      if (trkIds.size() < 1)
        continue;

      registry.fill(HIST("hCutFlow"), 0); // 0: All Cand

      // --- Event-level cut: exactly reqMatchMFT matched tracks of requested type ---
      int nValidTracks = 0;
      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        // TrackType の合致と、MCH-MFTの有効なマッチング(IDやChi2)を持つトラックをカウント
        if (tr.trackType() == reqTrackType && tr.chi2MatchMCHMFT() > 0) {
          nValidTracks++;
        }
      }
      if (nValidTracks != reqMatchMFT)
        continue;
      registry.fill(HIST("hCutFlow"), 1); // 1: Pass Exact MatchMFT

      // --- Track-level cuts ---
      std::vector<int32_t> goodTrkIds;
      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        if (passTrackCuts(tr)) {
          goodTrkIds.push_back(idx);
        }
      }

      // --- Dimuon pair loop ---
      for (size_t i = 0; i < goodTrkIds.size(); ++i) {
        auto tr1 = fwdTracks.iteratorAt(goodTrkIds[i]);
        const auto& mc1 = tr1.udMcParticle();

        for (size_t j = i + 1; j < goodTrkIds.size(); ++j) {
          auto tr2 = fwdTracks.iteratorAt(goodTrkIds[j]);
          const auto& mc2 = tr2.udMcParticle();

          registry.fill(HIST("hCutFlow"), 8); // 8: Pair All

          // Unlike-sign cut
          if (tr1.sign() * tr2.sign() >= 0)
            continue;
          registry.fill(HIST("hCutFlow"), 9); // 9: Pair Unlike-sign

          // Both tracks must be MC muons
          if (std::abs(mc1.pdgCode()) != kMuonPDG ||
              std::abs(mc2.pdgCode()) != kMuonPDG)
            continue;
          registry.fill(HIST("hCutFlow"), 10); // 10: Pair Both MC Muon

          // Build four-vectors
          TLorentzVector recoVec1, recoVec2;
          recoVec1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
          recoVec2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
          TLorentzVector pairReco = recoVec1 + recoVec2;

          TLorentzVector mcVec1, mcVec2;
          mcVec1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);
          mcVec2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);
          TLorentzVector pairMC = mcVec1 + mcVec2;

          float pairPtMC = pairMC.Pt();
          float pairPtReco = pairReco.Pt();
          float pairPt2MC = pairPtMC * pairPtMC;
          float pairPt2Reco = pairPtReco * pairPtReco;

          // =================================================================
          // PreCut fill  (after unlike-sign + both-MC-muon, before kin. cuts)
          // =================================================================
          registry.fill(HIST("hPairPtMC_PreCut"), pairPtMC);
          registry.fill(HIST("hPairPtReco_PreCut"), pairPtReco);
          registry.fill(HIST("hPairPt2MC_PreCut"), pairPt2MC);
          registry.fill(HIST("hPairPt2Reco_PreCut"), pairPt2Reco);

          // =================================================================
          // Pair kinematic cuts (applied on Reco)
          // =================================================================
          if (pairReco.Pt() >= pairPtMax)
            continue;
          registry.fill(HIST("hCutFlow"), 11); // 11: Pair Pass Pt

          if (pairReco.Rapidity() <= pairRapidityMin ||
              pairReco.Rapidity() >= pairRapidityMax)
            continue;
          registry.fill(HIST("hCutFlow"), 12); // 12: Pair Pass Rapidity

          if (pairReco.M() <= pairMassMin || pairReco.M() >= pairMassMax)
            continue;
          registry.fill(HIST("hCutFlow"), 13); // 13: Pair Pass Mass

          // =================================================================
          // PostCut fill  (all cuts passed)
          // =================================================================
          registry.fill(HIST("hPairPtMC_PostCut"), pairPtMC);
          registry.fill(HIST("hPairPtReco_PostCut"), pairPtReco);
          registry.fill(HIST("hPairPt2MC_PostCut"), pairPt2MC);
          registry.fill(HIST("hPairPt2Reco_PostCut"), pairPt2Reco);

          // pT response matrix
          registry.fill(HIST("hResponseMatrixPairPt"), pairPtMC, pairPtReco);
          // pT relative residual
          if (pairPtMC > 0.f) {
            registry.fill(HIST("hPairPtResoVsPtMC"), pairPtMC,
                          (pairPtReco - pairPtMC) / pairPtMC);
          }

          // pT^2 response matrix
          registry.fill(HIST("hResponseMatrixPairPt2"), pairPt2MC, pairPt2Reco);
          // pT^2 relative residual
          if (pairPt2MC > 0.f) {
            registry.fill(HIST("hPairPt2ResoVsPt2MC"), pairPt2MC,
                          (pairPt2Reco - pairPt2MC) / pairPt2MC);
          }
        } // j
      } // i

      (void)candId;
      (void)eventCandidates;
    }
  }
  PROCESS_SWITCH(UPCMuonPairResolution, processMcReco,
                 "Fill pair pT and pT^2 resolution histograms from MC reco data", true);
};

// ============================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCMuonPairResolution>(cfgc, TaskName{"my-upc-muon-pair-resolution"}),
  };
}
