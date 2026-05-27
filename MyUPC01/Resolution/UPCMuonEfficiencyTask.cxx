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

/// \file   UPCMuonEfficiencyTask.cxx
/// \brief  Efficiency analysis for dimuon pairs and single muons in UPC photoproduction.
///         Produces: 1D Efficiency histograms Eta Phi pT for single and dimuons
/// \author Takuma Matsumoto

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include <Framework/HistogramSpec.h>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TVector2.h"

#include <cmath>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ===========================================================================================

struct UPCMuonEfficiencyTask {
  Service<o2::framework::O2DatabasePDG> pdg;

  // Track types
  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  // Histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ===========================================================================
  // Configurables
  // ===========================================================================
  // Track types (setting for Global muon or Stand alone(under development))
  Configurable<int> reqTrackType{
    "reqTrackType",
    static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack),
    "0: GlobalMuon, 3: MuonStandalone"};

  // --- Event selection ---
  // Number of muon tracks in the events
  Configurable<int> reqMatchMFT{"reqMatchMFT", 2, "Required number of valid tracks per event"};

  // --- Acceptance region ---
  // Single Track (For MC true and MC reco)
  Configurable<float> etaMin{"etaMin", -3.5f, "Minimum eta"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum eta"};
  Configurable<float> pTmin{"pTmin", 0.3f, "Minimum pT [GeV/c]"}; // single pT > 0.3 GeV
  Configurable<float> pTmax{"pTmax", 10.f, "Maximum pT [GeV/c]"}; // Dont use for cut, only for Histogram bin edge

  // --- For Reco informations ---
  // track QA
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};

  // For Global track
  Configurable<float> maxChi2GlobalTrack{"maxChi2", 100.f, "Maximum allowed global track Chi2"};
  Configurable<float> maxChi2MatchMCHMFT{"maxChi2MatchMCHMFT", 100.f, "Maximum allowed MCH-MFT match Chi2"};

  // Dimuon Pair (For MC true and MC reco)
  Configurable<int> reqMuonPair{"reqMuonPair", 1, "Required number of muon pairs per event"}; // which is better use reqMatch MFT? or reqMuonPair?
  Configurable<float> pairRapMin{"pairRapMin", -3.5f, "Minimum pair rapidity"};
  Configurable<float> pairRapMax{"pairRapMax", -2.5f, "Maximum pair rapidity"};
  Configurable<float> pairMassMin{"pairMassMin", 1.0f, "Minimum pair mass [GeV/c^2]"};
  Configurable<float> pairMassMax{"pairMassMax", 10.0f, "Maximum pair mass [GeV/c^2]"};

  Configurable<float> pairPtMax{"pairPtMax", 10.0f, "Maximum pair pT [GeV/c]"}; // pair pT < 0.25 GeV coherent / 10 GeV for incoherent

  // --- Histogram binning ---
  // Phi
  Configurable<int> nBinsPhi{"nBinsPhi", 360, "Number of bins on phi axis"};
  // Eta and Rapidity
  Configurable<int> nBinsEta{"nBinsEta", 300, "Number of bins on eta axis"};
  // pT
  Configurable<int> nBinsPt{"nBinsPt", 1000, "Number of bins on pT axis"};
  Configurable<float> ptTrackMax{"ptTrackMax", 10.0f, "Upper edge of track pT axis [GeV/c]"};

  // Mass
  Configurable<int> nBinsMass{"nBinsMass", 300, "Number of bins on mass axis"};
  Configurable<float> massAxisMin{"massAxisMin", 1.0f, "Lower edge of mass axis [GeV/c^2]"};
  Configurable<float> massAxisMax{"massAxisMax", 10.0f, "Upper edge of mass axis [GeV/c^2]"};

  // Resolution
  // Under development

  // Muon PDG mass
  static constexpr int kMuonPDG = 13;
  float mMu = 0.105658f; // [GeV/c^2]

  // ===========================================================================
  // Histogram initialisation helpers
  // ===========================================================================
  void initControlHistograms()
  {
    // Cut flow
    auto hCutFlowMC = registry.add<TH1>("hCutFlowMC", "Selection Cut Flow;;Counts", HistType::kTH1I, {{18, 0., 18.}});
    auto hCutFlowReco = registry.add<TH1>("hCutFlowReco", "Selection Cut Flow;;Counts", HistType::kTH1I, {{18, 0., 18.}});
    registry.add<TH1>("hMuonMultMC", "Number of MC muons in acceptance per event;Number of muons;Events", HistType::kTH1I, {{10, 0., 10.}});
    registry.add<TH1>("hMuonMultReco", "Number of Reco muons in acceptance per event;Number of muons;Events", HistType::kTH1I, {{10, 0., 10.}});
    TString CutNames[18] = {
      "0: All Cand", // MC and Reco
      "1: Has Requested Track Type",
      "2: Pass Exact 2 Tracks( -> Pair ALL)",
      "3: Pass Exact MatchMFT",
      "4: Track All",
      "5: Track Pass Type",
      "6: Track Pass rAbs",
      "7: Track Pass pDCA",
      "8: Track Pass MatchMFT",
      "9: Track Pass Eta",    // MC and Reco
      "10: Track Pass pT",    // MC and Reco
      "11: Pair All",         // MC and Reco
      "12: Pair Unlike-sign", // MC and Reco
      "13: Pair Both MC Muon",
      "14: Pair Pass Pt",       // MC and Reco
      "15: Pair Pass Rapidity", // MC and Reco
      "16: Pair Pass Mass",     // MC and Reco
      "17: Filler"};
    for (int i = 0; i < 18; i++) {
      hCutFlowMC->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
      hCutFlowReco->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
    }
  }

  void initSingleTrackHistograms()
  {
    const AxisSpec axPdgCode{10000, -5000.f, 5000.f, "PDG Code"};

    registry.add("hPdgCodeMCAll", "PDG Codes of all MC particles", kTH1I, {axPdgCode});
    registry.add("hPdgCodeReco", "PDG Codes of Reco matched particles", kTH1I, {axPdgCode});

    const AxisSpec axPhiMCAll{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC All} (rad)"};
    const AxisSpec axPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC} (rad)"};
    const AxisSpec axPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{reco} (rad)"};
    const AxisSpec axPhiEff{120, 0.f, 1.2f, "#phi^{reco} / #phi^{MC}"};
    // const AxisSpec axPhiRes{200, -0.1f, 0.1f, "#phi^{reco} - #phi^{MC} (rad)"};

    const AxisSpec axEtaMCAll{nBinsEta, -10.f, 10.f, "#eta^{MC All}"};
    const AxisSpec axEtaMC{nBinsEta, etaMin, etaMax, "#eta^{MC}"};
    const AxisSpec axEtaReco{nBinsEta, etaMin, etaMax, "#eta^{reco}"};
    const AxisSpec axEtaEff{120, 0.f, 1.2f, "#eta^{reco} / #eta^{MC}"};
    // const AxisSpec axEtaRes{200, -0.1f, 0.1f, "#eta^{reco} - #eta^{MC}"};

    const AxisSpec axPtMCAll{nBinsPt, 0.f, 10.f, "pT^{MC All} (GeV/c)"};
    const AxisSpec axPtMC{nBinsPt, pTmin, pTmax, "pT^{MC} (GeV/c)"};
    const AxisSpec axPtReco{nBinsPt, pTmin, pTmax, "pT^{reco} (GeV/c)"};
    const AxisSpec axpTEff{120, 0.f, 1.2f, "pT^{reco} / pT^{MC}"};
    // const AxisSpec axPtRes{200, -0.1f, 0.1f, "pT^{reco} - pT^{MC} (GeV/c)"};

    // Histograms
    // 1D histograms
    registry.add("hTrkPhiMCAll", "Single track #phi MC All", kTH1F, {axPhiMCAll});
    registry.add("hTrkPhiMC", "Single track #phi MC", kTH1F, {axPhiMC});
    registry.add("hTrkPhiReco", "Single track #phi reco", kTH1F, {axPhiReco});
    registry.add("hTrkPhiEff", "Single track #phi eff", kTH1F, {axPhiEff});

    registry.add("hTrkEtaMCAll", "Single track #eta MC All", kTH1F, {axEtaMCAll});
    registry.add("hTrkEtaMC", "Single track #eta MC", kTH1F, {axEtaMC});
    registry.add("hTrkEtaReco", "Single track #eta reco", kTH1F, {axEtaReco});
    registry.add("hTrkEtaEff", "Single track #eta eff", kTH1F, {axEtaEff});

    registry.add("hTrkPtMCAll", "Single track pT MC All", kTH1F, {axPtMCAll});
    registry.add("hTrkPtMC", "Single track pT MC", kTH1F, {axPtMC});
    registry.add("hTrkPtReco", "Single track pT reco", kTH1F, {axPtReco});
    registry.add("hTrkPtEff", "Single track pT eff", kTH1F, {axpTEff});

    // 2D histograms
    // under development
  }

  void initDimuonHistograms()
  {
    const AxisSpec axPairPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{pair, MC} (rad)"};
    const AxisSpec axPairPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{pair, reco} (rad)"};
    const AxisSpec axPairPhiEff{120, 0.f, 1.2f, "#phi^{pair, reco} / #phi^{pair, MC}"};

    const AxisSpec axPairRapidityMC{nBinsEta, pairRapMin, pairRapMax, "y^{pair, MC}"};
    const AxisSpec axPairRapidityReco{nBinsEta, pairRapMin, pairRapMax, "y^{pair, reco}"};
    const AxisSpec axPairRapidityEff{120, 0.f, 1.2f, "y^{pair, reco} / y^{pair, MC}"};

    const AxisSpec axPairPtMC{nBinsPt, 0.f, pairPtMax, "pT^{pair, MC} (GeV/c)"};
    const AxisSpec axPairPtReco{nBinsPt, 0.f, pairPtMax, "pT^{pair, reco} (GeV/c)"};
    const AxisSpec axPairPtEff{120, 0.f, 1.2f, "pT^{pair, reco} / pT^{pair, MC}"};

    // Not need
    const AxisSpec axPairMassMC{nBinsMass, massAxisMin, massAxisMax, "m^{pair, MC} (GeV/c^2)"};
    const AxisSpec axPairMassReco{nBinsMass, massAxisMin, massAxisMax, "m^{pair, reco} (GeV/c^2)"};

    registry.add("hPairPhiMC", "Pair phi MC", kTH1F, {axPairPhiMC});
    registry.add("hPairPhiReco", "Pair phi reco", kTH1F, {axPairPhiReco});
    registry.add("hPairPhiEff", "Pair phi eff", kTH1F, {axPairPhiEff});

    registry.add("hPairRapidityMC", "Pair rapidity MC", kTH1F, {axPairRapidityMC});
    registry.add("hPairRapidityReco", "Pair rapidity reco", kTH1F, {axPairRapidityReco});
    registry.add("hPairRapidityEff", "Pair rapidity eff", kTH1F, {axPairRapidityEff});

    registry.add("hPairPtMC", "Pair pT MC", kTH1F, {axPairPtMC});
    registry.add("hPairPtReco", "Pair pT reco", kTH1F, {axPairPtReco});
    registry.add("hPairPtEff", "Pair pT eff", kTH1F, {axPairPtEff});

    registry.add("hPairMassMC", "Pair mass MC", kTH1F, {axPairMassMC});
    registry.add("hPairMassReco", "Pair mass reco", kTH1F, {axPairMassReco});

    // 2D histograms
    // under development
  }

  // ===========================================================================
  // Initialization
  // ===========================================================================
  void init(InitContext&)
  {
    auto* particle = TDatabasePDG::Instance()->GetParticle(kMuonPDG);
    mMu = particle ? particle->Mass() : 0.105658f;

    initControlHistograms();
    initSingleTrackHistograms();
    initDimuonHistograms();
  }

  // ===========================================================================
  // Core Methods
  // ===========================================================================

  // Collect reco track indices per collision (requires MC label)
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

  // Reco-level quality cuts for a single forward track
  template <typename TTrack>
  bool passRecoTrackCuts(const TTrack& tr)
  {
    registry.fill(HIST("hCutFlowReco"), 4); // 4: Track All

    float rAbs = tr.rAtAbsorberEnd();
    if (rAbs < rAbsMin || rAbs > rAbsMax)
      return false;
    registry.fill(HIST("hCutFlowReco"), 6); // 6: Track Pass rAbs

    float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
    if (tr.pDca() > pDcaMax)
      return false;
    registry.fill(HIST("hCutFlowReco"), 7); // 7: Track Pass pDCA

    if (tr.chi2() > maxChi2GlobalTrack)
      return false;

    if (reqTrackType == static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      float chi2MFT = tr.chi2MatchMCHMFT();
      if (chi2MFT < 0.f || chi2MFT > maxChi2MatchMCHMFT)
        return false;
    }
    registry.fill(HIST("hCutFlowReco"), 8); // 8: Track Pass MatchMFT

    TLorentzVector recoVec;
    recoVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    if (recoVec.Eta() <= etaMin || recoVec.Eta() >= etaMax)
      return false;
    registry.fill(HIST("hCutFlowReco"), 9); // 9: Track PassEta

    if (recoVec.Pt() < pTmin)
      return false;
    registry.fill(HIST("hCutFlowReco"), 10); // 10: Track Pass pT

    return true;
  }

  // ===========================================================================
  // Main process functions
  // ===========================================================================

  // ---------------------------------------------------------------------------
  // processMCTrue: MC-truth level
  //   Loops over UDMcParticles, selects muons within single-track acceptance,
  //   then forms unlike-sign pairs per MC collision and applies dimuon cuts.
  // ---------------------------------------------------------------------------
  void processMCTrue(o2::aod::UDMcCollisions const& mcCollisions,
                     o2::aod::UDMcParticles const& mcParticles)
  {
    // Group muon indices per MC collision
    std::unordered_map<int32_t, std::vector<int32_t>> muonsPerMcColl;

    for (const auto& mc : mcParticles) {
      registry.fill(HIST("hCutFlowMC"), 0); // 0: All Muon
      registry.fill(HIST("hPdgCodeMCAll"), mc.pdgCode());
      if (std::abs(mc.pdgCode()) != kMuonPDG)
        continue;
      registry.fill(HIST("hCutFlowMC"), 1); // 1: PDG=mu

      TLorentzVector v;
      v.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);
      registry.fill(HIST("hTrkPhiMCAll"), v.Phi());
      registry.fill(HIST("hTrkEtaMCAll"), v.Eta());
      registry.fill(HIST("hTrkPtMCAll"), v.Pt());

      // Single-track acceptance cuts (MC truth)
      if (v.Eta() <= etaMin || v.Eta() >= etaMax)
        continue;
      registry.fill(HIST("hCutFlowMC"), 9); // 9: Muon PassEta

      if (v.Pt() < pTmin)
        continue;
      registry.fill(HIST("hCutFlowMC"), 10); // 10: Muon Pass Pt

      // Fill single-track MC histograms
      registry.fill(HIST("hTrkPhiMC"), v.Phi());
      registry.fill(HIST("hTrkEtaMC"), v.Eta());
      registry.fill(HIST("hTrkPtMC"), v.Pt());

      muonsPerMcColl[mc.udMcCollisionId()].push_back(mc.globalIndex());
    }

    // Dimuon pair loop (MC truth)
    for (const auto& item : muonsPerMcColl) {
      const auto& ids = item.second;

      registry.fill(HIST("hMuonMultMC"), static_cast<float>(ids.size()));

      for (size_t i = 0; i < ids.size(); ++i) {
        auto mc1 = mcParticles.iteratorAt(ids[i]);
        TLorentzVector v1;
        v1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);

        for (size_t j = i + 1; j < ids.size(); ++j) {
          auto mc2 = mcParticles.iteratorAt(ids[j]);
          TLorentzVector v2;
          v2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);

          registry.fill(HIST("hCutFlowMC"), 11); // 11: Pair All

          // Unlike-sign (opposite PDG codes: +13 and -13)
          if (mc1.pdgCode() * mc2.pdgCode() >= 0)
            continue;
          registry.fill(HIST("hCutFlowMC"), 12); // 12: Pair Unlike-sign

          // Both are muons (guaranteed by outer loop, count for cutflow)
          registry.fill(HIST("hCutFlowMC"), 13); // 13: Pair Both MC Muon

          TLorentzVector pair = v1 + v2;

          // Pair kinematic cuts (using MC variables)
          if (pair.Pt() >= pairPtMax)
            continue;
          registry.fill(HIST("hCutFlowMC"), 14); // 14: Pair Pass Pt

          if (pair.Rapidity() <= pairRapMin || pair.Rapidity() >= pairRapMax)
            continue;
          registry.fill(HIST("hCutFlowMC"), 15); // 15: Pair Pass Rapidity

          if (pair.M() <= pairMassMin || pair.M() >= pairMassMax)
            continue;
          registry.fill(HIST("hCutFlowMC"), 16); // 16: Pair Pass Mass

          // Fill pair MC histograms
          registry.fill(HIST("hPairPhiMC"), pair.Phi());
          registry.fill(HIST("hPairRapidityMC"), pair.Rapidity());
          registry.fill(HIST("hPairPtMC"), pair.Pt());
          registry.fill(HIST("hPairMassMC"), pair.M());
        }
      }
    }

    (void)mcCollisions;
  }

  // ---------------------------------------------------------------------------
  // processReco: MC-reco level
  //   Loops over reco forward tracks with MC labels, applies detector quality
  //   cuts + kinematic cuts, fills Reco histograms.
  // ---------------------------------------------------------------------------
  void processReco(CandidatesFwd const& eventCandidates,
                   CompleteFwdTracks const& fwdTracks,
                   o2::aod::UDMcCollisions const&,
                   o2::aod::UDMcParticles const&)
  {
    // Group reco track indices by collision
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    for (const auto& item : tracksPerCand) {
      int32_t candId = item.first;
      const auto& trkIds = item.second;

      if (trkIds.empty())
        continue;

      registry.fill(HIST("hCutFlowReco"), 0); // 0: All Cand

      // Filter by requested track type
      std::vector<int32_t> candidateTrkIds;
      candidateTrkIds.reserve(trkIds.size());

      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        if (tr.trackType() == reqTrackType) {
          candidateTrkIds.push_back(idx);
        }
      }

      if (candidateTrkIds.empty())
        continue;
      registry.fill(HIST("hCutFlowReco"), 1); // 1: Has Requested Track Type

      // Apply reco quality + kinematic cuts per track
      std::vector<int32_t> goodTrkIds;
      goodTrkIds.reserve(reqMatchMFT);

      for (auto idx : candidateTrkIds) {
        auto tr = fwdTracks.iteratorAt(idx);

        if (!passRecoTrackCuts(tr))
          continue;
        if (!tr.has_udMcParticle())
          continue;
        const auto& mc = tr.udMcParticle();

        registry.fill(HIST("hPdgCodeReco"), mc.pdgCode());
        if (std::abs(mc.pdgCode()) != kMuonPDG)
          continue;

        goodTrkIds.push_back(idx);

        // Fill single-track Reco histograms
        TLorentzVector vReco;
        vReco.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
        registry.fill(HIST("hTrkPhiReco"), vReco.Phi());
        registry.fill(HIST("hTrkEtaReco"), vReco.Eta());
        registry.fill(HIST("hTrkPtReco"), vReco.Pt());
      }
      registry.fill(HIST("hMuonMultReco"), static_cast<float>(goodTrkIds.size()));

      //if (candidateTrkIds.size() != static_cast<size_t>(reqMatchMFT))
        //continue;
      registry.fill(HIST("hCutFlowReco"), 2); // 2: Pass Exact 2 Tracks

      if (goodTrkIds.size() != static_cast<size_t>(reqMatchMFT))
        continue;
      registry.fill(HIST("hCutFlowReco"), 3); // 3: Pass Exact MatchMFT

      // Dimuon pair loop (Reco)
      for (size_t i = 0; i < goodTrkIds.size(); ++i) {
        auto tr1 = fwdTracks.iteratorAt(goodTrkIds[i]);
        TLorentzVector vReco1;
        vReco1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);

        for (size_t j = i + 1; j < goodTrkIds.size(); ++j) {
          auto tr2 = fwdTracks.iteratorAt(goodTrkIds[j]);
          TLorentzVector vReco2;
          vReco2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);

          registry.fill(HIST("hCutFlowReco"), 11); // 11: Pair All

          // Unlike-sign
          if (tr1.sign() * tr2.sign() >= 0)
            continue;
          registry.fill(HIST("hCutFlowReco"), 12); // 12: Pair Unlike-sign

          // Both tracks are MC muons (already filtered above)
          registry.fill(HIST("hCutFlowReco"), 13); // 13: Pair Both MC Muon

          TLorentzVector pairReco = vReco1 + vReco2;

          // Pair kinematic cuts (using reco variables)
          if (pairReco.Pt() >= pairPtMax)
            continue;
          registry.fill(HIST("hCutFlowReco"), 14); // 14: Pair Pass Pt

          if (pairReco.Rapidity() <= pairRapMin || pairReco.Rapidity() >= pairRapMax)
            continue;
          registry.fill(HIST("hCutFlowReco"), 15); // 15: Pair Pass Rapidity

          if (pairReco.M() <= pairMassMin || pairReco.M() >= pairMassMax)
            continue;
          registry.fill(HIST("hCutFlowReco"), 16); // 16: Pair Pass Mass

          // Fill pair Reco histograms
          registry.fill(HIST("hPairPhiReco"), pairReco.Phi());
          registry.fill(HIST("hPairRapidityReco"), pairReco.Rapidity());
          registry.fill(HIST("hPairPtReco"), pairReco.Pt());
          registry.fill(HIST("hPairMassReco"), pairReco.M());
        }
      }

      (void)candId;
      (void)eventCandidates;
    }
  }

  PROCESS_SWITCH(UPCMuonEfficiencyTask, processMCTrue, "Fill MC-truth single-track and dimuon histograms", true);
  PROCESS_SWITCH(UPCMuonEfficiencyTask, processReco, "Fill reco single-track and dimuon histograms", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCMuonEfficiencyTask>(cfgc, TaskName{"upc-muon-efficiency"}),
  };
}
