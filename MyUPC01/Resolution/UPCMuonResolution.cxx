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
/// \brief  Resolution analysis for forward muons (pT, eta, phi, mass) in UPC events.
///         Includes configurable kinematics cuts, cut-flow, and post-cut single muon resolutions.
/// \author Takuma Matsumoto (Modified)

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TString.h"

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

  // ===========================================================================
  // Configurables (Cut conditions managed here)
  // ===========================================================================
  // Event Level
  Configurable<int>   reqMatchMID{"reqMatchMID", 2, "Required number of MCH-MID matched tracks (use 2 for exact match)"};
  
  // Single Track Level
  Configurable<float> etaMin{"etaMin", -4.0f, "Minimum pseudorapidity for muon tracks"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum pseudorapidity for muon tracks"};
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};
  
  // Dimuon Pair Level
  Configurable<float> pairPtMax{"pairPtMax", 0.25f, "Maximum dimuon pair pT [GeV/c]"};
  Configurable<float> pairRapidityMin{"pairRapidityMin", -4.0f, "Minimum dimuon pair rapidity"};
  Configurable<float> pairRapidityMax{"pairRapidityMax", -2.5f, "Maximum dimuon pair rapidity"};
  Configurable<float> pairMassMin{"pairMassMin", 1.0f, "Minimum dimuon pair mass [GeV/c^2]"};
  Configurable<float> pairMassMax{"pairMassMax", 10.0f, "Maximum dimuon pair mass [GeV/c^2]"};

  // Histogram Binning Settings
  Configurable<int>   nBinsPt{"nBinsPt", 1000, "Number of bins on pT axis"};
  Configurable<int>   nBinsPt2{"nBinsPt2", 1000, "Number of bins on pT^2 axis"};
  Configurable<float> ptMax{"ptMax", 10.0f, "Maximum single pT [GeV/c]"};
  Configurable<float> pt2Max{"pt2Max", 1.0f, "Maximum single pT^2 [GeV^2/c^2]"};
  Configurable<int>   nBinsMass{"nBinsMass", 200, "Number of bins on Mass axis"};

  static constexpr int kMuonPDG = 13;
  float mMu = 0.0f; 

  // ---------------------------------------------------------------------------
  void init(InitContext&)
  {
    auto* particle = TDatabasePDG::Instance()->GetParticle(kMuonPDG);
    mMu = particle ? particle->Mass() : 0.105658f; 

    // --- Cut Flow Histogram ---
    auto hCutFlow = registry.add<TH1>("hCutFlow", "Selection Cut Flow;;Counts", HistType::kTH1I, {{15, 0., 15.}});
    TString CutNames[15] = {
      "0: All Cand",             // 0
      "1: Pass Exact MatchMID",  // 1
      "2: Track All",            // 2
      "3: Track Pass rAbs",      // 3
      "4: Track Pass pDCA",      // 4
      "5: Track Pass MatchMID",  // 5
      "6: Track Pass Eta",       // 6
      "7: Pair All",             // 7
      "8: Pair Unlike-sign",     // 8
      "9: Pair Both MC Muon",    // 9
      "10: Pair Pass Pt",        // 10
      "11: Pair Pass Rapidity",  // 11
      "12: Pair Pass Mass",      // 12
      "13: Filler",              // 13
      "14: Filler"               // 14
    };
    for (int i = 0; i < 15; i++) {
      hCutFlow->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
    }

    // --- Axis definitions ---
    const AxisSpec axisPtTrue{nBinsPt, 0.f, ptMax, "#it{p}_{T}^{true} (GeV/#it{c})"};
    const AxisSpec axisPtReco{nBinsPt, 0.f, ptMax, "#it{p}_{T}^{reco} (GeV/#it{c})"};
    const AxisSpec axisPtRelRes{200, -0.5f, 0.5f, "(#it{p}_{T}^{reco} - #it{p}_{T}^{true}) / #it{p}_{T}^{true}"};
    const AxisSpec axisPairPtPreCut{500, 0.f, 5.0f, "#it{p}_{T,#mu#mu} (GeV/#it{c})"};
    const AxisSpec axisPairPtReco{nBinsPt, 0.f, pairPtMax, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};

    const AxisSpec axisEtaTrue{100, -4.5f, -2.0f, "#eta^{true}"};
    const AxisSpec axisEtaReco{100, -4.5f, -2.0f, "#eta^{reco}"};
    const AxisSpec axisEtaRes{100, -0.25f, 0.25f, "#eta^{reco} - #eta^{true}"};

    const AxisSpec axisPhiTrue{100, -TMath::Pi(), TMath::Pi(), "#varphi^{true} (rad)"};
    const AxisSpec axisPhiReco{100, -TMath::Pi(), TMath::Pi(), "#varphi^{reco} (rad)"};
    const AxisSpec axisPhiRes{100, -0.25f, 0.25f, "#varphi^{reco} - #varphi^{true}"};

    const AxisSpec axisMassTrue{nBinsMass, pairMassMin, pairMassMax, "#it{M}_{#mu#mu}^{true} (GeV/#it{c}^{2})"};
    const AxisSpec axisMassReco{nBinsMass, pairMassMin, pairMassMax, "#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})"};
    const AxisSpec axisMassRelRes{200, -0.2f, 0.2f, "(#it{M}^{reco} - #it{M}^{true}) / #it{M}^{true}"};

    const AxisSpec axisPt2True{nBinsPt2, 0.f, pt2Max, "#it{p}_{T,true}^{2} (GeV^{2}/#it{c}^{2})"};
    const AxisSpec axisPt2RelRes{nBinsPt2, -5.0f, 30.0f, "(#it{p}_{T,reco}^{2} - #it{p}_{T,true}^{2}) / #it{p}_{T,true}^{2}"};

    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});

    // --- 1D & 2D Histograms (Single Muon - Before Pair Cuts) ---
    registry.add("hPtTrue",  "pT True",  kTH1F, {axisPtTrue});
    registry.add("hPtReco",  "pT Reco",  kTH1F, {axisPtReco});
    registry.add("hEtaTrue", "Eta True", kTH1F, {axisEtaTrue});
    registry.add("hEtaReco", "Eta Reco", kTH1F, {axisEtaReco});
    registry.add("hPhiTrue", "Phi True", kTH1F, {axisPhiTrue});
    registry.add("hPhiReco", "Phi Reco", kTH1F, {axisPhiReco});

    registry.add("hPtResoVsPtTrue",  "pT Resolution",  kTH2F, {axisPtTrue, axisPtRelRes});
    registry.add("hEtaResoVspTTrue", "Eta Resolution", kTH2F, {axisPtTrue, axisEtaRes});
    
    registry.add("hPhiResoVsPtTrue",     "Phi Resolution (All)",      kTH2F, {axisPtTrue, axisPhiRes});
    registry.add("hPhiResoVsPtTrue_Pos", "Phi Resolution (Positive)", kTH2F, {axisPtTrue, axisPhiRes});
    registry.add("hPhiResoVsPtTrue_Neg", "Phi Resolution (Negative)", kTH2F, {axisPtTrue, axisPhiRes});
    
    registry.add("hResponseMatrixPt",  "pT Response",  kTH2F, {axisPtTrue, axisPtReco});
    registry.add("hResponseMatrixEta", "Eta Response", kTH2F, {axisEtaTrue, axisEtaReco});
    registry.add("hResponseMatrixPhi", "Phi Response", kTH2F, {axisPhiTrue, axisPhiReco});

    // --- 2D Histograms (Single Muon - AFTER Pair Cuts) ---
    registry.add("hPtResoVsPtTrue_PostCut",  "pT Resolution (Post Pair Cuts)",  kTH2F, {axisPtTrue, axisPtRelRes});
    registry.add("hEtaResoVspTTrue_PostCut", "Eta Resolution (Post Pair Cuts)", kTH2F, {axisPtTrue, axisEtaRes});
    registry.add("hPhiResoVsPtTrue_PostCut", "Phi Resolution (Post Pair Cuts)", kTH2F, {axisPtTrue, axisPhiRes});

    // --- Histograms (Unlike-sign Dimuon Mass) ---
    registry.add("hMassTrue", "Dimuon Mass True", kTH1F, {axisMassTrue});
    registry.add("hMassReco", "Dimuon Mass Reco", kTH1F, {axisMassReco});
    registry.add("hMassResoVsMassTrue", "Mass Resolution", kTH2F, {axisMassTrue, axisMassRelRes});
    registry.add("hResponseMatrixMass", "Mass Response",   kTH2F, {axisMassTrue, axisMassReco});

    // --- Dimuon Pair pT ---
    registry.add("hPairPtMC_PreCut", "Dimuon Pair pT True (Pre Pair Cuts)", kTH1F, {axisPairPtPreCut});
    registry.add("hPairPtReco_PreCut", "Dimuon Pair pT Reco (Pre Pair Cuts)", kTH1F, {axisPairPtPreCut});
    registry.add("hPairPtReco_PostCut", "Dimuon Pair pT Reco (Post Pair Cuts)", kTH1F, {axisPairPtReco});
    registry.add("hSingleMuonPtReco_PostCut", "Single Muon pT Reco (Post Pair Cuts)", kTH1F, {axisPtReco});

    registry.add("hPairPt2ResoVsPt2True", "pT2 Resolution", kTH2F, {axisPt2True, axisPt2RelRes});
  }

  // ---------------------------------------------------------------------------
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
  template <typename TTrack>
  bool passTrackCuts(const TTrack& tr)
  {
    registry.fill(HIST("hCutFlow"), 2); // 2: Track All

    float rAbs = tr.rAtAbsorberEnd();
    if (rAbs < rAbsMin || rAbs > rAbsMax) return false;
    registry.fill(HIST("hCutFlow"), 3); // 3: Track Pass rAbs

    float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
    if (tr.pDca() > pDcaMax) return false;
    registry.fill(HIST("hCutFlow"), 4); // 4: Track Pass pDCA

    if (tr.chi2MatchMCHMID() <= 0) return false;
    registry.fill(HIST("hCutFlow"), 5); // 5: Track Pass MatchMID

    TLorentzVector recoVec;
    recoVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    if (recoVec.Eta() <= etaMin || recoVec.Eta() >= etaMax) return false;
    registry.fill(HIST("hCutFlow"), 6); // 6: Track Pass Eta

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
      if (trkIds.size() < 1) continue;

      registry.fill(HIST("hCutFlow"), 0); // 0: All Cand

      int nMchMid = 0;
      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        if (tr.chi2MatchMCHMID() > 0) nMchMid++;
      }
      if (nMchMid != reqMatchMID) continue;
      
      registry.fill(HIST("hCutFlow"), 1); // 1: Pass Exact MatchMID

      std::vector<int32_t> goodTrkIds;
      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        if (passTrackCuts(tr)) {
          goodTrkIds.push_back(idx);
        }
      }

      // 1. Single Muon Processing (Before Pair Cuts)
      for (auto idx : goodTrkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        const auto& mc = tr.udMcParticle();
        if (std::abs(mc.pdgCode()) == kMuonPDG) {
          TLorentzVector recoVec, trueVec;
          recoVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
          trueVec.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);

          if (trueVec.Pt() > 0.f) {
            float ptTrue = trueVec.Pt(), ptReco = recoVec.Pt();
            float etaTrue = trueVec.Eta(), etaReco = recoVec.Eta();
            float phiTrue = trueVec.Phi(), phiReco = recoVec.Phi();

            registry.fill(HIST("hPtTrue"), ptTrue);
            registry.fill(HIST("hPtReco"), ptReco);
            registry.fill(HIST("hEtaTrue"), etaTrue);
            registry.fill(HIST("hEtaReco"), etaReco);
            registry.fill(HIST("hPhiTrue"), phiTrue);
            registry.fill(HIST("hPhiReco"), phiReco);

            registry.fill(HIST("hPtResoVsPtTrue"), ptTrue, (ptReco - ptTrue) / ptTrue);
            registry.fill(HIST("hResponseMatrixPt"), ptTrue, ptReco);
            registry.fill(HIST("hEtaResoVspTTrue"), ptTrue, etaReco - etaTrue);
            registry.fill(HIST("hResponseMatrixEta"), etaTrue, etaReco);

            float dPhi = phiReco - phiTrue;
            if (dPhi > TMath::Pi()) dPhi -= 2 * TMath::Pi();
            if (dPhi < -TMath::Pi()) dPhi += 2 * TMath::Pi();
            
            registry.fill(HIST("hPhiResoVsPtTrue"), ptTrue, dPhi);
            if (tr.sign() > 0) registry.fill(HIST("hPhiResoVsPtTrue_Pos"), ptTrue, dPhi);
            else if (tr.sign() < 0) registry.fill(HIST("hPhiResoVsPtTrue_Neg"), ptTrue, dPhi);
            
            registry.fill(HIST("hResponseMatrixPhi"), phiTrue, phiReco);
          }
        }
      }

      // 2. Dimuon Processing
      for (size_t i = 0; i < goodTrkIds.size(); ++i) {
        auto tr1 = fwdTracks.iteratorAt(goodTrkIds[i]);
        const auto& mc1 = tr1.udMcParticle();

        for (size_t j = i + 1; j < goodTrkIds.size(); ++j) {
          auto tr2 = fwdTracks.iteratorAt(goodTrkIds[j]);
          const auto& mc2 = tr2.udMcParticle();

          registry.fill(HIST("hCutFlow"), 7); // 7: Pair All

          if (tr1.sign() * tr2.sign() >= 0) continue;
          registry.fill(HIST("hCutFlow"), 8); // 8: Pair Unlike-sign

          if (std::abs(mc1.pdgCode()) != kMuonPDG || std::abs(mc2.pdgCode()) != kMuonPDG) continue;
          registry.fill(HIST("hCutFlow"), 9); // 9: Pair Both MC Muon

          TLorentzVector recoVec1, recoVec2;
          recoVec1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
          recoVec2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
          TLorentzVector pairReco = recoVec1 + recoVec2;

          TLorentzVector mcVec1, mcVec2;
          mcVec1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);
          mcVec2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);
          TLorentzVector pairMC = mcVec1 + mcVec2;

          registry.fill(HIST("hPairPtMC_PreCut"), pairMC.Pt());
          registry.fill(HIST("hPairPtReco_PreCut"), pairReco.Pt());

          // Kinematic Cuts on Pair
          if (pairReco.Pt() >= pairPtMax) continue;
          registry.fill(HIST("hCutFlow"), 10); // 10: Pair Pass Pt

          if (pairReco.Rapidity() <= pairRapidityMin || pairReco.Rapidity() >= pairRapidityMax) continue;
          registry.fill(HIST("hCutFlow"), 11); // 11: Pair Pass Rapidity

          if (pairReco.M() <= pairMassMin || pairReco.M() >= pairMassMax) continue;
          registry.fill(HIST("hCutFlow"), 12); // 12: Pair Pass Mass

          // ===================================================================
          // If passed all cuts, fill Post-Cut Single Muon Resolutions
          // ===================================================================
          registry.fill(HIST("hPairPtReco_PostCut"), pairReco.Pt());
          auto fillPostCutReso = [&](const auto& tr, const auto& mc) {
            TLorentzVector rVec, tVec;
            rVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
            tVec.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);
            if (tVec.Pt() > 0.f) {
              float pT_T = tVec.Pt(), pT_R = rVec.Pt();
              float dEta = rVec.Eta() - tVec.Eta();
              float dPhi = rVec.Phi() - tVec.Phi();
              if (dPhi > TMath::Pi()) dPhi -= 2 * TMath::Pi();
              if (dPhi < -TMath::Pi()) dPhi += 2 * TMath::Pi();

              registry.fill(HIST("hPtResoVsPtTrue_PostCut"), pT_T, (pT_R - pT_T) / pT_T);
              registry.fill(HIST("hEtaResoVspTTrue_PostCut"), pT_T, dEta);
              registry.fill(HIST("hPhiResoVsPtTrue_PostCut"), pT_T, dPhi);
              registry.fill(HIST("hSingleMuonPtReco_PostCut"), pT_R);
            }
          };
          fillPostCutReso(tr1, mc1);
          fillPostCutReso(tr2, mc2);

          // Calculate True pair and fill Mass Resolution
          TLorentzVector trueVec1, trueVec2;
          trueVec1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);
          trueVec2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);
          float massTrue = (trueVec1 + trueVec2).M();
          float massReco = pairReco.M();

          registry.fill(HIST("hMassTrue"), massTrue);
          registry.fill(HIST("hMassReco"), massReco);
          if (massTrue > 0.f) {
            registry.fill(HIST("hMassResoVsMassTrue"), massTrue, (massReco - massTrue) / massTrue);
          }
          registry.fill(HIST("hResponseMatrixMass"), massTrue, massReco);

          float pairPt2True = pairMC.Pt() * pairMC.Pt();
          float pairPt2Reco = pairReco.Pt() * pairReco.Pt();

          if (pairPt2True > 0.f) {
            float pairPt2RelRes = (pairPt2Reco - pairPt2True) / pairPt2True;
            registry.fill(HIST("hPairPt2ResoVsPt2True"), pairPt2True, pairPt2RelRes);
}
        }
      }

      (void)candId; 
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