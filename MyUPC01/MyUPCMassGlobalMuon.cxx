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

/// \file   MyUPCMassGlobalMuon.cxx
/// \brief  Combined data/MC dimuon mass analysis for GlobalMuon (MFT-MCH-MID) tracks.
///         - Data/MC switch: PROCESS_SWITCH (processDataNoZdc / processDataWithZdc / processMcReco / processMcTruth)
///         - ZDC usage switch: PROCESS_SWITCH (processDataWithZdc), since the
///           O2udzdcreduce table is only subscribed to when that process is enabled.
///           A runtime Configurable cannot be used for this: the table must be
///           absent from the function signature entirely when the input file
///           does not contain it (e.g. some MC productions).
///         - Mandatory histograms: invariant mass, pair pT, pair pT^2
///           (exactly one mu+ and one mu- per candidate, hence unlike-sign only by construction)
///         - All other histograms are QA for checking the effect of each cut
/// \author Takuma
/// \date 2026

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr int kMuonPDG = 13;
constexpr float kMaxZDCTime = 2.0f; // [ns]
}

struct MyUPCMassGlobalMuonTask {

  Service<o2::framework::O2DatabasePDG> pdg;

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ===========================================================================
  // Configurables
  // ===========================================================================

  // --- Optional-detector switches ---
  Configurable<bool> useV0ACut{"useV0ACut", true, "Apply V0A same-BC amplitude veto"};
  Configurable<float> cutV0AAmp{"cutV0AAmp", 100.0f, "Max V0A amplitude in same BC [a.u.]"};

  // --- Track type (this task targets GlobalMuon tracks) ---
  Configurable<int> reqTrackType{"reqTrackType",
                                  static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack),
                                  "Required forward-track type (0: GlobalMuonTrack, 3: MuonStandaloneTrack)"};

  // --- Single-track quality cuts ---
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};
  Configurable<float> maxChi2{"maxChi2", 100.f, "Maximum global-track chi2"};
  Configurable<float> maxChi2MatchMCHMFT{"maxChi2MatchMCHMFT", 100.f, "Maximum MCH-MFT match chi2 (GlobalMuon only)"};

  // --- Single-track kinematic acceptance ---
  Configurable<float> etaMin{"etaMin", -3.5f, "Minimum single-muon pseudorapidity"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum single-muon pseudorapidity"};
  Configurable<float> pTmin{"pTmin", 0.3f, "Minimum single-muon pT [GeV/c]"};

  // --- Dimuon pair kinematic cuts ---
  Configurable<float> pairPtMax{"pairPtMax", 0.25f, "Maximum dimuon pair pT [GeV/c]"};
  Configurable<float> pairRapidityMin{"pairRapidityMin", -3.5f, "Minimum dimuon pair rapidity"};
  Configurable<float> pairRapidityMax{"pairRapidityMax", -2.5f, "Maximum dimuon pair rapidity"};
  Configurable<float> pairMassMin{"pairMassMin", 1.0f, "Minimum dimuon pair mass [GeV/c^2]"};
  Configurable<float> pairMassMax{"pairMassMax", 10.0f, "Maximum dimuon pair mass [GeV/c^2]"};

  // --- ZDC neutron topology (processData, only used if useZDC == true) ---
  Configurable<float> cutZNAEnergy{"cutZNAEnergy", 1.0f, "ZNA energy threshold for neutron tag [TeV]"};
  Configurable<float> cutZNCEnergy{"cutZNCEnergy", 1.0f, "ZNC energy threshold for neutron tag [TeV]"};
  Configurable<int> targetTopology{"targetTopology", -1, "Required neutron topology (0=0n0n,1=Xn0n,2=0nXn,3=XnXn, -1=No Cut)"};

  // --- Mandatory histogram binning: mass / pT / pT^2 ---
  Configurable<int> nBinsMass{"nBinsMass", 900, "Number of bins on the invariant-mass axis"};
  Configurable<float> massAxisMin{"massAxisMin", 0.0f, "Lower edge of the invariant-mass axis [GeV/c^2]"};
  Configurable<float> massAxisMax{"massAxisMax", 10.0f, "Upper edge of the invariant-mass axis [GeV/c^2]"};

  Configurable<int> nBinsPt{"nBinsPt", 500, "Number of bins on the pair pT axis"};
  Configurable<float> ptAxisMax{"ptAxisMax", 1.0f, "Upper edge of the pair pT axis [GeV/c]"};

  Configurable<int> nBinsPt2{"nBinsPt2", 500, "Number of bins on the pair pT^2 axis"};
  Configurable<float> pt2AxisMax{"pt2AxisMax", 1.0f, "Upper edge of the pair pT^2 axis [GeV^2/c^2]"};

  // --- QA histogram binning (cut-check histograms) ---
  Configurable<int> nBinsPhi{"nBinsPhi", 120, "Number of bins on phi axes (QA)"};
  Configurable<int> nBinsEtaQA{"nBinsEtaQA", 80, "Number of bins on eta/rapidity axes (QA)"};
  Configurable<int> nBinsPtMuon{"nBinsPtMuon", 250, "Number of bins on single-muon pT axis (QA)"};
  Configurable<float> ptMuonAxisMax{"ptMuonAxisMax", 5.0f, "Upper edge of single-muon pT axis (QA) [GeV/c]"};

  float mMu = 0.0f;

  // ===========================================================================
  // ZDC helper
  // ===========================================================================
  struct ZDCinfo {
    float timeA = -999.f;
    float timeC = -999.f;
    float enA = -999.f;
    float enC = -999.f;
    int znClass = -1; // 0=0n0n, 1=Xn0n, 2=0nXn, 3=XnXn, -1=unknown
  };

  // ===========================================================================
  // Initialization
  // ===========================================================================
  void init(InitContext&)
  {
    mMu = pdg->Mass(kMuonPDG);

    // --- Mandatory axes: invariant mass, pair pT, pair pT^2 ---
    const AxisSpec axisMass{nBinsMass, massAxisMin, massAxisMax, "#it{M}_{#mu#mu} (GeV/#it{c}^{2})"};
    const AxisSpec axisPt{nBinsPt, 0., ptAxisMax, "#it{p}_{T,#mu#mu} (GeV/#it{c})"};
    const AxisSpec axisPt2{nBinsPt2, 0., pt2AxisMax, "#it{p}_{T,#mu#mu}^{2} (GeV^{2}/#it{c}^{2})"};

    // --- Mandatory histograms (unlike-sign pairs only, by construction) ---
    registry.add("hMassUnlike", "Invariant mass, unlike-sign pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPtUnlike", "Pair #it{p}_{T}, unlike-sign pairs;;#counts", kTH1D, {axisPt});
    registry.add("hPt2Unlike", "Pair #it{p}_{T}^{2}, unlike-sign pairs;;#counts", kTH1D, {axisPt2});

    // --- MC-truth counterparts (filled by processMcTruth only) ---
    registry.add("hMassTruth", "Invariant mass, MC truth unlike-sign pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPtTruth", "Pair #it{p}_{T}, MC truth unlike-sign pairs;;#counts", kTH1D, {axisPt});
    registry.add("hPt2Truth", "Pair #it{p}_{T}^{2}, MC truth unlike-sign pairs;;#counts", kTH1D, {axisPt2});

    // --- Event-level QA ---
    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("eventCounter", "Processed events", kTH1F, {axisCounter});
    registry.add("hNumContrib", "Number of contributors;N_{contrib};#counts", kTH1D, {{100, -0.5, 99.5}});
    registry.add("hNTracksTotal", "Forward tracks per candidate;N_{tracks};#counts", kTH1I, {{20, -0.5, 19.5}});
    registry.add("hNGoodTracks", "GlobalMuon tracks passing quality cuts per candidate;N_{tracks};#counts", kTH1I, {{20, -0.5, 19.5}});

    // --- Cut-flow ---
    auto hCutFlow = registry.add<TH1>("hCutFlow", "Cut flow;;#counts", HistType::kTH1I, {{13, 0., 13.}});
    TString cutNames[13] = {
      "AllCand", "V0A_pass", "ZDCTopology",
      "Track_All", "Track_Type", "Track_Quality",
      "HasTwoGoodTracks", "Pair_EtaCut",
      "Pair_PtCutMuon", "Pair_PtCut", "Pair_RapidityCut",
      "Pair_MassCut", "Filled"};
    for (int i = 0; i < 13; ++i) {
      hCutFlow->GetXaxis()->SetBinLabel(i + 1, cutNames[i].Data());
    }

    // --- Track-type distribution (QA) ---
    auto hTrackType = registry.add<TH1>("hTrackType", "Track type distribution;Track type;#counts", HistType::kTH1I, {{5, -0.5, 4.5}});
    hTrackType->GetXaxis()->SetBinLabel(1, "GlobalMuon");
    hTrackType->GetXaxis()->SetBinLabel(2, "OtherMatch");
    hTrackType->GetXaxis()->SetBinLabel(3, "GlobalFwd");
    hTrackType->GetXaxis()->SetBinLabel(4, "MuonStandalone");
    hTrackType->GetXaxis()->SetBinLabel(5, "MCHStandalone");

    // --- Single-track quality QA (filled before the cut, for cut-tuning) ---
    registry.add("hRAbs", "R at absorber end;R_{abs} (cm);#counts", kTH1F, {{120, 0., 120.}});
    registry.add("hPDCA", "p#times DCA;p#times DCA (GeV/#it{c}#upointcm);#counts", kTH1F, {{200, 0., 1000.}});
    registry.add("hChi2", "Global track #chi^{2};#chi^{2};#counts", kTH1F, {{200, 0., 200.}});
    registry.add("hChi2MatchMCHMFT", "MCH-MFT match #chi^{2};#chi^{2};#counts", kTH1F, {{220, -10., 200.}});

    // --- Single-muon kinematics QA ---
    const AxisSpec axisEtaMuon{nBinsEtaQA, -5.0, 0.0, "#eta_{#mu}"};
    const AxisSpec axisPhi{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi (rad)"};
    const AxisSpec axisPtMuon{nBinsPtMuon, 0., ptMuonAxisMax, "#it{p}_{T,#mu}} (GeV/#it{c})"};
    registry.add("hEtaMuon", "Single muon #eta;;#counts", kTH1F, {axisEtaMuon});
    registry.add("hPhiMuon", "Single muon #phi;;#counts", kTH1F, {axisPhi});
    registry.add("hPtMuon", "Single muon #it{p}_{T};;#counts", kTH1F, {axisPtMuon});

    // --- Pair kinematics QA ---
    const AxisSpec axisRapidityPair{nBinsEtaQA, pairRapidityMin, pairRapidityMax, "#it{y}_{#mu#mu}"};
    registry.add("hRapidityPair", "Pair rapidity;;#counts", kTH1F, {axisRapidityPair});
    registry.add("hPhiPair", "Pair #phi;;#counts", kTH1F, {axisPhi});

    // --- V0A QA ---
    registry.add("hV0AAmp", "Max V0A amplitude in same BC;Amplitude (a.u.);#counts", kTH1D, {{500, 0., 500.}});

    // --- ZDC QA (only if processDataWithZdc is enabled) ---
    if (doprocessDataWithZdc) {
      const AxisSpec axisZDCEnergy{1000, -2.5, 199.5, "E_{ZN} (TeV)"};
      const AxisSpec axisTopology{4, -0.5, 3.5, "0n0n, Xn0n, 0nXn, XnXn"};
      registry.add("hEnergyZNA", "Common ZNA energy;;#counts", kTH1F, {axisZDCEnergy});
      registry.add("hEnergyZNC", "Common ZNC energy;;#counts", kTH1F, {axisZDCEnergy});
      registry.add("hEnergyZNAvsZNC", "ZNA vs ZNC energy;E_{ZNC};E_{ZNA}", kTH2D, {axisZDCEnergy, axisZDCEnergy});
      auto hTopologyCounter = registry.add<TH1>("hTopologyCounter", "Neutron topology;;#counts", HistType::kTH1I, {{4, 0., 4.}});
      TString topologyClasses[4] = {"0n0n", "Xn0n", "0nXn", "XnXn"};
      for (int i = 0; i < 4; ++i) {
        hTopologyCounter->GetXaxis()->SetBinLabel(i + 1, topologyClasses[i].Data());
      }
    }
  }

  // ===========================================================================
  // Helpers
  // ===========================================================================

  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks const& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  // V0A same-BC amplitude veto. Always fills the QA histogram; the cut itself (Table producer do same cut)
  // is applied only if useV0ACut == true.
  template <typename TCand>
  bool passV0ACut(TCand const& cand)
  {
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();

    float maxAmpV0A = 0.0f;
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) == 0) {
        maxAmpV0A = std::max(maxAmpV0A, ampsV0A[i]);
      }
    }
    registry.fill(HIST("hV0AAmp"), maxAmpV0A);

    return maxAmpV0A < cutV0AAmp;
  }

  // Builds a per-candidate ZDC info map from the ZDC table (data only).
  std::unordered_map<int32_t, ZDCinfo> buildZdcMap(o2::aod::UDZdcsReduced const& zdcs)
  {
    std::unordered_map<int32_t, ZDCinfo> zdcMap;
    for (const auto& zdc : zdcs) {
      int32_t candId = zdc.udCollisionId();
      if (candId < 0) {
        continue;
      }

      ZDCinfo info;
      info.timeA = zdc.timeZNA();
      info.timeC = zdc.timeZNC();
      info.enA = zdc.energyCommonZNA();
      info.enC = zdc.energyCommonZNC();

      bool hasTimeA = (!std::isinf(info.timeA) && std::abs(info.timeA) < kMaxZDCTime);
      bool hasTimeC = (!std::isinf(info.timeC) && std::abs(info.timeC) < kMaxZDCTime);
      bool isNeutronA = hasTimeA && (info.enA > cutZNAEnergy);
      bool isNeutronC = hasTimeC && (info.enC > cutZNCEnergy);

      info.znClass = 0; // 0n0n
      if (isNeutronA && !isNeutronC) {
        info.znClass = 1; // Xn0n
      } else if (!isNeutronA && isNeutronC) {
        info.znClass = 2; // 0nXn
      } else if (isNeutronA && isNeutronC) {
        info.znClass = 3; // XnXn
      }

      zdcMap[candId] = info;
    }
    return zdcMap;
  }

  // Selects tracks of the requested track type passing the single-track
  // quality cuts (rAbs, pDCA, chi2, chi2MatchMCHMFT). Fills the corresponding
  // QA / cut-flow histograms along the way.
  template <typename TTracks>
  std::vector<int32_t> selectGoodTracks(std::vector<int32_t> const& trkIds, TTracks const& tracks)
  {
    std::vector<int32_t> good;
    good.reserve(trkIds.size());

    for (auto idx : trkIds) {
      auto tr = tracks.iteratorAt(idx);
      registry.fill(HIST("hCutFlow"), 3); // Track_All
      registry.fill(HIST("hTrackType"), static_cast<float>(tr.trackType()));

      if (tr.trackType() != reqTrackType) {
        continue;
      }
      registry.fill(HIST("hCutFlow"), 4); // Track_Type

      float rAbs = tr.rAtAbsorberEnd();
      registry.fill(HIST("hRAbs"), rAbs);
      if (rAbs < rAbsMin || rAbs > rAbsMax) {
        continue;
      }

      float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
      registry.fill(HIST("hPDCA"), tr.pDca());
      if (tr.pDca() > pDcaMax) {
        continue;
      }

      registry.fill(HIST("hChi2"), tr.chi2());
      if (tr.chi2() > maxChi2) {
        continue;
      }

      if (reqTrackType == static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
        float chi2MFT = tr.chi2MatchMCHMFT();
        registry.fill(HIST("hChi2MatchMCHMFT"), chi2MFT);
        if (chi2MFT < 0.f || chi2MFT > maxChi2MatchMCHMFT) {
          continue;
        }
      }

      registry.fill(HIST("hCutFlow"), 5); // Track_Quality
      good.push_back(idx);
    }
    return good;
  }

  // Fills the mandatory mass/pT/pT^2 histograms plus the pair-level QA
  // histograms, applying the single-muon and pair kinematic cuts.
  // Shared between processData and processMcReco. The pair is guaranteed to
  // be unlike-sign (one mu+ and one mu-) by the HasTwoGoodTracks selection
  // in processCandidates.
  template <typename TTrack>
  void fillPairHistograms(TTrack const& tr1, TTrack const& tr2)
  {
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);

    registry.fill(HIST("hEtaMuon"), p1.Eta());
    registry.fill(HIST("hEtaMuon"), p2.Eta());
    registry.fill(HIST("hPhiMuon"), p1.Phi());
    registry.fill(HIST("hPhiMuon"), p2.Phi());
    registry.fill(HIST("hPtMuon"), p1.Pt());
    registry.fill(HIST("hPtMuon"), p2.Pt());

    if (p1.Eta() <= etaMin || p1.Eta() >= etaMax || p2.Eta() <= etaMin || p2.Eta() >= etaMax) {
      return;
    }
    registry.fill(HIST("hCutFlow"), 7); // Pair_EtaCut

    if (p1.Pt() < pTmin || p2.Pt() < pTmin) {
      return;
    }
    registry.fill(HIST("hCutFlow"), 8); // Pair_PtCutMuon

    TLorentzVector p = p1 + p2;
    float pt2 = p.Pt() * p.Pt();

    registry.fill(HIST("hRapidityPair"), p.Rapidity());
    registry.fill(HIST("hPhiPair"), p.Phi());

    if (p.Pt() >= pairPtMax) {
      return;
    }
    registry.fill(HIST("hCutFlow"), 9); // Pair_PtCut

    if (p.Rapidity() <= pairRapidityMin || p.Rapidity() >= pairRapidityMax) {
      return;
    }
    registry.fill(HIST("hCutFlow"), 10); // Pair_RapidityCut

    if (p.M() <= pairMassMin || p.M() >= pairMassMax) {
      return;
    }
    registry.fill(HIST("hCutFlow"), 11); // Pair_MassCut

    registry.fill(HIST("hCutFlow"), 12); // Filled
    registry.fill(HIST("hMassUnlike"), p.M());
    registry.fill(HIST("hPtUnlike"), p.Pt());
    registry.fill(HIST("hPt2Unlike"), pt2);
  }

  // Common candidate loop, shared between processData and processMcReco.
  // `zdcMap`/`applyZdc` are no-ops unless useZDC == true (and applyZdc == true).
  template <typename TTracks>
  void processCandidates(CandidatesFwd const& eventCandidates, TTracks const& tracks,
                          std::unordered_map<int32_t, ZDCinfo> const& zdcMap, bool applyZdc)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, tracks);

    for (const auto& item : tracksPerCand) {
      if (item.second.size() < 2) {
        continue; // need at least 2 tracks to form a pair
      }

      registry.fill(HIST("hCutFlow"), 0); // AllCand
      registry.fill(HIST("hNTracksTotal"), static_cast<float>(item.second.size()));

      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      registry.fill(HIST("hNumContrib"), cand.numContrib());

      if (useV0ACut && !passV0ACut(cand)) {
        continue;
      }
      registry.fill(HIST("hCutFlow"), 1); // V0A_pass

      if (applyZdc) {
        ZDCinfo zdc;
        auto it = zdcMap.find(candID);
        if (it != zdcMap.end()) {
          zdc = it->second;
        }
        registry.fill(HIST("hEnergyZNA"), zdc.enA);
        registry.fill(HIST("hEnergyZNC"), zdc.enC);
        registry.fill(HIST("hEnergyZNAvsZNC"), zdc.enC, zdc.enA);
        if (zdc.znClass >= 0) {
          registry.fill(HIST("hTopologyCounter"), zdc.znClass);
        }
        if (targetTopology > -1 && zdc.znClass != targetTopology) {
          continue;
        }
      }
      registry.fill(HIST("hCutFlow"), 2); // ZDCTopology

      auto goodTracks = selectGoodTracks(item.second, tracks);
      registry.fill(HIST("hNGoodTracks"), static_cast<float>(goodTracks.size()));
      if (goodTracks.size() != 2) {
        continue;
      }

      auto tr1 = tracks.iteratorAt(goodTracks[0]);
      auto tr2 = tracks.iteratorAt(goodTracks[1]);
      if (tr1.sign() + tr2.sign() != 0) {
        continue; // require exactly one mu+ and one mu-
      }
      registry.fill(HIST("hCutFlow"), 6); // HasTwoGoodTracks

      fillPairHistograms(tr1, tr2);
    }
  }

  // ===========================================================================
  // Process: real data, without ZDC (O2udzdcreduce table not subscribed)
  // ===========================================================================
  void processDataNoZdc(CandidatesFwd const& eventCandidates,
                         ForwardTracks const& fwdTracks)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    processCandidates(eventCandidates, fwdTracks, std::unordered_map<int32_t, ZDCinfo>{}, false);
  }
  PROCESS_SWITCH(MyUPCMassGlobalMuonTask, processDataNoZdc, "Process real data without ZDC information", true);

  // ===========================================================================
  // Process: real data, with ZDC (O2udzdcreduce table subscribed)
  // ===========================================================================
  void processDataWithZdc(CandidatesFwd const& eventCandidates,
                           ForwardTracks const& fwdTracks,
                           o2::aod::UDZdcsReduced const& zdcs)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    std::unordered_map<int32_t, ZDCinfo> zdcMap = buildZdcMap(zdcs);
    processCandidates(eventCandidates, fwdTracks, zdcMap, true);
  }
  PROCESS_SWITCH(MyUPCMassGlobalMuonTask, processDataWithZdc, "Process real data with ZDC neutron-topology selection/QA", false);

  // ===========================================================================
  // Process: MC reco (no ZDC table required)
  // ===========================================================================
  void processMcReco(CandidatesFwd const& eventCandidates,
                      CompleteFwdTracks const& fwdTracks,
                      o2::aod::UDMcCollisions const&,
                      o2::aod::UDMcParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    processCandidates(eventCandidates, fwdTracks, std::unordered_map<int32_t, ZDCinfo>{}, false);
  }
  PROCESS_SWITCH(MyUPCMassGlobalMuonTask, processMcReco, "Process MC reco tracks (ZDC not used)", false);

  // ===========================================================================
  // Process: MC truth (fills hMassTruth / hPtTruth / hPt2Truth)
  // ===========================================================================
  void processMcTruth(o2::aod::UDMcCollisions const& /*mcCollisions*/,
                       o2::aod::UDMcParticles const& mcParticles)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> muonsPerMcColl;

    for (const auto& mc : mcParticles) {
      if (std::abs(mc.pdgCode()) != kMuonPDG) {
        continue;
      }

      TLorentzVector v;
      v.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);
      if (v.Eta() <= etaMin || v.Eta() >= etaMax) {
        continue;
      }
      if (v.Pt() < pTmin) {
        continue;
      }

      muonsPerMcColl[mc.udMcCollisionId()].push_back(mc.globalIndex());
    }

    for (const auto& item : muonsPerMcColl) {
      const auto& ids = item.second;
      for (size_t i = 0; i < ids.size(); ++i) {
        auto mc1 = mcParticles.iteratorAt(ids[i]);
        TLorentzVector v1;
        v1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);

        for (size_t j = i + 1; j < ids.size(); ++j) {
          auto mc2 = mcParticles.iteratorAt(ids[j]);
          if (mc1.pdgCode() * mc2.pdgCode() >= 0) {
            continue; // require unlike-sign muon pair
          }

          TLorentzVector v2;
          v2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);
          TLorentzVector p = v1 + v2;

          if (p.Pt() >= pairPtMax) {
            continue;
          }
          if (p.Rapidity() <= pairRapidityMin || p.Rapidity() >= pairRapidityMax) {
            continue;
          }
          if (p.M() <= pairMassMin || p.M() >= pairMassMax) {
            continue;
          }

          registry.fill(HIST("hMassTruth"), p.M());
          registry.fill(HIST("hPtTruth"), p.Pt());
          registry.fill(HIST("hPt2Truth"), p.Pt() * p.Pt());
        }
      }
    }
  }
  PROCESS_SWITCH(MyUPCMassGlobalMuonTask, processMcTruth, "Fill MC-truth mass/pT/pT^2 histograms", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyUPCMassGlobalMuonTask>(cfgc, TaskName{"my-upc-mass-global-muon"})};
}
