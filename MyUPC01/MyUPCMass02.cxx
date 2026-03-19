// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file MyUPCMass02.cxx
/// \brief Custom task to compute invariant mass of dimuon at forward rapidity for UPC events.
/// Integrates exact kinematic cuts from MyUPCMass01, tree/MC structures from FwdMuonsUPC,
/// and adds event counter with > 1.0 GeV mass cut.

#include "PWGUD/DataModel/UDTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <unordered_map>
#include <vector>

// =========================================================================
// Table declarations for saving tree with info on data
// =========================================================================
namespace dimu
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(PhiAv, phiAv, float);
DECLARE_SOA_COLUMN(PhiCh, phiCh, float);
DECLARE_SOA_COLUMN(EnergyP, energyP, float);
DECLARE_SOA_COLUMN(Pxp, pxp, float);
DECLARE_SOA_COLUMN(Pyp, pyp, float);
DECLARE_SOA_COLUMN(Pzp, pzp, float);
DECLARE_SOA_COLUMN(Ptp, ptp, float);
DECLARE_SOA_COLUMN(Etap, etap, float);
DECLARE_SOA_COLUMN(Phip, phip, float);
DECLARE_SOA_COLUMN(TrackTypep, trackTypep, int);
DECLARE_SOA_COLUMN(EnergyN, energyN, float);
DECLARE_SOA_COLUMN(Pxn, pxn, float);
DECLARE_SOA_COLUMN(Pyn, pyn, float);
DECLARE_SOA_COLUMN(Pzn, pzn, float);
DECLARE_SOA_COLUMN(Ptn, ptn, float);
DECLARE_SOA_COLUMN(Etan, etan, float);
DECLARE_SOA_COLUMN(Phin, phin, float);
DECLARE_SOA_COLUMN(TrackTypen, trackTypen, int);
DECLARE_SOA_COLUMN(Tzna, tzna, float);
DECLARE_SOA_COLUMN(Ezna, ezna, float);
DECLARE_SOA_COLUMN(Tznc, tznc, float);
DECLARE_SOA_COLUMN(Eznc, eznc, float);
DECLARE_SOA_COLUMN(Nclass, nclass, int);
} // namespace dimu

namespace o2::aod
{
DECLARE_SOA_TABLE(DiMu, "AOD", "DIMU",
                  dimu::RunNumber, dimu::M, dimu::Energy, dimu::Px, dimu::Py, dimu::Pz, dimu::Pt, dimu::Rap, dimu::Phi,
                  dimu::PhiAv, dimu::PhiCh,
                  dimu::EnergyP, dimu::Pxp, dimu::Pyp, dimu::Pzp, dimu::Ptp, dimu::Etap, dimu::Phip, dimu::TrackTypep,
                  dimu::EnergyN, dimu::Pxn, dimu::Pyn, dimu::Pzn, dimu::Ptn, dimu::Etan, dimu::Phin, dimu::TrackTypen,
                  dimu::Tzna, dimu::Ezna, dimu::Tznc, dimu::Eznc, dimu::Nclass);
} // namespace o2::aod

namespace gendimu
{
DECLARE_SOA_COLUMN(GenM, genM, float);
DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenRap, genRap, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenPhiAv, genPhiAv, float);
DECLARE_SOA_COLUMN(GenPhiCh, genPhiCh, float);
DECLARE_SOA_COLUMN(GenPtp, genPtp, float);
DECLARE_SOA_COLUMN(GenEtap, genEtap, float);
DECLARE_SOA_COLUMN(GenPhip, genPhip, float);
DECLARE_SOA_COLUMN(GenPtn, genPtn, float);
DECLARE_SOA_COLUMN(GenEtan, genEtan, float);
DECLARE_SOA_COLUMN(GenPhin, genPhin, float);
} // namespace gendimu

namespace o2::aod
{
DECLARE_SOA_TABLE(GenDimu, "AOD", "GENDIMU",
                  gendimu::GenM, gendimu::GenPt, gendimu::GenRap, gendimu::GenPhi,
                  gendimu::GenPhiAv, gendimu::GenPhiCh,
                  gendimu::GenPtp, gendimu::GenEtap, gendimu::GenPhip,
                  gendimu::GenPtn, gendimu::GenEtan, gendimu::GenPhin);
} // namespace o2::aod

namespace recodimu
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(PhiAv, phiAv, float);
DECLARE_SOA_COLUMN(PhiCh, phiCh, float);
DECLARE_SOA_COLUMN(Ptp, ptp, float);
DECLARE_SOA_COLUMN(Etap, etap, float);
DECLARE_SOA_COLUMN(Phip, phip, float);
DECLARE_SOA_COLUMN(TrackTypep, trackTypep, int);
DECLARE_SOA_COLUMN(Ptn, ptn, float);
DECLARE_SOA_COLUMN(Etan, etan, float);
DECLARE_SOA_COLUMN(Phin, phin, float);
DECLARE_SOA_COLUMN(TrackTypen, trackTypen, int);
DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenRap, genRap, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenPtp, genPtp, float);
DECLARE_SOA_COLUMN(GenEtap, genEtap, float);
DECLARE_SOA_COLUMN(GenPhip, genPhip, float);
DECLARE_SOA_COLUMN(GenPtn, genPtn, float);
DECLARE_SOA_COLUMN(GenEtan, genEtan, float);
DECLARE_SOA_COLUMN(GenPhin, genPhin, float);
} // namespace recodimu

namespace o2::aod
{
DECLARE_SOA_TABLE(RecoDimu, "AOD", "RECODIMU",
                  recodimu::RunNumber, recodimu::M, recodimu::Pt, recodimu::Rap, recodimu::Phi,
                  recodimu::PhiAv, recodimu::PhiCh,
                  recodimu::Ptp, recodimu::Etap, recodimu::Phip, recodimu::TrackTypep,
                  recodimu::Ptn, recodimu::Etan, recodimu::Phin, recodimu::TrackTypen,
                  recodimu::GenPt, recodimu::GenRap, recodimu::GenPhi,
                  recodimu::GenPtp, recodimu::GenEtap, recodimu::GenPhip,
                  recodimu::GenPtn, recodimu::GenEtan, recodimu::GenPhin);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Constants
const float kEtaMinGeneral = -4.0;
const float kEtaMaxGeneral = -2.5;
const int kReqMatchMIDTracks = 2;
const int kReqMatchMFTTracks = 2;
const int kMaxChi2MFTMatch = 30;
const float kMaxZDCTime = 2.;
const int kMuonPDG = 13;

struct MyUPCMass02 {

  Service<o2::framework::O2DatabasePDG> pdg;

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  Produces<o2::aod::DiMu> dimuSel;
  Produces<o2::aod::GenDimu> dimuGen;
  Produces<o2::aod::RecoDimu> dimuReco;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcGenRegistry{"mcGenRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcRecoRegistry{"mcRecoRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // --- Configurables ---
  static constexpr double Pi = o2::constants::math::PI;
  Configurable<int> nBinsPt{"nBinsPt", 250, "N bins in pT histo"};
  Configurable<float> lowPt{"lowPt", 0.0f, "lower limit in pT histo [GeV/c]"};
  Configurable<float> highPt{"highPt", 0.25f, "upper limit in pT histo [GeV/c]"};
  Configurable<int> nBinsMass{"nBinsMass", 9000, "N bins in InvMass histo (10 MeV/c^2 per bin)"};
  Configurable<float> lowMass{"lowMass", 1.0f, "lower limit in mass histo [GeV/c^2]"};
  Configurable<float> highMass{"highMass", 10.0f, "upper limit in mass histo [GeV/c^2]"};
  Configurable<int> nBinsRapidity{"nBinsRapidity", 150, "N bins in rapidity histo"};
  Configurable<float> lowRapidity{"lowRapidity", -4.0f, "lower limit in rapidity histo"};
  Configurable<float> highRapidity{"highRapidity", -2.5f, "upper limit in rapidity histo"};

  Configurable<int> myTrackType{"myTrackType", 3, "My track type (0 = MFT, 1 = MCH-MID)"};
  Configurable<float> cutZNAEnergy{"cutZNAEnergy", 1.0f, "ZNA energy threshold [TeV]"};
  Configurable<float> cutZNCEnergy{"cutZNCEnergy", 1.0f, "ZNC energy threshold [TeV]"};
  Configurable<int> targetTopology{"targetTopology", -1, "Required neutron topology (1=0n0n, 2=Xn0n, 3=0nXn, 4=XnXn, -1=No Cut)"};

  void init(InitContext&)
  {
    const AxisSpec axisPt{nBinsPt, lowPt, highPt, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisMass{nBinsMass, lowMass, highMass, "m_{#mu#mu} GeV/#it{c}^{2}"};
    const AxisSpec axisRapidity{nBinsRapidity, lowRapidity, highRapidity, "Rapidity"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisV0A{500, 0.0, 500.0, "Max V0A Amplitude (Same BC) [a.u.]"};

    // Add event counter
    registry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});
    mcGenRegistry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});
    mcRecoRegistry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "Selection Cut Flow;;Candidates", HistType::kTH1I, {{15, 0., 15.}});
    TString SelectionCuts[15] = {"AllPairs", "V0A_pass", "TopologySelected", "TrackEnd", "pDCA", "MatchMID", "MatchMFT", "EtaCut", "PairPtCut", "PairRapCut", "MassCut", "UnlikeSign", "LikeSign", "SavedToTree", "Filler"};
    for (int i = 0; i < 15; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    // QA Histogram for V0A
    registry.add("hV0A_Amp", "Maximum V0A Amplitude in Same BC;;#counts", kTH1D, {axisV0A});

    registry.add("hMassUnlike", "Invariant mass of Unlike-sign pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPtUnlike", "Transverse momentum of Unlike-sign pairs;;#counts", kTH1D, {axisPt});
    registry.add("hMassLike", "Invariant mass of Like-sign pairs;;#counts", kTH1D, {axisMass});
    registry.add("hPtLike", "Transverse momentum of Like-sign pairs;;#counts", kTH1D, {axisPt});
    registry.add("hRapidity", "Rapidty of muon pairs;;#counts", kTH1D, {axisRapidity});

    registry.add("hMassVsRapidityUnlike", "Invariant mass vs Rapidity of Unlike-sign pairs;Rapidity;m_{#mu#mu} GeV/#it{c}^{2}", kTH2D, {axisRapidity, axisMass});

    mcGenRegistry.add("hMass", "Invariant mass of muon pairs;;#counts", kTH1D, {axisMass});
    mcRecoRegistry.add("hMass", "Invariant mass of muon pairs;;#counts", kTH1D, {axisMass});
  }

  float particleMass(int pid) { return pdg->Mass(pid); }

  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId >= 0)
        tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  template <typename TTracks>
  void collectMcCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udMcCollisionId();
      if (candId >= 0)
        tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  template <typename TTracks>
  void collectRecoCandID(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId >= 0 && tr.has_udMcParticle()) {
        tracksPerCand[candId].push_back(tr.globalIndex());
        tracksPerCand[candId].push_back(tr.udMcParticle().globalIndex());
      }
    }
  }

  struct ZDCinfo {
    float timeA;
    float timeC;
    float enA;
    float enC;
    int32_t znClass;
  };

  void computePhiAnis(TLorentzVector p1, TLorentzVector p2, int sign1, float& phiAverage, float& phiCharge)
  {
    TLorentzVector tSum, tDiffAv, tDiffCh;
    tSum = p1 + p2;
    float halfUnity = 0.5;
    if (sign1 > 0) {
      tDiffCh = p1 - p2;
      tDiffAv = (gRandom->Rndm() > halfUnity) ? (p1 - p2) : (p2 - p1);
    } else {
      tDiffCh = p2 - p1;
      tDiffAv = (gRandom->Rndm() > halfUnity) ? (p2 - p1) : (p1 - p2);
    }
    phiAverage = tSum.DeltaPhi(tDiffAv);
    phiCharge = tSum.DeltaPhi(tDiffCh);
  }

  template <typename TTrack>
  bool checkPDCA(const TTrack& tr)
  {
    float rAbs = tr.rAtAbsorberEnd();
    float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
    return tr.pDca() <= pDcaMax;
  }

  // --- Real Data Processing Logic ---
  void processCandidate(CandidatesFwd::iterator const& cand,
                        ForwardTracks::iterator const& tr1, ForwardTracks::iterator const& tr2,
                        ZDCinfo const& zdc)
  {
    registry.fill(HIST("hSelectionCounter"), 0); // 0: AllPairs

    // V0A Amplitude cut: A(FV0) < 100 a.u. in the same BC & QA Histogram
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    bool skipV0A = false;
    float maxAmpV0A = 0.0f;
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) == 0) { // Same BC
        if (ampsV0A[i] > maxAmpV0A) {
          maxAmpV0A = ampsV0A[i];
        }
        if (ampsV0A[i] >= 100.0f) {
          skipV0A = true;
        }
      }
    }

    // Fill QA histogram before cutting
    registry.fill(HIST("hV0A_Amp"), maxAmpV0A);

    if (skipV0A)
      return;

    registry.fill(HIST("hSelectionCounter"), 1); // 1: V0A_pass

    if (tr1.rAtAbsorberEnd() < 17.6f || tr1.rAtAbsorberEnd() > 89.5f ||
        tr2.rAtAbsorberEnd() < 17.6f || tr2.rAtAbsorberEnd() > 89.5f)
      return;
    registry.fill(HIST("hSelectionCounter"), 2); // 2: TrackEnd

    if (targetTopology > -1 && zdc.znClass != targetTopology)
      return;
    registry.fill(HIST("hSelectionCounter"), 3); // 3: TopologySelected

    if (!checkPDCA(*tr1) || !checkPDCA(*tr2))
      return;
    registry.fill(HIST("hSelectionCounter"), 4); // 4: pDCA

    if (tr1.chi2MatchMCHMID() <= 0 || tr2.chi2MatchMCHMID() <= 0)
      return;
    registry.fill(HIST("hSelectionCounter"), 5); // 5: MatchMID

    if (myTrackType == 0) {
      if (tr1.chi2MatchMCHMFT() <= 0 || tr1.chi2MatchMCHMFT() >= kMaxChi2MFTMatch ||
          tr2.chi2MatchMCHMFT() <= 0 || tr2.chi2MatchMCHMFT() >= kMaxChi2MFTMatch)
        return;
    }
    registry.fill(HIST("hSelectionCounter"), 6); // 6: MatchMFT

    TLorentzVector p1, p2;
    auto mMu = particleMass(kMuonPDG);
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;

    if (p1.Eta() <= kEtaMinGeneral || p1.Eta() >= kEtaMaxGeneral ||
        p2.Eta() <= kEtaMinGeneral || p2.Eta() >= kEtaMaxGeneral)
      return;
    registry.fill(HIST("hSelectionCounter"), 7); // 7: EtaCut

    if (p.Pt() >= highPt)
      return;
    registry.fill(HIST("hSelectionCounter"), 8); // 8: PairPtCut

    if (p.Rapidity() <= lowRapidity || p.Rapidity() >= highRapidity)
      return;
    registry.fill(HIST("hSelectionCounter"), 9); // 9: PairRapCut

    if (p.M() <= lowMass || p.M() >= highMass)
      return;
    registry.fill(HIST("hSelectionCounter"), 10); // 10: MassCut

    bool isUnlikeSign = (tr1.sign() + tr2.sign()) == 0;
    if (isUnlikeSign) {
      registry.fill(HIST("hSelectionCounter"), 11); // 11: UnlikeSign
      registry.fill(HIST("hMassUnlike"), p.M());
      registry.fill(HIST("hPtUnlike"), p.Pt());
      registry.fill(HIST("hRapidity"), p.Rapidity());
      registry.fill(HIST("hMassVsRapidityUnlike"), p.Rapidity(), p.M());
    } else {
      registry.fill(HIST("hSelectionCounter"), 12); // 12: LikeSign
      registry.fill(HIST("hMassLike"), p.M());
      registry.fill(HIST("hPtLike"), p.Pt());
    }

    if (isUnlikeSign) {
      registry.fill(HIST("hSelectionCounter"), 13); // 13: SavedToTree
      float phiAverage = 0, phiCharge = 0;
      computePhiAnis(p1, p2, tr1.sign(), phiAverage, phiCharge);

      auto& pPos = (tr1.sign() > 0) ? p1 : p2;
      auto& pNeg = (tr1.sign() > 0) ? p2 : p1;

      dimuSel(cand.runNumber(),
              p.M(), p.E(), p.Px(), p.Py(), p.Pz(), p.Pt(), p.Rapidity(), p.Phi(),
              phiAverage, phiCharge,
              pPos.E(), pPos.Px(), pPos.Py(), pPos.Pz(), pPos.Pt(), pPos.PseudoRapidity(), pPos.Phi(), static_cast<int>(myTrackType),
              pNeg.E(), pNeg.Px(), pNeg.Py(), pNeg.Pz(), pNeg.Pt(), pNeg.PseudoRapidity(), pNeg.Phi(), static_cast<int>(myTrackType),
              zdc.timeA, zdc.enA, zdc.timeC, zdc.enC, zdc.znClass);
    }
  }

  void processData(CandidatesFwd const& eventCandidates,
                   o2::aod::UDZdcsReduced const& ZDCs,
                   ForwardTracks const& fwdTracks)
  {
    // Fill the event counter for each processed event block
    registry.fill(HIST("eventCounter"), 0.5);

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    std::unordered_map<int32_t, ZDCinfo> zdcPerCand;
    for (const auto& zdc : ZDCs) {
      int32_t candId = zdc.udCollisionId();
      if (candId < 0)
        continue;

      float tA = zdc.timeZNA();
      float tC = zdc.timeZNC();
      float eA = zdc.energyCommonZNA();
      float eC = zdc.energyCommonZNC();

      bool hasTimeA = (!std::isinf(tA) && std::abs(tA) < kMaxZDCTime);
      bool hasTimeC = (!std::isinf(tC) && std::abs(tC) < kMaxZDCTime);
      bool hasEnergyA = (eA > cutZNAEnergy);
      bool hasEnergyC = (eC > cutZNCEnergy);

      bool isNeutronA = hasTimeA && hasEnergyA;
      bool isNeutronC = hasTimeC && hasEnergyC;

      int zClass = 1; // 1 = 0n0n
      if (isNeutronA && !isNeutronC)
        zClass = 2; // 2 = Xn0n
      else if (!isNeutronA && isNeutronC)
        zClass = 3; // 3 = 0nXn
      else if (isNeutronA && isNeutronC)
        zClass = 4; // 4 = XnXn

      zdcPerCand[candId] = {tA, tC, eA, eC, zClass};
    }

    for (const auto& item : tracksPerCand) {
      if (item.second.size() < 2)
        continue;

      // exactly two MCH-MID tracks cut
      int nMatchMID = 0;
      for (size_t k = 0; k < item.second.size(); ++k) {
        auto tr = fwdTracks.iteratorAt(item.second[k]);
        if (tr.chi2MatchMCHMID() > 0) {
          nMatchMID++;
        }
      }
      if (nMatchMID != 2)
        continue;

      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      ZDCinfo zdc = zdcPerCand.count(candID) ? zdcPerCand.at(candID) : ZDCinfo{-999, -999, -999, -999, -1};

      for (size_t i = 0; i < item.second.size() - 1; ++i) {
        for (size_t j = i + 1; j < item.second.size(); ++j) {
          auto tr1 = fwdTracks.iteratorAt(item.second[i]);
          auto tr2 = fwdTracks.iteratorAt(item.second[j]);

          if (tr1.chi2MatchMCHMID() > 0 && tr2.chi2MatchMCHMID() > 0) {
            processCandidate(cand, tr1, tr2, zdc);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(MyUPCMass02, processData, "", true);

  // --- MC Truth Logic ---
  void processMcGenCand(aod::UDMcCollisions::iterator const& /*mcCand*/,
                        aod::UDMcParticles::iterator const& McPart1, aod::UDMcParticles::iterator const& McPart2)
  {
    if (std::abs(McPart1.pdgCode()) != kMuonPDG || std::abs(McPart2.pdgCode()) != kMuonPDG)
      return;
    if (McPart1.pdgCode() + McPart2.pdgCode() != 0)
      return;

    TLorentzVector p1, p2;
    auto mMu = particleMass(kMuonPDG);
    p1.SetXYZM(McPart1.px(), McPart1.py(), McPart1.pz(), mMu);
    p2.SetXYZM(McPart2.px(), McPart2.py(), McPart2.pz(), mMu);
    TLorentzVector p = p1 + p2;

    if (p.M() <= lowMass || p.M() >= highMass)
      return;
    if (p.Pt() >= highPt)
      return;
    if (p.Rapidity() <= lowRapidity || p.Rapidity() >= highRapidity)
      return;

    float phiAverage = 0, phiCharge = 0;
    computePhiAnis(p1, p2, -McPart1.pdgCode(), phiAverage, phiCharge);
    mcGenRegistry.fill(HIST("hMass"), p.M());

    auto& pPos = (McPart1.pdgCode() < 0) ? p1 : p2;
    auto& pNeg = (McPart1.pdgCode() < 0) ? p2 : p1;

    dimuGen(p.M(), p.Pt(), p.Rapidity(), p.Phi(),
            phiAverage, phiCharge,
            pPos.Pt(), pPos.PseudoRapidity(), pPos.Phi(),
            pNeg.Pt(), pNeg.PseudoRapidity(), pNeg.Phi());
  }

  void processMcGen(aod::UDMcCollisions const& mccollisions, aod::UDMcParticles const& McParts)
  {
    mcGenRegistry.fill(HIST("eventCounter"), 0.5);

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectMcCandIDs(tracksPerCand, McParts);

    for (const auto& item : tracksPerCand) {
      if (item.second.size() < 2)
        continue;
      auto cand = mccollisions.iteratorAt(item.first);

      for (size_t i = 0; i < item.second.size() - 1; ++i) {
        for (size_t j = i + 1; j < item.second.size(); ++j) {
          auto tr1 = McParts.iteratorAt(item.second[i]);
          auto tr2 = McParts.iteratorAt(item.second[j]);
          processMcGenCand(cand, tr1, tr2);
        }
      }
    }
  }
  PROCESS_SWITCH(MyUPCMass02, processMcGen, "", false);

  // --- MC Reco Logic ---
  void processMcRecoCand(CandidatesFwd::iterator const& cand,
                         CompleteFwdTracks::iterator const& tr1, aod::UDMcParticles::iterator const& McPart1,
                         CompleteFwdTracks::iterator const& tr2, aod::UDMcParticles::iterator const& McPart2)
  {
    if ((tr1.sign() + tr2.sign()) != 0)
      return;

    TLorentzVector p1, p2;
    auto mMu = particleMass(kMuonPDG);
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    TLorentzVector p = p1 + p2;

    if (p.M() <= lowMass || p.M() >= highMass)
      return;
    if (p.Pt() >= highPt)
      return;
    if (p.Rapidity() <= lowRapidity || p.Rapidity() >= highRapidity)
      return;

    TLorentzVector p1Mc, p2Mc;
    p1Mc.SetXYZM(McPart1.px(), McPart1.py(), McPart1.pz(), mMu);
    p2Mc.SetXYZM(McPart2.px(), McPart2.py(), McPart2.pz(), mMu);
    TLorentzVector pMc = p1Mc + p2Mc;

    float phiAverage = 0, phiCharge = 0;
    computePhiAnis(p1, p2, tr1.sign(), phiAverage, phiCharge);
    mcRecoRegistry.fill(HIST("hMass"), p.M());

    auto& pPos = (tr1.sign() > 0) ? p1 : p2;
    auto& pNeg = (tr1.sign() > 0) ? p2 : p1;
    auto& pPosMc = (tr1.sign() > 0) ? p1Mc : p2Mc;
    auto& pNegMc = (tr1.sign() > 0) ? p2Mc : p1Mc;

    dimuReco(cand.runNumber(),
             p.M(), p.Pt(), p.Rapidity(), p.Phi(), phiAverage, phiCharge,
             pPos.Pt(), pPos.PseudoRapidity(), pPos.Phi(), static_cast<int>(myTrackType),
             pNeg.Pt(), pNeg.PseudoRapidity(), pNeg.Phi(), static_cast<int>(myTrackType),
             pMc.Pt(), pMc.Rapidity(), pMc.Phi(),
             pPosMc.Pt(), pPosMc.PseudoRapidity(), pPosMc.Phi(),
             pNegMc.Pt(), pNegMc.PseudoRapidity(), pNegMc.Phi());
  }

  void processMcReco(CandidatesFwd const& eventCandidates,
                     CompleteFwdTracks const& fwdTracks,
                     aod::UDMcCollisions const&,
                     aod::UDMcParticles const& McParts)
  {
    mcRecoRegistry.fill(HIST("eventCounter"), 0.5);

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCandAll;
    collectRecoCandID(tracksPerCandAll, fwdTracks);

    for (const auto& item : tracksPerCandAll) {
      size_t nPairsData = item.second.size() / 2;
      if (nPairsData < 2)
        continue;

      // exactly two MCH-MID tracks cut
      int nMatchMID = 0;
      for (size_t k = 0; k < nPairsData; ++k) {
        auto tr = fwdTracks.iteratorAt(item.second[k * 2]);
        if (tr.chi2MatchMCHMID() > 0) {
          nMatchMID++;
        }
      }
      if (nMatchMID != 2)
        continue;

      auto cand = eventCandidates.iteratorAt(item.first);

      for (size_t i = 0; i < nPairsData - 1; ++i) {
        for (size_t j = i + 1; j < nPairsData; ++j) {
          auto tr1 = fwdTracks.iteratorAt(item.second[i * 2]);
          auto trMc1 = McParts.iteratorAt(item.second[i * 2 + 1]);
          auto tr2 = fwdTracks.iteratorAt(item.second[j * 2]);
          auto trMc2 = McParts.iteratorAt(item.second[j * 2 + 1]);

          if (tr1.chi2MatchMCHMID() > 0 && tr2.chi2MatchMCHMID() > 0) {
            processMcRecoCand(cand, tr1, trMc1, tr2, trMc2);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(MyUPCMass02, processMcReco, "", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyUPCMass02>(cfgc, TaskName{"my-upc-mass-02"}),
  };
}
