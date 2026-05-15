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
/// \brief  Resolution analysis for dimuon pairs and single muons in UPC photoproduction.
///         Produces: 1D Pair pT, 1D Pair pT^2, and 2D response matrix.
///         Track type distribution, separated cut-flows, and single/pair resolutions are added.
/// \author Takuma Matsumoto

#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"

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
  // Configurables (設定パラメータ群)
  // ===========================================================================

  // --- Event & Track Selection Level ---
  Configurable<int> reqMatchMFT{"reqMatchMFT", 2, "Required number of valid tracks per event"};
  Configurable<int> reqTrackType{
    "reqTrackType",
    static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack),
    "Track type required for the analysis (0: GlobalMuonTrack, 3: MuonStandaloneTrack)"};

  Configurable<float> etaMin{"etaMin", -4.0f, "Minimum pseudorapidity for muon tracks"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum pseudorapidity for muon tracks"};
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};
  Configurable<float> maxChi2{"maxChi2", 100.f, "Maximum allowed global track Chi2"};
  Configurable<float> maxChi2MatchMCHMFT{"maxChi2MatchMCHMFT", 100.f, "Maximum allowed MCH-MFT match Chi2"};

  // --- Dimuon Pair Kinematic Cuts ---
  Configurable<float> pairPtMax{"pairPtMax", 0.25f, "Maximum dimuon pair pT [GeV/c]"};
  Configurable<float> pairRapidityMin{"pairRapidityMin", -4.0f, "Minimum dimuon pair rapidity"};
  Configurable<float> pairRapidityMax{"pairRapidityMax", -2.5f, "Maximum dimuon pair rapidity"};
  Configurable<float> pairMassMin{"pairMassMin", 1.0f, "Minimum dimuon pair mass [GeV/c^2]"};
  Configurable<float> pairMassMax{"pairMassMax", 10.0f, "Maximum dimuon pair mass [GeV/c^2]"};

  // --- Histogram Binning: Dimuon Pair ---
  Configurable<int> nBinsPt{"nBinsPt", 1000, "Number of bins on pair pT axis"};
  Configurable<float> ptMax{"ptMax", 0.5f, "Upper edge of pair pT axis [GeV/c]"};
  Configurable<int> nBinsPt2{"nBinsPt2", 5000, "Number of bins on pair pT^2 axis"};
  Configurable<float> pt2Max{"pt2Max", 0.5f, "Upper edge of pair pT^2 axis [GeV^2/c^2]"};

  // --- Histogram Binning: Mass ---
  Configurable<int> nBinsMass{"nBinsMass", 300, "Number of bins on mass axis"};
  Configurable<float> massAxisMin{"massAxisMin", 2.0f, "Lower edge of mass axis [GeV/c^2]"};
  Configurable<float> massAxisMax{"massAxisMax", 5.0f, "Upper edge of mass axis [GeV/c^2]"};

  // --- Histogram Binning: Single Track ---
  Configurable<int> nBinsPhi{"nBinsPhi", 360, "Number of bins on phi axis"};
  Configurable<int> nBinsEta{"nBinsEta", 300, "Number of bins on eta axis"};
  Configurable<int> nBinsP{"nBinsP", 1000, "Number of bins on p axis"};
  Configurable<float> pTrackMax{"pTrackMax", 50.0f, "Upper edge of track p axis [GeV/c]"};
  Configurable<float> ptTrackMax{"ptTrackMax", 10.0f, "Upper edge of track pT axis [GeV/c]"};
  Configurable<float> pxpyTrackMax{"pxpyTrackMax", 10.0f, "Upper edge of track px/py axis [GeV/c]"};

  static constexpr int kMuonPDG = 13;
  float mMu = 0.0f;

  // ===========================================================================
  // Initialization
  // ===========================================================================
  void init(InitContext&)
  {
    auto* particle = TDatabasePDG::Instance()->GetParticle(kMuonPDG);
    mMu = particle ? particle->Mass() : 0.105658f;

    initControlHistograms();
    initSingleTrackHistograms();
    initPairHistograms();
  }

  void initControlHistograms()
  {
    auto hCutFlow = registry.add<TH1>("hCutFlow", "Selection Cut Flow;;Counts", HistType::kTH1I, {{18, 0., 18.}});
    TString CutNames[18] = {
      "0: All Cand", "1: Has Requested Track Type", "2: Pass Exact 2 Tracks", "3: Pass Exact MatchMFT",
      "4: Track All", "5: Track Pass Type", "6: Track Pass rAbs", "7: Track Pass pDCA",
      "8: Track Pass MatchMFT", "9: Track Pass Eta", "10: Pair All", "11: Pair Unlike-sign",
      "12: Pair Both MC Muon", "13: Pair Pass Pt", "14: Pair Pass Rapidity", "15: Pair Pass Mass",
      "16: Filler", "17: Filler"};
    for (int i = 0; i < 18; i++) {
      hCutFlow->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
    }

    auto hTrackType = registry.add<TH1>("hTrackType", "Track Type Distribution;Track Type;Counts", HistType::kTH1I, {{5, -0.5, 4.5}});
    hTrackType->GetXaxis()->SetBinLabel(1, "GlobalMuon");     // 0
    hTrackType->GetXaxis()->SetBinLabel(2, "OtherMatch");     // 1
    hTrackType->GetXaxis()->SetBinLabel(3, "GlobalFwd");      // 2
    hTrackType->GetXaxis()->SetBinLabel(4, "MuonStandalone"); // 3
    hTrackType->GetXaxis()->SetBinLabel(5, "MCHStandalone");  // 4

    // Number of tracks per event
    registry.add<TH1>("hNTracksTotal", "Total Forward Tracks per Event;N Tracks;Events", HistType::kTH1I, {{10, -0.5, 10.5}});
    registry.add<TH1>("hNGlobalMuons", "Global Muon Tracks per Event;N Global Muons;Events", HistType::kTH1I, {{20, -0.5, 19.5}});

    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});
  }

  void initSingleTrackHistograms()
  {
    const AxisSpec axPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC} (rad)"};
    const AxisSpec axPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{reco} (rad)"};
    const AxisSpec axPhiRes{200, -0.1f, 0.1f, "#phi^{reco} - #phi^{MC} (rad)"};

    const AxisSpec axEtaMC{nBinsEta, etaMin, etaMax, "#eta^{MC}"};
    const AxisSpec axEtaReco{nBinsEta, etaMin, etaMax, "#eta^{reco}"};
    const AxisSpec axEtaRes{200, -0.1f, 0.1f, "#eta^{reco} - #eta^{MC}"};

    const AxisSpec axPMC{nBinsP, 0.f, pTrackMax, "#it{p}^{MC} (GeV/#it{c})"};
    const AxisSpec axPReco{nBinsP, 0.f, pTrackMax, "#it{p}^{reco} (GeV/#it{c})"};

    const AxisSpec axPtMC{nBinsP, 0.f, ptTrackMax, "#it{p}_{T}^{MC} (GeV/#it{c})"};
    const AxisSpec axPtReco{nBinsP, 0.f, ptTrackMax, "#it{p}_{T}^{reco} (GeV/#it{c})"};

    const AxisSpec axPxMC{nBinsP, -pxpyTrackMax, pxpyTrackMax, "#it{p}_{x}^{MC} (GeV/#it{c})"};
    const AxisSpec axPxReco{nBinsP, -pxpyTrackMax, pxpyTrackMax, "#it{p}_{x}^{reco} (GeV/#it{c})"};

    const AxisSpec axPyMC{nBinsP, -pxpyTrackMax, pxpyTrackMax, "#it{p}_{y}^{MC} (GeV/#it{c})"};
    const AxisSpec axPyReco{nBinsP, -pxpyTrackMax, pxpyTrackMax, "#it{p}_{y}^{reco} (GeV/#it{c})"};

    const AxisSpec axPzMC{nBinsP, -pTrackMax, 0.f, "#it{p}_{z}^{MC} (GeV/#it{c})"};
    const AxisSpec axPzReco{nBinsP, -pTrackMax, 0.f, "#it{p}_{z}^{reco} (GeV/#it{c})"};

    const AxisSpec axPRelRes{200, -0.5f, 0.5f, "(#it{p}^{reco} - #it{p}^{MC}) / #it{p}^{MC}"};

    // 1D distributions
    registry.add("hTrkPhiMC", "Track #phi MC", kTH1F, {axPhiMC});
    registry.add("hTrkPhiReco", "Track #phi Reco", kTH1F, {axPhiReco});
    registry.add("hTrkEtaMC", "Track #eta MC", kTH1F, {axEtaMC});
    registry.add("hTrkEtaReco", "Track #eta Reco", kTH1F, {axEtaReco});

    registry.add("hTrkPMC", "Track #it{p} MC", kTH1F, {axPMC});
    registry.add("hTrkPReco", "Track #it{p} Reco", kTH1F, {axPReco});
    registry.add("hTrkPtMC", "Track #it{p}_{T} MC", kTH1F, {axPtMC});
    registry.add("hTrkPtReco", "Track #it{p}_{T} Reco", kTH1F, {axPtReco});
    registry.add("hTrkPxMC", "Track #it{p}_{x} MC", kTH1F, {axPxMC});
    registry.add("hTrkPxReco", "Track #it{p}_{x} Reco", kTH1F, {axPxReco});
    registry.add("hTrkPyMC", "Track #it{p}_{y} MC", kTH1F, {axPyMC});
    registry.add("hTrkPyReco", "Track #it{p}_{y} Reco", kTH1F, {axPyReco});
    registry.add("hTrkPzMC", "Track #it{p}_{z} MC", kTH1F, {axPzMC});
    registry.add("hTrkPzReco", "Track #it{p}_{z} Reco", kTH1F, {axPzReco});

    // 2D Response / Residuals
    registry.add("hResoPhi", "Track #phi Resolution", kTH2F, {axPhiMC, axPhiRes});
    registry.add("hResoEta", "Track #eta Resolution", kTH2F, {axEtaMC, axEtaRes});
    registry.add("hResoP", "Track #it{p} Resolution", kTH2F, {axPMC, axPRelRes});
    registry.add("hResoPt", "Track #it{p}_{T} Resolution", kTH2F, {axPtMC, axPRelRes});
    registry.add("hResoPx", "Track #it{p}_{x} Resolution", kTH2F, {axPxMC, axPRelRes});
    registry.add("hResoPy", "Track #it{p}_{y} Resolution", kTH2F, {axPyMC, axPRelRes});
    registry.add("hResoPz", "Track #it{p}_{z} Resolution", kTH2F, {axPzMC, axPRelRes});
  }

  // ---------------------------------------------------------------------------
  // ヒストグラム初期化: Dimuonペアレベル
  void initPairHistograms()
  {
    const AxisSpec axisPairPtMC{nBinsPt, 0.f, ptMax, "#it{p}_{T,#mu#mu}^{MC} (GeV/#it{c})"};
    const AxisSpec axisPairPtReco{nBinsPt, 0.f, ptMax, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axisPairPt2MC{nBinsPt2, 0.f, pt2Max, "#it{p}_{T,#mu#mu}^{2,MC} (GeV^{2}/#it{c}^{2})"};
    const AxisSpec axisPairPt2Reco{nBinsPt2, 0.f, pt2Max, "#it{p}_{T,#mu#mu}^{2,reco} (GeV^{2}/#it{c}^{2})"};

    const AxisSpec axPairRapMC{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}_{#mu#mu}^{MC}"};
    const AxisSpec axPairRapReco{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}_{#mu#mu}^{reco}"};
    const AxisSpec axPairRapRes{200, -0.1f, 0.1f, "#it{y}_{#mu#mu}^{reco} - #it{y}_{#mu#mu}^{MC}"};

    const AxisSpec axPairEtaMC{nBinsEta, -10.0f, 0.0f, "#eta_{#mu#mu}^{MC}"};
    const AxisSpec axPairEtaReco{nBinsEta, -10.0f, 0.0f, "#eta_{#mu#mu}^{reco}"};
    const AxisSpec axPairEtaRes{200, -2.0f, 2.0f, "#eta_{#mu#mu}^{reco} - #eta_{#mu#mu}^{MC}"};

    const AxisSpec axPairPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi_{#mu#mu}^{MC} (rad)"};
    const AxisSpec axPairPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi_{#mu#mu}^{reco} (rad)"};
    const AxisSpec axPairPhiRes{200, -0.1f, 0.1f, "#phi_{#mu#mu}^{reco} - #phi_{#mu#mu}^{MC} (rad)"};

    const AxisSpec axPairMassMC{nBinsMass, massAxisMin, massAxisMax, "#it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassReco{nBinsMass, massAxisMin, massAxisMax, "#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassRes{200, -0.5f, 0.5f, "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};

    // 1D Histograms
    registry.add("hPairPtMC_PreCut", "Dimuon Pair #it{p}_{T} MC (Pre Kin. Cuts)", kTH1F, {axisPairPtMC});
    registry.add("hPairPtReco_PreCut", "Dimuon Pair #it{p}_{T} Reco (Pre Kin. Cuts)", kTH1F, {axisPairPtReco});
    registry.add("hPairPtMC_PostCut", "Dimuon Pair #it{p}_{T} MC (Post All Cuts)", kTH1F, {axisPairPtMC});
    registry.add("hPairPtReco_PostCut", "Dimuon Pair #it{p}_{T} Reco (Post All Cuts)", kTH1F, {axisPairPtReco});

    registry.add("hPairPt2MC_PreCut", "Dimuon Pair #it{p}_{T}^{2} MC (Pre Kin. Cuts)", kTH1F, {axisPairPt2MC});
    registry.add("hPairPt2Reco_PreCut", "Dimuon Pair #it{p}_{T}^{2} Reco (Pre Kin. Cuts)", kTH1F, {axisPairPt2Reco});
    registry.add("hPairPt2MC_PostCut", "Dimuon Pair #it{p}_{T}^{2} MC (Post All Cuts)", kTH1F, {axisPairPt2MC});
    registry.add("hPairPt2Reco_PostCut", "Dimuon Pair #it{p}_{T}^{2} Reco (Post All Cuts)", kTH1F, {axisPairPt2Reco});

    registry.add("hPairRapMC_PostCut", "Dimuon Pair Rapidity MC (Post All Cuts)", kTH1F, {axPairRapMC});
    registry.add("hPairRapReco_PostCut", "Dimuon Pair Rapidity Reco (Post All Cuts)", kTH1F, {axPairRapReco});
    registry.add("hPairEtaMC_PostCut", "Dimuon Pair #eta MC (Post All Cuts)", kTH1F, {axPairEtaMC});
    registry.add("hPairEtaReco_PostCut", "Dimuon Pair #eta Reco (Post All Cuts)", kTH1F, {axPairEtaReco});
    registry.add("hPairPhiMC_PostCut", "Dimuon Pair #phi MC (Post All Cuts)", kTH1F, {axPairPhiMC});
    registry.add("hPairPhiReco_PostCut", "Dimuon Pair #phi Reco (Post All Cuts)", kTH1F, {axPairPhiReco});

    registry.add("hPairMassMC_PostCut", "Dimuon Pair Mass MC (Post All Cuts)", kTH1F, {axPairMassMC});
    registry.add("hPairMassReco_PostCut", "Dimuon Pair Mass Reco (Post All Cuts)", kTH1F, {axPairMassReco});

    // 2D Response Matrices
    registry.add("hResponseMatrixPairPt", "Pair #it{p}_{T} Response Matrix", kTH2F, {axisPairPtMC, axisPairPtReco});
    registry.add("hResponseMatrixPairPt2", "Pair #it{p}_{T}^{2} Response Matrix", kTH2F, {axisPairPt2MC, axisPairPt2Reco});
    registry.add("hResponseMatrixPairRap", "Pair Rapidity Response Matrix", kTH2F, {axPairRapMC, axPairRapReco});
    registry.add("hResponseMatrixPairEta", "Pair #eta Response Matrix", kTH2F, {axPairEtaMC, axPairEtaReco});
    registry.add("hResponseMatrixPairPhi", "Pair #phi Response Matrix", kTH2F, {axPairPhiMC, axPairPhiReco});
    registry.add("hResponseMatrixPairMass", "Pair Mass Response Matrix", kTH2F, {axPairMassMC, axPairMassReco});

    // 2D Residuals (Resolutions)
    const AxisSpec axisPtRelRes{200, -1.0f, 1.0f, "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};
    const AxisSpec axisPt2RelRes{200, -5.0f, 30.0f, "(#it{p}_{T}^{2,reco} - #it{p}_{T}^{2,MC}) / #it{p}_{T}^{2,MC}"};
    registry.add("hPairPtResoVsPtMC", "Pair #it{p}_{T} Resolution vs #it{p}_{T}^{MC}", kTH2F, {axisPairPtMC, axisPtRelRes});
    registry.add("hPairPt2ResoVsPt2MC", "Pair #it{p}_{T}^{2} Resolution vs #it{p}_{T}^{2,MC}", kTH2F, {axisPairPt2MC, axisPt2RelRes});
    registry.add("hPairRapResoVsRapMC", "Pair Rapidity Resolution vs #it{y}^{MC}", kTH2F, {axPairRapMC, axPairRapRes});
    registry.add("hPairEtaResoVsEtaMC", "Pair #eta Resolution vs #eta^{MC}", kTH2F, {axPairEtaMC, axPairEtaRes});
    registry.add("hPairPhiResoVsPhiMC", "Pair #phi Resolution vs #phi^{MC}", kTH2F, {axPairPhiMC, axPairPhiRes});
    registry.add("hPairMassResoVsMassMC", "Pair Mass Resolution vs Mass MC", kTH2F, {axPairMassMC, axPairMassRes});
  }

  // ===========================================================================
  // Core Methods
  // ===========================================================================

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
    registry.fill(HIST("hCutFlow"), 4); // 4: Track All

    float rAbs = tr.rAtAbsorberEnd();
    if (rAbs < rAbsMin || rAbs > rAbsMax)
      return false;
    registry.fill(HIST("hCutFlow"), 6); // 6: Track Pass rAbs

    float pDcaMax = (rAbs < 26.5f) ? 350.0f : 200.0f;
    if (tr.pDca() > pDcaMax)
      return false;
    registry.fill(HIST("hCutFlow"), 7); // 7: Track Pass pDCA

    if (tr.chi2() > maxChi2)
      return false;

    if (reqTrackType == static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      float chi2MFT = tr.chi2MatchMCHMFT();
      if (chi2MFT < 0.f || chi2MFT > maxChi2MatchMCHMFT) {
        return false;
      }
    }
    registry.fill(HIST("hCutFlow"), 8); // 8: Track Pass MatchMFT

    TLorentzVector recoVec;
    recoVec.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    if (recoVec.Eta() <= etaMin || recoVec.Eta() >= etaMax)
      return false;
    registry.fill(HIST("hCutFlow"), 9); // 9: Track Pass Eta

    return true;
  }

  // ---------------------------------------------------------------------------
  template <typename TTrack, typename TMcParticle>
  void fillSingleTrackAnalysis(const TTrack& tr, const TMcParticle& mc)
  {
    TLorentzVector vReco, vMC;
    vReco.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    vMC.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);

    registry.fill(HIST("hTrkPhiMC"), vMC.Phi());
    registry.fill(HIST("hTrkPhiReco"), vReco.Phi());
    registry.fill(HIST("hTrkEtaMC"), vMC.Eta());
    registry.fill(HIST("hTrkEtaReco"), vReco.Eta());

    registry.fill(HIST("hTrkPMC"), vMC.P());
    registry.fill(HIST("hTrkPReco"), vReco.P());
    registry.fill(HIST("hTrkPtMC"), vMC.Pt());
    registry.fill(HIST("hTrkPtReco"), vReco.Pt());
    registry.fill(HIST("hTrkPxMC"), vMC.Px());
    registry.fill(HIST("hTrkPxReco"), vReco.Px());
    registry.fill(HIST("hTrkPyMC"), vMC.Py());
    registry.fill(HIST("hTrkPyReco"), vReco.Py());
    registry.fill(HIST("hTrkPzMC"), vMC.Pz());
    registry.fill(HIST("hTrkPzReco"), vReco.Pz());

    registry.fill(HIST("hResoPhi"), vMC.Phi(), TVector2::Phi_mpi_pi(vReco.Phi() - vMC.Phi()));
    registry.fill(HIST("hResoEta"), vMC.Eta(), vReco.Eta() - vMC.Eta());

    if (vMC.P() > 0.f)
      registry.fill(HIST("hResoP"), vMC.P(), (vReco.P() - vMC.P()) / vMC.P());
    if (vMC.Pt() > 0.f)
      registry.fill(HIST("hResoPt"), vMC.Pt(), (vReco.Pt() - vMC.Pt()) / vMC.Pt());
    if (std::abs(vMC.Px()) > 0.f)
      registry.fill(HIST("hResoPx"), vMC.Px(), (vReco.Px() - vMC.Px()) / std::abs(vMC.Px()));
    if (std::abs(vMC.Py()) > 0.f)
      registry.fill(HIST("hResoPy"), vMC.Py(), (vReco.Py() - vMC.Py()) / std::abs(vMC.Py()));
    if (std::abs(vMC.Pz()) > 0.f)
      registry.fill(HIST("hResoPz"), vMC.Pz(), (vReco.Pz() - vMC.Pz()) / std::abs(vMC.Pz()));
  }

  // ---------------------------------------------------------------------------
  template <typename TTrack, typename TMcParticle>
  void processPairAnalysis(const TTrack& tr1, const TMcParticle& mc1,
                           const TTrack& tr2, const TMcParticle& mc2)
  {
    registry.fill(HIST("hCutFlow"), 10); // 10: Pair All

    if (tr1.sign() * tr2.sign() >= 0)
      return;
    registry.fill(HIST("hCutFlow"), 11); // 11: Pair Unlike-sign

    if (std::abs(mc1.pdgCode()) != kMuonPDG || std::abs(mc2.pdgCode()) != kMuonPDG)
      return;
    registry.fill(HIST("hCutFlow"), 12); // 12: Pair Both MC Muon

    TLorentzVector recoVec1, recoVec2, mcVec1, mcVec2;
    recoVec1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    recoVec2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    mcVec1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);
    mcVec2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);

    TLorentzVector pairReco = recoVec1 + recoVec2;
    TLorentzVector pairMC = mcVec1 + mcVec2;

    float pairPtMC = pairMC.Pt();
    float pairPtReco = pairReco.Pt();
    float pairPt2MC = pairPtMC * pairPtMC;
    float pairPt2Reco = pairPtReco * pairPtReco;
    float pairRapMC = pairMC.Rapidity();
    float pairRapReco = pairReco.Rapidity();
    float pairEtaMC = pairMC.Eta();
    float pairEtaReco = pairReco.Eta();
    float pairPhiMC = pairMC.Phi();
    float pairPhiReco = pairReco.Phi();
    float pairMassMC = pairMC.M();
    float pairMassReco = pairReco.M();

    registry.fill(HIST("hPairPtMC_PreCut"), pairPtMC);
    registry.fill(HIST("hPairPtReco_PreCut"), pairPtReco);
    registry.fill(HIST("hPairPt2MC_PreCut"), pairPt2MC);
    registry.fill(HIST("hPairPt2Reco_PreCut"), pairPt2Reco);

    if (pairReco.Pt() >= pairPtMax)
      return;
    registry.fill(HIST("hCutFlow"), 13); // 13: Pair Pass Pt

    if (pairReco.Rapidity() <= pairRapidityMin || pairReco.Rapidity() >= pairRapidityMax)
      return;
    registry.fill(HIST("hCutFlow"), 14); // 14: Pair Pass Rapidity

    if (pairReco.M() <= pairMassMin || pairReco.M() >= pairMassMax)
      return;
    registry.fill(HIST("hCutFlow"), 15); // 15: Pair Pass Mass

    registry.fill(HIST("hPairPtMC_PostCut"), pairPtMC);
    registry.fill(HIST("hPairPtReco_PostCut"), pairPtReco);
    registry.fill(HIST("hPairPt2MC_PostCut"), pairPt2MC);
    registry.fill(HIST("hPairPt2Reco_PostCut"), pairPt2Reco);

    registry.fill(HIST("hPairRapMC_PostCut"), pairRapMC);
    registry.fill(HIST("hPairRapReco_PostCut"), pairRapReco);
    registry.fill(HIST("hResponseMatrixPairRap"), pairRapMC, pairRapReco);
    registry.fill(HIST("hPairRapResoVsRapMC"), pairRapMC, pairRapReco - pairRapMC);

    registry.fill(HIST("hPairEtaMC_PostCut"), pairEtaMC);
    registry.fill(HIST("hPairEtaReco_PostCut"), pairEtaReco);
    registry.fill(HIST("hResponseMatrixPairEta"), pairEtaMC, pairEtaReco);
    registry.fill(HIST("hPairEtaResoVsEtaMC"), pairEtaMC, pairEtaReco - pairEtaMC);

    registry.fill(HIST("hPairPhiMC_PostCut"), pairPhiMC);
    registry.fill(HIST("hPairPhiReco_PostCut"), pairPhiReco);
    registry.fill(HIST("hResponseMatrixPairPhi"), pairPhiMC, pairPhiReco);
    registry.fill(HIST("hPairPhiResoVsPhiMC"), pairPhiMC, TVector2::Phi_mpi_pi(pairPhiReco - pairPhiMC));

    registry.fill(HIST("hPairMassMC_PostCut"), pairMassMC);
    registry.fill(HIST("hPairMassReco_PostCut"), pairMassReco);
    registry.fill(HIST("hResponseMatrixPairMass"), pairMassMC, pairMassReco);
    registry.fill(HIST("hPairMassResoVsMassMC"), pairMassMC, pairMassReco - pairMassMC);

    registry.fill(HIST("hResponseMatrixPairPt"), pairPtMC, pairPtReco);
    if (pairPtMC > 0.f) {
      registry.fill(HIST("hPairPtResoVsPtMC"), pairPtMC, (pairPtReco - pairPtMC) / pairPtMC);
    }

    registry.fill(HIST("hResponseMatrixPairPt2"), pairPt2MC, pairPt2Reco);
    if (pairPt2MC > 0.f) {
      registry.fill(HIST("hPairPt2ResoVsPt2MC"), pairPt2MC, (pairPt2Reco - pairPt2MC) / pairPt2MC);
    }
  }

  // ===========================================================================
  // Main Processing Loop
  // ===========================================================================
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

      // --- NEW: イベント内の全前方トラック数をカウント ---
      registry.fill(HIST("hNTracksTotal"), static_cast<float>(trkIds.size()));

      if (trkIds.empty())
        continue;

      registry.fill(HIST("hCutFlow"), 0); // 0: All Cand

      std::vector<int32_t> candidateTrkIds;
      candidateTrkIds.reserve(trkIds.size());

      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        registry.fill(HIST("hTrackType"), static_cast<float>(tr.trackType()));

        if (tr.trackType() == reqTrackType) {
          candidateTrkIds.push_back(idx);
        }
      }

      // --- NEW: イベント内の GlobalMuon トラック数をカウント ---
      registry.fill(HIST("hNGlobalMuons"), static_cast<float>(candidateTrkIds.size()));

      if (candidateTrkIds.empty())
        continue;
      registry.fill(HIST("hCutFlow"), 1); // 1: Has Requested Track Type

      if (candidateTrkIds.size() != reqMatchMFT)
        continue;
      registry.fill(HIST("hCutFlow"), 2); // 2: Pass Exact 2 Tracks

      std::vector<int32_t> goodTrkIds;
      goodTrkIds.reserve(reqMatchMFT);

      for (auto idx : candidateTrkIds) {
        auto tr = fwdTracks.iteratorAt(idx);

        if (passTrackCuts(tr)) {
          goodTrkIds.push_back(idx);
          const auto& mcParticle = tr.udMcParticle();
          if (std::abs(mcParticle.pdgCode()) == kMuonPDG) {
            fillSingleTrackAnalysis(tr, mcParticle);
          }
        }
      }

      if (goodTrkIds.size() != reqMatchMFT)
        continue;
      registry.fill(HIST("hCutFlow"), 3); // 3: Pass Exact MatchMFT

      for (size_t i = 0; i < goodTrkIds.size(); ++i) {
        auto tr1 = fwdTracks.iteratorAt(goodTrkIds[i]);
        const auto& mc1 = tr1.udMcParticle();

        for (size_t j = i + 1; j < goodTrkIds.size(); ++j) {
          auto tr2 = fwdTracks.iteratorAt(goodTrkIds[j]);
          const auto& mc2 = tr2.udMcParticle();

          processPairAnalysis(tr1, mc1, tr2, mc2);
        }
      }

      (void)candId;
      (void)eventCandidates;
    }
  }

  PROCESS_SWITCH(UPCMuonPairResolution, processMcReco, "Fill single and pair resolution histograms from MC reco data", true);
};

// ============================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCMuonPairResolution>(cfgc, TaskName{"my-upc-muon-pair-resolution"}),
  };
}
