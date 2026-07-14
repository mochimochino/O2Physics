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

/// \file   UPCMuonAnalysisMC.cxx
/// \brief  Unified resolution and efficiency analysis for dimuon pairs and single muons
///         in UPC photoproduction.
///         Produces:
///           - Resolution: 1D Pair pT, pT^2; 2D response matrices and residuals for
///             single tracks and dimuon pairs.
///           - Efficiency: 1D MC-truth and reco histograms (phi, eta, pT) for single
///             tracks and dimuon pairs; ratio (eff) histograms.
///           - High-pT tail diagnostics (run separately on Coherent/Incoherent MC
///             samples and compared offline): pair Delta-pT (reco-MC) vs reco mass,
///             |Delta-phi| between the two reco muons vs pair reco pT, and single-track
///             eta/phi residual vs MC p (MFT-MCH matching bias).
///         The MC-reco processing loop is shared between both analyses.
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

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ============================================================================

struct UPCMuonAnalysisMC {

  Service<o2::framework::O2DatabasePDG> pdg;

  // Track type aliases
  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;
  using CompleteFwdTracks = soa::Join<ForwardTracks, o2::aod::UDMcFwdTrackLabels>;

  // Histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ===========================================================================
  // Configurables
  // ===========================================================================

  // --- Track type ---
  Configurable<int> reqMatchMFT{"reqMatchMFT", 2, "Required number of valid tracks per event"};
  Configurable<int> reqTrackType{
    "reqTrackType",
    static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack),
    "Track type required for the analysis (0: GlobalMuonTrack, 3: MuonStandaloneTrack)"};

  // --- Single track acceptance (applied to both MC truth and reco) ---
  Configurable<float> etaMin{"etaMin", -3.6f, "Minimum pseudorapidity for muon tracks"};
  Configurable<float> etaMax{"etaMax", -2.5f, "Maximum pseudorapidity for muon tracks"};
  Configurable<float> pTmin{"pTmin", 0.3f, "Minimum single muon pT [GeV/c]"};

  // --- Reco quality cuts ---
  Configurable<float> rAbsMin{"rAbsMin", 17.6f, "Minimum R at absorber end [cm]"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5f, "Maximum R at absorber end [cm]"};
  Configurable<float> maxChi2{"maxChi2", 100.f, "Maximum allowed global track Chi2"};
  Configurable<float> maxChi2MatchMCHMFT{"maxChi2MatchMCHMFT", 100.f, "Maximum allowed MCH-MFT match Chi2"};

  // --- Dimuon pair kinematic cuts (applied to both MC truth and reco) ---
  Configurable<float> pairPtMax{"pairPtMax", 0.25f, "Maximum dimuon pair pT [GeV/c]"};
  Configurable<float> pairRapidityMin{"pairRapidityMin", -3.6f, "Minimum dimuon pair rapidity"};
  Configurable<float> pairRapidityMax{"pairRapidityMax", -2.5f, "Maximum dimuon pair rapidity"};
  Configurable<float> pairMassMin{"pairMassMin", 1.0f, "Minimum dimuon pair mass [GeV/c^2]"};
  Configurable<float> pairMassMax{"pairMassMax", 10.0f, "Maximum dimuon pair mass [GeV/c^2]"};

  // --- J/psi rapidity window for Acc x Eff denominator ---
  Configurable<float> jpsiRapMin{"jpsiRapMin", -4.5f,
                                 "Min J/psi rapidity for Acc x Eff denominator"};
  Configurable<float> jpsiRapMax{"jpsiRapMax", -2.0f,
                                 "Max J/psi rapidity for Acc x Eff denominator"};

  // --- Histogram binning: pair pT / pT^2 ---
  Configurable<int> nBinsPt{"nBinsPt", 1000, "Number of bins on pair pT axis"};
  Configurable<float> ptMax{"ptMax", 0.5f, "Upper edge of pair pT axis [GeV/c]"};
  Configurable<int> nBinsPt2{"nBinsPt2", 5000, "Number of bins on pair pT^2 axis"};
  Configurable<float> pt2Max{"pt2Max", 0.5f, "Upper edge of pair pT^2 axis [GeV^2/c^2]"};

  // --- Histogram binning: mass ---
  Configurable<int> nBinsMass{"nBinsMass", 300, "Number of bins on mass axis"};
  Configurable<float> massAxisMin{"massAxisMin", 2.0f, "Lower edge of mass axis [GeV/c^2]"};
  Configurable<float> massAxisMax{"massAxisMax", 5.0f, "Upper edge of mass axis [GeV/c^2]"};

  // --- Histogram binning: single track ---
  Configurable<int> nBinsPhi{"nBinsPhi", 360, "Number of bins on phi axis"};
  Configurable<int> nBinsEta{"nBinsEta", 300, "Number of bins on eta/rapidity axis"};
  Configurable<int> nBinsP{"nBinsP", 1000, "Number of bins on p axis"};
  Configurable<float> pTrackMax{"pTrackMax", 50.0f, "Upper edge of track p axis [GeV/c]"};
  Configurable<float> ptTrackMax{"ptTrackMax", 10.0f, "Upper edge of track pT axis [GeV/c]"};
  Configurable<float> pxpyTrackMax{"pxpyTrackMax", 10.0f, "Upper edge of track px/py axis [GeV/c]"};

  // --- 3D histogram binning (Mass x Rapidity x pT, filled before pair kinematic cuts) ---
  Configurable<int> nBinsMass3D{"nBinsMass3D", 200, "Mass bins for 3D histogram"};
  Configurable<int> nBinsRap3D{"nBinsRap3D", 50, "Rapidity bins for 3D histogram"};
  Configurable<float> rapMin3D{"rapMin3D", -5.0f, "Min rapidity for 3D histogram"};
  Configurable<float> rapMax3D{"rapMax3D", 0.0f, "Max rapidity for 3D histogram"};
  Configurable<int> nBinsPt3D{"nBinsPt3D", 100, "pT bins for 3D histogram"};
  Configurable<float> ptMax3D{"ptMax3D", 1.0f, "Upper edge of pair pT axis for 3D histogram [GeV/c]"};

  // --- High-pT tail diagnostics (coherent pT mismeasurement tail / Delta-phi edge / MFT-MCH matching bias) ---
  // All filled before the pair pT/rapidity/mass cuts so the full high-pT tail is visible.
  Configurable<int> nBinsDeltaPt{"nBinsDeltaPt", 200, "Bins on pair #Delta#it{p}_{T} (reco - MC) axis"};
  Configurable<float> deltaPtMax{"deltaPtMax", 2.0f, "Max |#Delta#it{p}_{T}| for the coherent high-#it{p}_{T} tail check [GeV/c]"};
  Configurable<int> nBinsPtTail{"nBinsPtTail", 200, "Bins on pair #it{p}_{T} axis used for the Delta-phi tail check"};
  Configurable<float> ptTailMax{"ptTailMax", 2.0f, "Upper edge of pair #it{p}_{T} axis used for the Delta-phi tail check [GeV/c]"};
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 180, "Bins on |#Delta#phi| axis between the two reco muons"};

  static constexpr int kMuonPDG = 13;
  static constexpr int kJpsiPDG = 443;
  float mMu = 0.0f;

  // ===========================================================================
  // Initialization
  // ===========================================================================
  void init(InitContext&)
  {
    auto* particle = TDatabasePDG::Instance()->GetParticle(kMuonPDG);
    mMu = particle ? particle->Mass() : 0.105658f;

    initControlHistograms();
    initEfficiencySingleTrackHistograms();
    initEfficiencyPairHistograms();
    initPairMassRapPt3DHistograms();
    initResolutionSingleTrackHistograms();
    initResolutionPairHistograms();
  }

  // ---------------------------------------------------------------------------
  void initControlHistograms()
  {
    // --- Cut-flow for MC-truth level ---
    auto hCutFlowMC = registry.add<TH1>("hCutFlowMC", "MC Truth Cut Flow;;Counts",
                                        HistType::kTH1I, {{18, 0., 18.}});
    // --- Cut-flow for reco level ---
    auto hCutFlowReco = registry.add<TH1>("hCutFlowReco", "Reco Cut Flow;;Counts",
                                          HistType::kTH1I, {{18, 0., 18.}});

    TString CutNames[18] = {
      "0: All Cand",
      "1: Has Requested Track Type",
      "2: Pass Exact 2 Tracks",
      "3: Pass Exact MatchMFT",
      "4: Track All",
      "5: Track Pass Type",
      "6: Track Pass rAbs",
      "7: Track Pass pDCA",
      "8: Track Pass MatchMFT",
      "9: Track Pass Eta",
      "10: Track Pass pT",
      "11: Pair All",
      "12: Pair Unlike-sign",
      "13: Pair Both MC Muon",
      "14: Pair Pass Pt",
      "15: Pair Pass Rapidity",
      "16: Pair Pass Mass",
      "17: Filler"};
    for (int i = 0; i < 18; i++) {
      hCutFlowMC->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
      hCutFlowReco->GetXaxis()->SetBinLabel(i + 1, CutNames[i].Data());
    }

    // Track-type distribution (reco only)
    auto hTrackType = registry.add<TH1>("hTrackType", "Track Type Distribution;Track Type;Counts",
                                        HistType::kTH1I, {{5, -0.5, 4.5}});
    hTrackType->GetXaxis()->SetBinLabel(1, "GlobalMuon");
    hTrackType->GetXaxis()->SetBinLabel(2, "OtherMatch");
    hTrackType->GetXaxis()->SetBinLabel(3, "GlobalFwd");
    hTrackType->GetXaxis()->SetBinLabel(4, "MuonStandalone");
    hTrackType->GetXaxis()->SetBinLabel(5, "MCHStandalone");

    registry.add<TH1>("hNTracksTotal", "Total Forward Tracks per Event;N Tracks;Events",
                      HistType::kTH1I, {{10, -0.5, 10.5}});
    registry.add<TH1>("hNGlobalMuons", "Global Muon Tracks per Event;N Global Muons;Events",
                      HistType::kTH1I, {{20, -0.5, 19.5}});
    registry.add<TH1>("hMuonMultMC", "MC muons in acceptance per event;N muons;Events",
                      HistType::kTH1I, {{10, 0., 10.}});
    registry.add<TH1>("hMuonMultReco", "Reco muons in acceptance per event;N muons;Events",
                      HistType::kTH1I, {{10, 0., 10.}});

    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("eventCounter", "Processed Events", kTH1F, {axisCounter});
  }

  // ---------------------------------------------------------------------------
  // Efficiency histograms — single track (Dont use Acceptance efficiency denominator, this histograms are required to muon goes detector, not Jpsi)
  // ---------------------------------------------------------------------------
  void initEfficiencySingleTrackHistograms()
  {
    const AxisSpec axPdgCode{10000, -5000.f, 5000.f, "PDG Code"};
    registry.add("hPdgCodeMCAll", "PDG Codes of all MC particles", kTH1I, {axPdgCode});
    registry.add("hPdgCodeReco", "PDG Codes of Reco matched particles", kTH1I, {axPdgCode});

    const AxisSpec axPhiMCAll{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC All} (rad)"};
    const AxisSpec axPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC} (rad)"};
    const AxisSpec axPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{reco} (rad)"};

    const AxisSpec axEtaMCAll{nBinsEta, -10.f, 10.f, "#eta^{MC All}"};
    const AxisSpec axEtaMC{nBinsEta, etaMin, etaMax, "#eta^{MC}"};
    const AxisSpec axEtaReco{nBinsEta, etaMin, etaMax, "#eta^{reco}"};

    const AxisSpec axPtMCAll{nBinsP, 0.f, 10.f, "#it{p}_{T}^{MC All} (GeV/#it{c})"};
    const AxisSpec axPtMC{nBinsP, pTmin, ptTrackMax, "#it{p}_{T}^{MC} (GeV/#it{c})"};
    const AxisSpec axPtReco{nBinsP, pTmin, ptTrackMax, "#it{p}_{T}^{reco} (GeV/#it{c})"};

    registry.add("hEffTrkPhiMCAll", "Single track #phi MC All (eff)", kTH1F, {axPhiMCAll});
    registry.add("hEffTrkPhiMC", "Single track #phi MC (eff)", kTH1F, {axPhiMC});
    registry.add("hEffTrkPhiReco", "Single track #phi reco (eff)", kTH1F, {axPhiReco});

    registry.add("hEffTrkEtaMCAll", "Single track #eta MC All (eff)", kTH1F, {axEtaMCAll});
    registry.add("hEffTrkEtaMC", "Single track #eta MC (eff)", kTH1F, {axEtaMC});
    registry.add("hEffTrkEtaReco", "Single track #eta reco (eff)", kTH1F, {axEtaReco});

    registry.add("hEffTrkPtMCAll", "Single track #it{p}_{T} MC All (eff)", kTH1F, {axPtMCAll});
    registry.add("hEffTrkPtMC", "Single track #it{p}_{T} MC (eff)", kTH1F, {axPtMC});
    registry.add("hEffTrkPtReco", "Single track #it{p}_{T} reco (eff)", kTH1F, {axPtReco});

    // (eta,phi) maps, used to build 2D (Acc x Eff) maps and compare
    // Global-vs-Standalone track-type reconstruction locally in eta-phi space.
    // hEffTrkEtaPhiMC uses the same acceptance-cut MC-truth population as
    // hEffTrkEtaMC/hEffTrkPhiMC; hEffTrkEtaPhiReco uses the same reco-track
    // population as hEffTrkEtaReco/hEffTrkPhiReco (single reqTrackType per run).
    registry.add("hEffTrkEtaPhiMC", "Single track #eta vs #phi MC (eff)", kTH2F, {axEtaMC, axPhiMC});
    registry.add("hEffTrkEtaPhiReco", "Single track #eta vs #phi reco (eff)", kTH2F, {axEtaReco, axPhiReco});
  }

  // ---------------------------------------------------------------------------
  // Efficiency histograms — dimuon pair 
  // ---------------------------------------------------------------------------
  void initEfficiencyPairHistograms()
  {
    const AxisSpec axPairPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{pair,MC} (rad)"};
    const AxisSpec axPairPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{pair,reco} (rad)"};

    const AxisSpec axPairRapMC{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}^{pair,MC}"};
    const AxisSpec axPairRapReco{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}^{pair,reco}"};

    const AxisSpec axPairPtMC{nBinsPt, 0.f, ptMax, "#it{p}_{T}^{pair,MC} (GeV/#it{c})"};
    const AxisSpec axPairPtReco{nBinsPt, 0.f, ptMax, "#it{p}_{T}^{pair,reco} (GeV/#it{c})"};

    const AxisSpec axPairPt2MC{nBinsPt2, 0.f, pt2Max, "#it{p}_{T}^{2,pair,MC} (GeV^{2}/#it{c}^{2})"};
    const AxisSpec axPairPt2Reco{nBinsPt2, 0.f, pt2Max, "#it{p}_{T}^{2,pair,reco} (GeV^{2}/#it{c}^{2})"};

    const AxisSpec axPairMassMC{nBinsMass, massAxisMin, massAxisMax, "#it{M}^{pair,MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassReco{nBinsMass, massAxisMin, massAxisMax, "#it{M}^{pair,reco} (GeV/#it{c}^{2})"};

    // Dont use Acceptance efficiency denominator, this histograms are required to both muons goes detector, not Jpsi
    registry.add("hEffPairPhiMC", "Pair #phi MC (eff)", kTH1F, {axPairPhiMC});
    registry.add("hEffPairPhiReco", "Pair #phi reco (eff)", kTH1F, {axPairPhiReco});

    registry.add("hEffPairRapidityMC", "Pair rapidity MC (eff)", kTH1F, {axPairRapMC});
    registry.add("hEffPairRapidityReco", "Pair rapidity reco (eff)", kTH1F, {axPairRapReco});

    registry.add("hEffPairPtMC", "Pair #it{p}_{T} MC (eff)", kTH1F, {axPairPtMC});
    registry.add("hEffPairPtReco", "Pair #it{p}_{T} reco (eff)", kTH1F, {axPairPtReco});

    registry.add("hEffPairPt2MC", "Pair #it{p}_{T}^{2} MC (eff)", kTH1F, {axPairPt2MC});
    registry.add("hEffPairPt2Reco", "Pair #it{p}_{T}^{2} reco (eff)", kTH1F, {axPairPt2Reco});

    registry.add("hEffPairMassMC", "Pair mass MC (eff)", kTH1F, {axPairMassMC});
    registry.add("hEffPairMassReco", "Pair mass reco (eff)", kTH1F, {axPairMassReco});

    // --- Acc x Eff denominator: J/psi parent particle ---
    // Rapidity axis covers the wider generator window (jpsiRapMin to jpsiRapMax)
    const AxisSpec axJpsiRap{nBinsEta, jpsiRapMin, jpsiRapMax,
                             "#it{y}_{J/#psi}^{MC}"};
    // pT axis is intentionally identical to hEffPairPtReco for bin-by-bin division
    const AxisSpec axJpsiPt{nBinsPt, 0.f, ptMax,
                            "#it{p}_{T,J/#psi}^{MC} (GeV/#it{c})"};
    const AxisSpec axJpsiMass{200, 2.5f, 3.7f,
                              "#it{M}_{J/#psi}^{MC} (GeV/#it{c}^{2})"};

    registry.add("hAccEffDenomJpsiRap",  "J/#psi rapidity (Acc#timesEff denom)",
                 kTH1F, {axJpsiRap});
    registry.add("hAccEffDenomJpsiPt",   "J/#psi #it{p}_{T} (Acc#timesEff denom)",
                 kTH1F, {axJpsiPt});
    registry.add("hAccEffDenomJpsiMass", "J/#psi mass sanity check",
                 kTH1F, {axJpsiMass});
  }

  // ---------------------------------------------------------------------------
  // 3D histogram: Mass x Rapidity x pT (filled before pair kinematic cuts)
  // ---------------------------------------------------------------------------
  void initPairMassRapPt3DHistograms()
  {
    const AxisSpec axisMass3D{nBinsMass3D, pairMassMin, pairMassMax, "#it{M}_{#mu#mu} (GeV/#it{c}^{2})"};
    const AxisSpec axisRap3D{nBinsRap3D, rapMin3D, rapMax3D, "#it{y}_{#mu#mu}"};
    const AxisSpec axisPt3D{nBinsPt3D, 0.f, ptMax3D, "#it{p}_{T,#mu#mu} (GeV/#it{c})"};

    registry.add("hMassRapPt3DMC", "Mass vs rapidity vs pair pT, MC (before pair cuts);#it{M}_{#mu#mu} (GeV/#it{c}^{2});#it{y}_{#mu#mu};#it{p}_{T,#mu#mu} (GeV/#it{c})", kTH3D, {axisMass3D, axisRap3D, axisPt3D});
    registry.add("hMassRapPt3DReco", "Mass vs rapidity vs pair pT, reco (before pair cuts);#it{M}_{#mu#mu} (GeV/#it{c}^{2});#it{y}_{#mu#mu};#it{p}_{T,#mu#mu} (GeV/#it{c})", kTH3D, {axisMass3D, axisRap3D, axisPt3D});

    // --- Mass vs eta vs p (acceptance-edge check), filled before any single-muon cuts ---
    const AxisSpec axisEtaMuonQA{nBinsEta, etaMin, etaMax, "#eta_{#mu}"};
    const AxisSpec axisPMuonQA{nBinsP, 0.f, pTrackMax, "#it{p}_{#mu} (GeV/#it{c})"};
    registry.add("hMassEtaPMuon", "Pair mass vs single muon #eta vs #it{p} (both legs, reco);#it{M}_{#mu#mu} (GeV/#it{c}^{2});#eta_{#mu};#it{p}_{#mu} (GeV/#it{c})", kTH3D, {axisMass3D, axisEtaMuonQA, axisPMuonQA});
    // Per pair, only the leg closer to either acceptance edge (etaMin/etaMax) is filled here,
    // to check whether one decay lepton is piling up against the detector acceptance boundary.
    registry.add("hMassEtaPMuonEdge", "Pair mass vs single muon #eta vs #it{p}, leg closest to acceptance edge (reco);#it{M}_{#mu#mu} (GeV/#it{c}^{2});#eta_{#mu};#it{p}_{#mu} (GeV/#it{c})", kTH3D, {axisMass3D, axisEtaMuonQA, axisPMuonQA});

    // --- Coherent high-pT tail check: pair Delta-pT (reco - MC) vs reco mass, before pair cuts.
    // Run separately on Coherent and Incoherent MC samples; a diagonal mass-vs-pT tail here
    // signals a mismeasured leg dragging the pair pT up while the mass drops.
    const AxisSpec axisDeltaPt{nBinsDeltaPt, -deltaPtMax, deltaPtMax, "#it{p}_{T,#mu#mu}^{reco} - #it{p}_{T,#mu#mu}^{MC} (GeV/#it{c})"};
    registry.add("hDeltaPtVsMassReco", "Pair #Delta#it{p}_{T} (reco-MC) vs reco mass (before pair cuts);#it{p}_{T,#mu#mu}^{reco} - #it{p}_{T,#mu#mu}^{MC} (GeV/#it{c});#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})", kTH2F, {axisDeltaPt, axisMass3D});

    // --- Delta-phi between the two reco muons vs pair reco pT, before pair cuts.
    // Compare the high-pT slice of this histogram between Coherent and Incoherent samples:
    // an unnaturally sharp edge near pi in the Coherent high-pT tail indicates pair-pT smearing
    // artificially squeezing the two legs together/apart in phi.
    const AxisSpec axisPtTail{nBinsPtTail, 0.f, ptTailMax, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axisDeltaPhi{nBinsDeltaPhi, 0.f, TMath::Pi(), "|#Delta#phi^{reco}_{#mu#mu}| (rad)"};
    registry.add("hDeltaPhiVsPairPtReco", "|#Delta#phi| between reco muons vs pair #it{p}_{T} (before pair cuts);#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c});|#Delta#phi^{reco}_{#mu#mu}| (rad)", kTH2F, {axisPtTail, axisDeltaPhi});
  }

  // ---------------------------------------------------------------------------
  // Resolution histograms — single track
  // ---------------------------------------------------------------------------
  void initResolutionSingleTrackHistograms()
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

    // Note: reco-side phi/eta/pT marginal distributions are not duplicated here —
    // they are identical in content to hEffTrkPhiReco/hEffTrkEtaReco/hEffTrkPtReco
    // (same selection, same fill value), so only the MC-side truth distributions
    // of matched-and-reconstructed tracks are kept.
    registry.add("hResoTrkPhiMC", "Track #phi MC (reso)", kTH1F, {axPhiMC});
    registry.add("hResoTrkEtaMC", "Track #eta MC (reso)", kTH1F, {axEtaMC});

    registry.add("hResoTrkPMC", "Track #it{p} MC", kTH1F, {axPMC});
    registry.add("hResoTrkPReco", "Track #it{p} Reco", kTH1F, {axPReco});
    registry.add("hResoTrkPtMC", "Track #it{p}_{T} MC (reso)", kTH1F, {axPtMC});
    registry.add("hResoTrkPxMC", "Track #it{p}_{x} MC", kTH1F, {axPxMC});
    registry.add("hResoTrkPxReco", "Track #it{p}_{x} Reco", kTH1F, {axPxReco});
    registry.add("hResoTrkPyMC", "Track #it{p}_{y} MC", kTH1F, {axPyMC});
    registry.add("hResoTrkPyReco", "Track #it{p}_{y} Reco", kTH1F, {axPyReco});
    registry.add("hResoTrkPzMC", "Track #it{p}_{z} MC", kTH1F, {axPzMC});
    registry.add("hResoTrkPzReco", "Track #it{p}_{z} Reco", kTH1F, {axPzReco});

    registry.add("hResoPhi", "Track #phi Resolution", kTH2F, {axPhiMC, axPhiRes});
    registry.add("hResoEta", "Track #eta Resolution", kTH2F, {axEtaMC, axEtaRes});
    registry.add("hResoP", "Track #it{p} Resolution", kTH2F, {axPMC, axPRelRes});
    registry.add("hResoPt", "Track #it{p}_{T} Resolution", kTH2F, {axPtMC, axPRelRes});
    registry.add("hResoPx", "Track #it{p}_{x} Resolution", kTH2F, {axPxMC, axPRelRes});
    registry.add("hResoPy", "Track #it{p}_{y} Resolution", kTH2F, {axPyMC, axPRelRes});
    registry.add("hResoPz", "Track #it{p}_{z} Resolution", kTH2F, {axPzMC, axPRelRes});

    // --- MFT-MCH matching bias check: eta/phi residual vs MC p, low-p region is where
    // multiple-scattering-driven angle bias (track opening/closing) is expected to show up.
    registry.add("hResoEtaVsP", "Track #eta Resolution vs #it{p}^{MC} (MFT-MCH matching bias)", kTH2F, {axPMC, axEtaRes});
    registry.add("hResoPhiVsP", "Track #phi Resolution vs #it{p}^{MC} (MFT-MCH matching bias)", kTH2F, {axPMC, axPhiRes});
  }

  // ---------------------------------------------------------------------------
  // Resolution histograms — dimuon pair
  // ---------------------------------------------------------------------------
  void initResolutionPairHistograms()
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

    registry.add("hPairPtMC_PreCut", "Dimuon #it{p}_{T} MC (Pre Kin. Cuts)", kTH1F, {axisPairPtMC});
    registry.add("hPairPtReco_PreCut", "Dimuon #it{p}_{T} Reco (Pre Kin. Cuts)", kTH1F, {axisPairPtReco});
    registry.add("hPairPtMC_PostCut", "Dimuon #it{p}_{T} MC (Post All Cuts)", kTH1F, {axisPairPtMC});
    registry.add("hPairPtReco_PostCut", "Dimuon #it{p}_{T} Reco (Post All Cuts)", kTH1F, {axisPairPtReco});

    registry.add("hPairPt2MC_PreCut", "Dimuon #it{p}_{T}^{2} MC (Pre Kin. Cuts)", kTH1F, {axisPairPt2MC});
    registry.add("hPairPt2Reco_PreCut", "Dimuon #it{p}_{T}^{2} Reco (Pre Kin. Cuts)", kTH1F, {axisPairPt2Reco});
    registry.add("hPairPt2MC_PostCut", "Dimuon #it{p}_{T}^{2} MC (Post All Cuts)", kTH1F, {axisPairPt2MC});
    registry.add("hPairPt2Reco_PostCut", "Dimuon #it{p}_{T}^{2} Reco (Post All Cuts)", kTH1F, {axisPairPt2Reco});

    registry.add("hPairRapMC_PostCut", "Dimuon rapidity MC (Post All Cuts)", kTH1F, {axPairRapMC});
    registry.add("hPairRapReco_PostCut", "Dimuon rapidity Reco (Post All Cuts)", kTH1F, {axPairRapReco});
    registry.add("hPairEtaMC_PostCut", "Dimuon #eta MC (Post All Cuts)", kTH1F, {axPairEtaMC});
    registry.add("hPairEtaReco_PostCut", "Dimuon #eta Reco (Post All Cuts)", kTH1F, {axPairEtaReco});
    registry.add("hPairPhiMC_PostCut", "Dimuon #phi MC (Post All Cuts)", kTH1F, {axPairPhiMC});
    registry.add("hPairPhiReco_PostCut", "Dimuon #phi Reco (Post All Cuts)", kTH1F, {axPairPhiReco});
    registry.add("hPairMassMC_PostCut", "Dimuon mass MC (Post All Cuts)", kTH1F, {axPairMassMC});
    registry.add("hPairMassReco_PostCut", "Dimuon mass Reco (Post All Cuts)", kTH1F, {axPairMassReco});

    registry.add("hResponseMatrixPairPt", "Pair #it{p}_{T} Response Matrix", kTH2F, {axisPairPtMC, axisPairPtReco});
    registry.add("hResponseMatrixPairPt2", "Pair #it{p}_{T}^{2} Response Matrix", kTH2F, {axisPairPt2MC, axisPairPt2Reco});
    registry.add("hResponseMatrixPairRap", "Pair Rapidity Response Matrix", kTH2F, {axPairRapMC, axPairRapReco});
    registry.add("hResponseMatrixPairEta", "Pair #eta Response Matrix", kTH2F, {axPairEtaMC, axPairEtaReco});
    registry.add("hResponseMatrixPairPhi", "Pair #phi Response Matrix", kTH2F, {axPairPhiMC, axPairPhiReco});
    registry.add("hResponseMatrixPairMass", "Pair Mass Response Matrix", kTH2F, {axPairMassMC, axPairMassReco});

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
  // Core methods
  // ===========================================================================

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

    if (tr.chi2() > maxChi2)
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
    registry.fill(HIST("hCutFlowReco"), 9); // 9: Track Pass Eta

    if (recoVec.Pt() < pTmin)
      return false;
    registry.fill(HIST("hCutFlowReco"), 10); // 10: Track Pass pT

    return true;
  }

  template <typename TTrack, typename TMcParticle>
  void fillSingleTrackResolution(const TTrack& tr, const TMcParticle& mc)
  {
    TLorentzVector vReco, vMC;
    vReco.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    vMC.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);

    registry.fill(HIST("hResoTrkPhiMC"), vMC.Phi());
    registry.fill(HIST("hResoTrkEtaMC"), vMC.Eta());

    registry.fill(HIST("hResoTrkPMC"), vMC.P());
    registry.fill(HIST("hResoTrkPReco"), vReco.P());
    registry.fill(HIST("hResoTrkPtMC"), vMC.Pt());
    registry.fill(HIST("hResoTrkPxMC"), vMC.Px());
    registry.fill(HIST("hResoTrkPxReco"), vReco.Px());
    registry.fill(HIST("hResoTrkPyMC"), vMC.Py());
    registry.fill(HIST("hResoTrkPyReco"), vReco.Py());
    registry.fill(HIST("hResoTrkPzMC"), vMC.Pz());
    registry.fill(HIST("hResoTrkPzReco"), vReco.Pz());

    registry.fill(HIST("hResoPhi"), vMC.Phi(), TVector2::Phi_mpi_pi(vReco.Phi() - vMC.Phi()));
    registry.fill(HIST("hResoEta"), vMC.Eta(), vReco.Eta() - vMC.Eta());

    // MFT-MCH matching bias: eta/phi residual vs MC p (low-p emphasises multiple-scattering bias)
    registry.fill(HIST("hResoEtaVsP"), vMC.P(), vReco.Eta() - vMC.Eta());
    registry.fill(HIST("hResoPhiVsP"), vMC.P(), TVector2::Phi_mpi_pi(vReco.Phi() - vMC.Phi()));

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

  template <typename TTrack, typename TMcParticle>
  void fillPairAnalysis(const TTrack& tr1, const TMcParticle& mc1,
                        const TTrack& tr2, const TMcParticle& mc2)
  {
    registry.fill(HIST("hCutFlowReco"), 11); // 11: Pair All

    if (tr1.sign() * tr2.sign() >= 0)
      return;
    registry.fill(HIST("hCutFlowReco"), 12); // 12: Pair Unlike-sign

    if (std::abs(mc1.pdgCode()) != kMuonPDG || std::abs(mc2.pdgCode()) != kMuonPDG)
      return;
    registry.fill(HIST("hCutFlowReco"), 13); // 13: Pair Both MC Muon

    TLorentzVector recoVec1, recoVec2, mcVec1, mcVec2;
    recoVec1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mMu);
    recoVec2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mMu);
    mcVec1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);
    mcVec2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);

    TLorentzVector pairReco = recoVec1 + recoVec2;
    TLorentzVector pairMC = mcVec1 + mcVec2;

    registry.fill(HIST("hMassEtaPMuon"), pairReco.M(), recoVec1.Eta(), recoVec1.P());
    registry.fill(HIST("hMassEtaPMuon"), pairReco.M(), recoVec2.Eta(), recoVec2.P());

    // Fill only the leg closer to either acceptance edge (etaMin/etaMax),
    // to check whether one decay lepton is piling up against the boundary.
    float distToEdge1 = std::min(recoVec1.Eta() - etaMin, etaMax - recoVec1.Eta());
    float distToEdge2 = std::min(recoVec2.Eta() - etaMin, etaMax - recoVec2.Eta());
    if (distToEdge1 <= distToEdge2) {
      registry.fill(HIST("hMassEtaPMuonEdge"), pairReco.M(), recoVec1.Eta(), recoVec1.P());
    } else {
      registry.fill(HIST("hMassEtaPMuonEdge"), pairReco.M(), recoVec2.Eta(), recoVec2.P());
    }

    const float pairPtMC = pairMC.Pt();
    const float pairPtReco = pairReco.Pt();
    const float pairPt2MC = pairPtMC * pairPtMC;
    const float pairPt2Reco = pairPtReco * pairPtReco;
    const float pairRapMC = pairMC.Rapidity();
    const float pairRapReco = pairReco.Rapidity();
    const float pairEtaMC = pairMC.Eta();
    const float pairEtaReco = pairReco.Eta();
    const float pairPhiMC = pairMC.Phi();
    const float pairPhiReco = pairReco.Phi();
    const float pairMassMC = pairMC.M();
    const float pairMassReco = pairReco.M();

    registry.fill(HIST("hPairPtMC_PreCut"), pairPtMC);
    registry.fill(HIST("hPairPtReco_PreCut"), pairPtReco);
    registry.fill(HIST("hPairPt2MC_PreCut"), pairPt2MC);
    registry.fill(HIST("hPairPt2Reco_PreCut"), pairPt2Reco);

    registry.fill(HIST("hMassRapPt3DReco"), pairMassReco, pairRapReco, pairPtReco);

    // Coherent high-pT tail: Delta-pT (reco-MC) vs reco mass, before pair cuts.
    registry.fill(HIST("hDeltaPtVsMassReco"), pairPtReco - pairPtMC, pairMassReco);

    // Delta-phi between the two reco muons vs pair reco pT, before pair cuts.
    float deltaPhiReco = std::abs(TVector2::Phi_mpi_pi(recoVec1.Phi() - recoVec2.Phi()));
    registry.fill(HIST("hDeltaPhiVsPairPtReco"), pairPtReco, deltaPhiReco);

    if (pairReco.Pt() >= pairPtMax)
      return;
    registry.fill(HIST("hCutFlowReco"), 14); // 14: Pair Pass Pt

    if (pairReco.Rapidity() <= pairRapidityMin || pairReco.Rapidity() >= pairRapidityMax)
      return;
    registry.fill(HIST("hCutFlowReco"), 15); // 15: Pair Pass Rapidity

    if (pairReco.M() <= pairMassMin || pairReco.M() >= pairMassMax)
      return;
    registry.fill(HIST("hCutFlowReco"), 16); // 16: Pair Pass Mass

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
    if (pairPtMC > 0.f)
      registry.fill(HIST("hPairPtResoVsPtMC"), pairPtMC, (pairPtReco - pairPtMC) / pairPtMC);

    registry.fill(HIST("hResponseMatrixPairPt2"), pairPt2MC, pairPt2Reco);
    if (pairPt2MC > 0.f)
      registry.fill(HIST("hPairPt2ResoVsPt2MC"), pairPt2MC, (pairPt2Reco - pairPt2MC) / pairPt2MC);

    registry.fill(HIST("hEffPairPhiReco"), pairPhiReco);
    registry.fill(HIST("hEffPairRapidityReco"), pairRapReco);
    registry.fill(HIST("hEffPairPtReco"), pairPtReco);
    registry.fill(HIST("hEffPairPt2Reco"), pairPt2Reco);
    registry.fill(HIST("hEffPairMassReco"), pairMassReco);
  }

  // ===========================================================================
  // Process: MC truth
  // ===========================================================================
  void processMCTrue(o2::aod::UDMcCollisions const& mcCollisions,
                     o2::aod::UDMcParticles const& mcParticles)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> muonsPerMcColl;
    // All muons with no η cut — used to build the Acc×Eff denominator
    std::unordered_map<int32_t, std::vector<int32_t>> allMuonsPerMcColl;

    for (const auto& mc : mcParticles) {
      registry.fill(HIST("hCutFlowMC"), 0);
      registry.fill(HIST("hPdgCodeMCAll"), mc.pdgCode());

      if (std::abs(mc.pdgCode()) != kMuonPDG)
        continue;
      registry.fill(HIST("hCutFlowMC"), 1);

      TLorentzVector v;
      v.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);
      registry.fill(HIST("hEffTrkPhiMCAll"), v.Phi());
      registry.fill(HIST("hEffTrkEtaMCAll"), v.Eta());
      registry.fill(HIST("hEffTrkPtMCAll"), v.Pt());

      // Collect all muons before acceptance cuts for Acc×Eff denominator
      allMuonsPerMcColl[mc.udMcCollisionId()].push_back(mc.globalIndex());

      if (v.Eta() <= etaMin || v.Eta() >= etaMax)
        continue;
      registry.fill(HIST("hCutFlowMC"), 9);

      if (v.Pt() < pTmin)
        continue;
      registry.fill(HIST("hCutFlowMC"), 10);

      registry.fill(HIST("hEffTrkPhiMC"), v.Phi());
      registry.fill(HIST("hEffTrkEtaMC"), v.Eta());
      registry.fill(HIST("hEffTrkEtaPhiMC"), v.Eta(), v.Phi());
      registry.fill(HIST("hEffTrkPtMC"), v.Pt());

      muonsPerMcColl[mc.udMcCollisionId()].push_back(mc.globalIndex());
    }

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

          registry.fill(HIST("hCutFlowMC"), 11);

          if (mc1.pdgCode() * mc2.pdgCode() >= 0)
            continue;
          registry.fill(HIST("hCutFlowMC"), 12);
          registry.fill(HIST("hCutFlowMC"), 13);

          TLorentzVector pair = v1 + v2;

          registry.fill(HIST("hMassRapPt3DMC"), pair.M(), pair.Rapidity(), pair.Pt());

          if (pair.Pt() >= pairPtMax)
            continue;
          registry.fill(HIST("hCutFlowMC"), 14);

          if (pair.Rapidity() <= pairRapidityMin || pair.Rapidity() >= pairRapidityMax)
            continue;
          registry.fill(HIST("hCutFlowMC"), 15);

          if (pair.M() <= pairMassMin || pair.M() >= pairMassMax)
            continue;
          registry.fill(HIST("hCutFlowMC"), 16);

          registry.fill(HIST("hEffPairPhiMC"), pair.Phi());
          registry.fill(HIST("hEffPairRapidityMC"), pair.Rapidity());
          registry.fill(HIST("hEffPairPtMC"), pair.Pt());
          registry.fill(HIST("hEffPairPt2MC"), pair.Pt() * pair.Pt());
          registry.fill(HIST("hEffPairMassMC"), pair.M());
        }
      }
    }

    // --- Acc x Eff denominator: J/ψ reconstructed from MC muon daughters (no η cut) ---
    for (const auto& item : allMuonsPerMcColl) {
      const auto& ids = item.second;
      for (size_t i = 0; i < ids.size(); ++i) {
        auto mc1 = mcParticles.iteratorAt(ids[i]);
        TLorentzVector v1;
        v1.SetXYZM(mc1.px(), mc1.py(), mc1.pz(), mMu);
        for (size_t j = i + 1; j < ids.size(); ++j) {
          auto mc2 = mcParticles.iteratorAt(ids[j]);
          if (mc1.pdgCode() * mc2.pdgCode() >= 0)
            continue; // require opposite charge (μ⁺μ⁻)
          TLorentzVector v2;
          v2.SetXYZM(mc2.px(), mc2.py(), mc2.pz(), mMu);
          TLorentzVector vJpsi = v1 + v2;
          float y = vJpsi.Rapidity();
          if (y <= jpsiRapMin || y >= jpsiRapMax)
            continue;
          registry.fill(HIST("hAccEffDenomJpsiRap"),  y);
          registry.fill(HIST("hAccEffDenomJpsiPt"),   vJpsi.Pt());
          registry.fill(HIST("hAccEffDenomJpsiMass"), vJpsi.M());
        }
      }
    }

    (void)mcCollisions;
  }

  // ===========================================================================
  // Process: MC reco
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

      registry.fill(HIST("hNTracksTotal"), static_cast<float>(trkIds.size()));

      if (trkIds.empty())
        continue;

      registry.fill(HIST("hCutFlowReco"), 0); // 0: All Cand

      // --- Filter by requested track type ---
      std::vector<int32_t> candidateTrkIds;
      candidateTrkIds.reserve(trkIds.size());

      for (auto idx : trkIds) {
        auto tr = fwdTracks.iteratorAt(idx);
        registry.fill(HIST("hTrackType"), static_cast<float>(tr.trackType()));
        if (tr.trackType() == reqTrackType) {
          candidateTrkIds.push_back(idx);
        }
      }

      registry.fill(HIST("hNGlobalMuons"), static_cast<float>(candidateTrkIds.size()));

      if (candidateTrkIds.empty())
        continue;
      registry.fill(HIST("hCutFlowReco"), 1); // 1: Has Requested Track Type

      // -----------------------------------------------------------------------
      // Single-track loop: quality cuts + kinematic acceptance.
      // Runs over ALL candidateTrkIds regardless of how many there are,
      // matching the original UPCMuonEfficiencyTask behaviour where single-
      // track histograms are filled before the 2-track requirement is applied.
      // -----------------------------------------------------------------------
      std::vector<int32_t> goodTrkIds;
      goodTrkIds.reserve(candidateTrkIds.size());

      for (auto idx : candidateTrkIds) {
        auto tr = fwdTracks.iteratorAt(idx);

        if (!passRecoTrackCuts(tr))
          continue;

        const auto& mc = tr.udMcParticle();
        registry.fill(HIST("hPdgCodeReco"), mc.pdgCode());

        if (std::abs(mc.pdgCode()) != kMuonPDG)
          continue;

        goodTrkIds.push_back(idx);

        // Efficiency: reco single-track histograms (filled here, before the
        // 2-track requirement below, consistent with the original task)
        TLorentzVector vReco;
        vReco.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
        registry.fill(HIST("hEffTrkPhiReco"), vReco.Phi());
        registry.fill(HIST("hEffTrkEtaReco"), vReco.Eta());
        registry.fill(HIST("hEffTrkEtaPhiReco"), vReco.Eta(), vReco.Phi());
        registry.fill(HIST("hEffTrkPtReco"), vReco.Pt());

        // Resolution: single-track histograms
        fillSingleTrackResolution(tr, mc);
      }

      registry.fill(HIST("hMuonMultReco"), static_cast<float>(goodTrkIds.size()));

      // --- 2-track requirement for cutflow step 2 ---
      if (candidateTrkIds.size() == static_cast<size_t>(reqMatchMFT))
        registry.fill(HIST("hCutFlowReco"), 2); // 2: Pass Exact 2 Tracks (candidate level)

      // --- Pair loop: only when exactly reqMatchMFT good muon tracks found ---
      if (goodTrkIds.size() != static_cast<size_t>(reqMatchMFT))
        continue;
      registry.fill(HIST("hCutFlowReco"), 3); // 3: Pass Exact MatchMFT

      for (size_t i = 0; i < goodTrkIds.size(); ++i) {
        auto tr1 = fwdTracks.iteratorAt(goodTrkIds[i]);
        const auto& mc1 = tr1.udMcParticle();

        for (size_t j = i + 1; j < goodTrkIds.size(); ++j) {
          auto tr2 = fwdTracks.iteratorAt(goodTrkIds[j]);
          const auto& mc2 = tr2.udMcParticle();

          fillPairAnalysis(tr1, mc1, tr2, mc2);
        }
      }

      (void)candId;
      (void)eventCandidates;
    }
  }

  PROCESS_SWITCH(UPCMuonAnalysisMC, processMCTrue,
                 "Fill MC-truth single-track and dimuon efficiency histograms", true);
  PROCESS_SWITCH(UPCMuonAnalysisMC, processMcReco,
                 "Fill reco single-track and dimuon resolution + efficiency histograms", true);
};

// ============================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCMuonAnalysisMC>(cfgc, TaskName{"upc-muon-analysis"}),
  };
}
