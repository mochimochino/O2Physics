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

/// \file   UPCMuonAnalysisAlignmentDebug.cxx
/// \brief  UPCMuonAnalysisMC + alignment-diagnostic histograms.
///         Added debug tasks on top of the base analysis (all original
///         histograms are preserved):
///
///         [Debug-1] (q/pT) residual per charge
///           (q/pT)^{reco} - (q/pT)^{MC} split into mu+ / mu-.
///           q/pT is the fitting variable of the track fitter and has
///           Gaussian errors, making rotation/twist between MFT and MCH
///           visible as a charge-dependent mean shift.
///
///         [Debug-2] Delta(Delta-phi) bimodal
///           (phi_{mu+} - phi_{mu-})^{reco} - (phi_{mu+} - phi_{mu-})^{MC}
///           for each unlike-sign muon pair.  A global phi rotation +/-delta
///           appears as a clean bimodal at +/-2*delta, independent of pT.
///
///         [Debug-3] (q/pT) residual vs total momentum p^{MC}
///           Material budget mismatch causes a curved (Bethe-Bloch) trend
///           at low p; a pure B-field scale error produces a flat offset.
///
///         [Debug-4] (q/pT) residual vs eta^{MC}
///           Material budget mismatch grows with |eta| as the absorber path
///           length increases; a pure B-field scale error is nearly flat.
///
///         [Debug-8] Pair mass resolution vs opening angle / DeltaR
///           M_{reco}-M_{MC} vs the 3D opening angle theta and vs
///           DeltaR=sqrt(Delta-eta^2+Delta-phi^2) between the two muons
///           (both MC-truth and reco binned). Since M^2 ~ p_1 p_2 theta^2,
///           a mass deficit that grows with the reco angle points to a
///           geometric (angle) mismeasurement rather than a momentum-scale
///           bias.
///
///         [Debug-9] Pair pT / mass resolution vs pair rapidity y^{MC}
///           Tests whether resolution degrades toward the edges of the
///           rapidity acceptance.
///
///         [Debug-10] Pair pT / mass resolution vs helicity angle cos(theta_HE)
///           cos(theta_HE) is the angle of mu+ in the dimuon rest frame
///           w.r.t. the dimuon's lab-frame flight direction. cos(theta_HE)~0
///           means the two muons split the pair momentum symmetrically;
///           cos(theta_HE)~+/-1 means one muon is soft, near the single-muon
///           pT/acceptance cut. Comparing the two regimes separates a
///           generic high-p resolution deficit from an acceptance-edge /
///           low-pT-cut distortion.
///
///         [Debug-11] Resolution vs number of MFT clusters
///           Single-track (q/pT) residual and pT resolution, and pair
///           mass/pT resolution (binned by min N_{MFT clusters} of the two
///           legs), vs the number of MFT clusters used in the track fit.
///           Requires the UD table producer to have been run with
///           fSaveMFTClusters=true so that UDFwdTracksCls is populated.
///
///         [Debug-12] Pair pT / mass vs number of MCH clusters
///           UDFwdTracksExtra_001::fNClusters (tr.nClusters(), no extra table
///           needed) binned by min(N_{MCH clusters}) of the two legs, vs
///           both the raw reco pair pT/mass distribution and the pT/mass
///           resolution (reco-MC), to separate a shift in the reco
///           distribution from a widening/shift of the resolution.
///
///         [Debug-13] Softer-leg total momentum p vs pair mass
///           Complements [Debug-8]: M^2 ~ p_1 p_2 theta^2, so this tests the
///           momentum half of the relation (the softer leg by total p^{MC},
///           not pT). Raw reco pair mass vs p^{MC}_{soft leg} separates an
///           acceptance/low-p edge effect; pair mass resolution vs the soft
///           leg's own p resolution tests directly whether its momentum
///           underestimation drives the mass deficit.
///
///         [Debug-14] MC-truth direct comparison + mass-deficit-tail gate
///           (1) theta^{reco}-theta^{MC} vs pair pT^{reco} (pre-kin. cuts):
///               tests whether an underestimated opening angle correlates
///               with an inflated reco pair pT.
///           (2) pT_mu^{reco}-pT_mu^{MC} (either leg) vs pair pT^{reco}
///               (pre-kin. cuts): separates one leg's pT being
///               overestimated (scenario A) from both legs' pT staying
///               accurate while the pair pT still rises (scenario B,
///               pointing back at (1)).
///           (3) Delta-theta / Delta-pT_mu distributions gated on the reco
///               mass tail [massTailMin,massTailMax] vs the inclusive
///               (post-cut) distributions, to see whether the mass-deficit
///               tail is populated by a distinct mismeasured sub-population.
///
///         [Debug-15] Verification method B: mass resolution normalized by
///         true pair pT
///           The pT^2-resolution cross term 2*p_{T,true}*Sigma(eps_x) scales
///           with the event's true pair pT, so it dominates the raw
///           (M_reco-M_MC) spread for incoherent events with a wide
///           p_{T,true} spread. Dividing the mass residual by
///           p_{T,mu mu}^{true} removes this scale dependence; the raw
///           (hPairMassResoVsPtMC) and normalized (hPairMassResoNormVsPtMC)
///           versions are both filled so the "before/after" can be compared
///           directly.
///
///         [Debug-16] Delta-pT_mu vs single-muon reco pT (Inclusive / MassTail)
///           hDebugDeltaPtMu_Inclusive/_MassTail ([Debug-14](3)) are 1D and
///           only show the overall shape of the momentum-underestimated
///           sub-population in the mass-deficit tail. This adds the same
///           Inclusive/MassTail split but 2D against the muon's own reco pT,
///           to localise *where* in reco-pT space that sub-population sits
///           (e.g. concentrated near the pTmin acceptance edge vs spread
///           across the full range).
///
///         [Debug-17] Delta-pT_mu vs single-muon phi (reco), and split by charge
///           (1) vs phi^{reco}: same Inclusive/MassTail gate, looking for a
///               detector-sector dependence (MFT disk boundary, MCH
///               quadrant) in the momentum-underestimated sub-population.
///           (2) split by mu+/mu-: if that sub-population is concentrated in
///               one charge, that points back at a charge-dependent bias
///               (e.g. an MFT/MCH rotation/twist, as in [Debug-1]).
///
///         [Debug-18] Approximate single-muon (x,y) position at the absorber end
///           True fitted (x,y) isn't available in the current UD skimmed
///           tables (UDFwdTracksProp/UDFwdTracksCovProp are not produced by
///           upcCandProducerGlobalMuon.cxx). Approximated instead from
///           already-available quantities as x~R_abs*cos(phi^{reco}),
///           y~R_abs*sin(phi^{reco}). Same Inclusive/MassTail gate as
///           [Debug-14](3), plus a third map restricted to legs whose own
///           pT is anomalously underestimated (< deltaPtMuAnomalousMax), to
///           directly localise the momentum-underestimated sub-population in
///           detector space.
///
///         [Debug-19] Delta-pT_mu vs pDca, and event vertex position
///           (1) vs pDca^{reco}: correlates the momentum residual directly
///               with pDca (both outputs of the same track fit), same
///               Inclusive/MassTail gate.
///           (2) Event (collision) vertex position (posX/posY/posZ from
///               UDCollisions, already joined as CandidatesFwd but
///               previously unused): a poorly-reconstructed vertex feeds
///               directly into the MCH track fit's DCA calculation, so a
///               vertex-position anomaly could explain a pDca/pT anomaly.
///               Same Inclusive/MassTail/AnomalousPt split as [Debug-18].
///
///         [Debug-20] Muon pT asymmetry (reco) vs pair pT (reco)
///           Exclusive UPC production predicts p_{T,mu1}^{true} ~
///           -p_{T,mu2}^{true} (near back-to-back), so a pair with
///           strongly asymmetric leg pT's is expected to show a larger
///           vector-sum (pair) pT than a balanced pair. Tests that
///           relationship with reco quantities only (no MC truth), same
///           Inclusive/MassTail gate as [Debug-16]/[Debug-17].
///
///         [Debug-21] Muon pair track-time difference
///           trackTime()/trackTimeRes() are already on UDFwdTracks (no extra
///           table join needed) but were unused until now. If the two legs
///           of a "pair" are actually mismatched to different bunch
///           crossings, that should show up as an anomalously large
///           |trackTime_mu1-trackTime_mu2|, which would also explain a wrong
///           vertex/cluster assignment feeding into the momentum anomaly
///           investigated in [Debug-14]-[Debug-20]. Step (1) here: plain
///           Inclusive/MassTail comparison; step (2) (vs Delta-pT_mu / pair
///           pT) to follow once a trend is seen.
///
///         [Debug-22] Pair mass vs single-muon pT / eta (continuous, both legs)
///           Unlike the massTailMin/massTailMax binary gate used throughout
///           [Debug-14] etc., this correlates the reco pair mass directly
///           and continuously with each leg's own reco pT/eta, to see
///           whether a mass distortion (not just the low-mass tail
///           specifically) tracks a particular single-muon pT or eta region
///           (e.g. low pT, or the eta acceptance edge) across the whole
///           mass range.
///
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

struct UPCMuonAnalysisAlignmentDebug {

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

  // --- [Debug-5] Pair resolution vs muon pT asymmetry ---
  Configurable<int> nBinsDeltaPtMuons{"nBinsDeltaPtMuons", 200,
                                      "Bins on (p_{T,mu+} - p_{T,mu-}) axis for the pair-resolution-vs-muon-pT-asymmetry check"};
  Configurable<float> deltaPtMuonsMax{"deltaPtMuonsMax", 5.0f,
                                      "Max |p_{T,mu+} - p_{T,mu-}| for the pair-resolution-vs-muon-pT-asymmetry check [GeV/c]"};

  // --- [Debug-8] Opening angle / DeltaR axes ---
  Configurable<float> openingAngleMax{"openingAngleMax", 1.0f,
                                      "Upper edge of pair opening-angle #theta axis [rad]"};
  Configurable<float> deltaRMax{"deltaRMax", 1.0f,
                                "Upper edge of pair #DeltaR = #sqrt{#Delta#eta^{2}+#Delta#phi^{2}} axis"};

  // --- [Debug-14] Mass-deficit tail gate ---
  Configurable<float> massTailMin{"massTailMin", 2.5f,
                                  "Lower edge of the reco mass tail gate used to isolate mass-deficit events [GeV/c^2]"};
  Configurable<float> massTailMax{"massTailMax", 3.0f,
                                  "Upper edge of the reco mass tail gate used to isolate mass-deficit events [GeV/c^2]"};

  // --- [Debug-18] Anomalous single-muon pT threshold ---
  Configurable<float> deltaPtMuAnomalousMax{"deltaPtMuAnomalousMax", -0.05f,
                                            "p_{T,mu}^{reco}-p_{T,mu}^{MC} threshold below which a leg is flagged as anomalously underestimated [GeV/c]"};

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
    initEfficiencySingleTrackHistograms();
    initEfficiencyPairHistograms();
    initResolutionSingleTrackHistograms();
    initResolutionPairHistograms();
    initDebugAlignmentHistograms();
    initDebugKinematicHistograms();
    initDebugMFTClusterHistograms();
    initDebugNClustersHistograms();
    initDebugSoftLegMomentumHistograms();
    initDebugMcTruthDirectComparisonHistograms();
    initDebugNormalizedMassResoHistograms();
    initDebugDeltaPtMuVsRecoPtHistograms();
    initDebugDeltaPtMuVsPhiAndSignHistograms();
    initDebugXYPositionHistograms();
    initDebugDcaVertexHistograms();
    initDebugDeltaPtMuonsRecoHistograms();
    initDebugDeltaTrackTimeHistograms();
    initDebugMassVsSingleMuonKinematicsHistograms();
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
  // Efficiency histograms — single track
  // ---------------------------------------------------------------------------
  void initEfficiencySingleTrackHistograms()
  {
    const AxisSpec axPdgCode{10000, -5000.f, 5000.f, "PDG Code"};
    registry.add("hPdgCodeMCAll", "PDG Codes of all MC particles", kTH1I, {axPdgCode});
    registry.add("hPdgCodeReco", "PDG Codes of Reco matched particles", kTH1I, {axPdgCode});

    const AxisSpec axPhiMCAll{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC All} (rad)"};
    const AxisSpec axPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{MC} (rad)"};
    const AxisSpec axPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{reco} (rad)"};
    const AxisSpec axPhiEff{120, 0.f, 1.2f, "Efficiency"};

    const AxisSpec axEtaMCAll{nBinsEta, -10.f, 10.f, "#eta^{MC All}"};
    const AxisSpec axEtaMC{nBinsEta, etaMin, etaMax, "#eta^{MC}"};
    const AxisSpec axEtaReco{nBinsEta, etaMin, etaMax, "#eta^{reco}"};
    const AxisSpec axEtaEff{120, 0.f, 1.2f, "Efficiency"};

    const AxisSpec axPtMCAll{nBinsP, 0.f, 10.f, "#it{p}_{T}^{MC All} (GeV/#it{c})"};
    const AxisSpec axPtMC{nBinsP, pTmin, ptTrackMax, "#it{p}_{T}^{MC} (GeV/#it{c})"};
    const AxisSpec axPtReco{nBinsP, pTmin, ptTrackMax, "#it{p}_{T}^{reco} (GeV/#it{c})"};
    const AxisSpec axPtEff{120, 0.f, 1.2f, "Efficiency"};

    registry.add("hEffTrkPhiMCAll", "Single track #phi MC All (eff)", kTH1F, {axPhiMCAll});
    registry.add("hEffTrkPhiMC", "Single track #phi MC (eff)", kTH1F, {axPhiMC});
    registry.add("hEffTrkPhiReco", "Single track #phi reco (eff)", kTH1F, {axPhiReco});
    registry.add("hEffTrkPhiEff", "Single track #phi efficiency", kTH1F, {axPhiEff});

    registry.add("hEffTrkEtaMCAll", "Single track #eta MC All (eff)", kTH1F, {axEtaMCAll});
    registry.add("hEffTrkEtaMC", "Single track #eta MC (eff)", kTH1F, {axEtaMC});
    registry.add("hEffTrkEtaReco", "Single track #eta reco (eff)", kTH1F, {axEtaReco});
    registry.add("hEffTrkEtaEff", "Single track #eta efficiency", kTH1F, {axEtaEff});

    registry.add("hEffTrkPtMCAll", "Single track #it{p}_{T} MC All (eff)", kTH1F, {axPtMCAll});
    registry.add("hEffTrkPtMC", "Single track #it{p}_{T} MC (eff)", kTH1F, {axPtMC});
    registry.add("hEffTrkPtReco", "Single track #it{p}_{T} reco (eff)", kTH1F, {axPtReco});
    registry.add("hEffTrkPtEff", "Single track #it{p}_{T} efficiency", kTH1F, {axPtEff});
  }

  // ---------------------------------------------------------------------------
  // Efficiency histograms — dimuon pair
  // ---------------------------------------------------------------------------
  void initEfficiencyPairHistograms()
  {
    const AxisSpec axPairPhiMC{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{pair,MC} (rad)"};
    const AxisSpec axPairPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi^{pair,reco} (rad)"};
    const AxisSpec axPairPhiEff{120, 0.f, 1.2f, "Efficiency"};

    const AxisSpec axPairRapMC{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}^{pair,MC}"};
    const AxisSpec axPairRapReco{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}^{pair,reco}"};
    const AxisSpec axPairRapEff{120, 0.f, 1.2f, "Efficiency"};

    const AxisSpec axPairPtMC{nBinsPt, 0.f, ptMax, "#it{p}_{T}^{pair,MC} (GeV/#it{c})"};
    const AxisSpec axPairPtReco{nBinsPt, 0.f, ptMax, "#it{p}_{T}^{pair,reco} (GeV/#it{c})"};
    const AxisSpec axPairPtEff{120, 0.f, 1.2f, "Efficiency"};

    const AxisSpec axPairPt2MC{nBinsPt2, 0.f, pt2Max, "#it{p}_{T}^{2,pair,MC} (GeV^{2}/#it{c}^{2})"};
    const AxisSpec axPairPt2Reco{nBinsPt2, 0.f, pt2Max, "#it{p}_{T}^{2,pair,reco} (GeV^{2}/#it{c}^{2})"};

    const AxisSpec axPairMassMC{nBinsMass, massAxisMin, massAxisMax, "#it{M}^{pair,MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassReco{nBinsMass, massAxisMin, massAxisMax, "#it{M}^{pair,reco} (GeV/#it{c}^{2})"};

    registry.add("hEffPairPhiMC", "Pair #phi MC (eff)", kTH1F, {axPairPhiMC});
    registry.add("hEffPairPhiReco", "Pair #phi reco (eff)", kTH1F, {axPairPhiReco});
    registry.add("hEffPairPhiEff", "Pair #phi efficiency", kTH1F, {axPairPhiEff});

    registry.add("hEffPairRapidityMC", "Pair rapidity MC (eff)", kTH1F, {axPairRapMC});
    registry.add("hEffPairRapidityReco", "Pair rapidity reco (eff)", kTH1F, {axPairRapReco});
    registry.add("hEffPairRapidityEff", "Pair rapidity efficiency", kTH1F, {axPairRapEff});

    registry.add("hEffPairPtMC", "Pair #it{p}_{T} MC (eff)", kTH1F, {axPairPtMC});
    registry.add("hEffPairPtReco", "Pair #it{p}_{T} reco (eff)", kTH1F, {axPairPtReco});
    registry.add("hEffPairPtEff", "Pair #it{p}_{T} efficiency", kTH1F, {axPairPtEff});

    registry.add("hEffPairPt2MC", "Pair #it{p}_{T}^{2} MC (eff)", kTH1F, {axPairPt2MC});
    registry.add("hEffPairPt2Reco", "Pair #it{p}_{T}^{2} reco (eff)", kTH1F, {axPairPt2Reco});

    registry.add("hEffPairMassMC", "Pair mass MC (eff)", kTH1F, {axPairMassMC});
    registry.add("hEffPairMassReco", "Pair mass reco (eff)", kTH1F, {axPairMassReco});
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

    registry.add("hResoTrkPhiMC", "Track #phi MC (reso)", kTH1F, {axPhiMC});
    registry.add("hResoTrkPhiReco", "Track #phi Reco (reso)", kTH1F, {axPhiReco});
    registry.add("hResoTrkEtaMC", "Track #eta MC (reso)", kTH1F, {axEtaMC});
    registry.add("hResoTrkEtaReco", "Track #eta Reco (reso)", kTH1F, {axEtaReco});

    registry.add("hResoTrkPMC", "Track #it{p} MC", kTH1F, {axPMC});
    registry.add("hResoTrkPReco", "Track #it{p} Reco", kTH1F, {axPReco});
    registry.add("hResoTrkPtMC", "Track #it{p}_{T} MC (reso)", kTH1F, {axPtMC});
    registry.add("hResoTrkPtReco", "Track #it{p}_{T} Reco (reso)", kTH1F, {axPtReco});
    registry.add("hResoTrkPxMC", "Track #it{p}_{x} MC", kTH1F, {axPxMC});
    registry.add("hResoTrkPxReco", "Track #it{p}_{x} Reco", kTH1F, {axPxReco});
    registry.add("hResoTrkPyMC", "Track #it{p}_{y} MC", kTH1F, {axPyMC});
    registry.add("hResoTrkPyReco", "Track #it{p}_{y} Reco", kTH1F, {axPyReco});
    registry.add("hResoTrkPzMC", "Track #it{p}_{z} MC", kTH1F, {axPzMC});
    registry.add("hResoTrkPzReco", "Track #it{p}_{z} Reco", kTH1F, {axPzReco});

    registry.add("hResoPhi", "Track #phi Resolution", kTH2F, {axPhiMC, axPhiRes});
    registry.add("hResoPhi_MuPlus", "Track #phi Resolution #mu^{+}", kTH2F, {axPhiMC, axPhiRes});
    registry.add("hResoPhi_MuMinus", "Track #phi Resolution #mu^{-}", kTH2F, {axPhiMC, axPhiRes});
    registry.add("hResoEta", "Track #eta Resolution", kTH2F, {axEtaMC, axEtaRes});
    registry.add("hResoP", "Track #it{p} Resolution", kTH2F, {axPMC, axPRelRes});
    registry.add("hResoPt", "Track #it{p}_{T} Resolution", kTH2F, {axPtMC, axPRelRes});
    registry.add("hResoPx", "Track #it{p}_{x} Resolution", kTH2F, {axPxMC, axPRelRes});
    registry.add("hResoPy", "Track #it{p}_{y} Resolution", kTH2F, {axPyMC, axPRelRes});
    registry.add("hResoPz", "Track #it{p}_{z} Resolution", kTH2F, {axPzMC, axPRelRes});
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

  // ---------------------------------------------------------------------------
  // [Debug] Alignment diagnostic histograms
  // ---------------------------------------------------------------------------
  void initDebugAlignmentHistograms()
  {
    // --- [Debug-1] (q/pT) residual per charge ---
    // A rotation/twist between MFT and MCH adds a spurious curvature that
    // has opposite sign for mu+ and mu-.  Because q/pT is the linear
    // fitting variable, its residual is Gaussian and the mean shift is
    // directly readable as the alignment bias in (GeV/c)^{-1}.
    const AxisSpec axQOverPtRes{400, -0.5f, 0.5f,
                                "(q/p_{T})^{reco}-(q/p_{T})^{MC} (GeV/c)^{-1}"};
    const AxisSpec axPtMCDbg{nBinsP, 0.f, ptTrackMax,
                             "p_{T}^{MC} (GeV/c)"};
    // p axis: Bethe-Bloch energy loss is strongest at low total momentum,
    // so material-budget mismatches appear as a curved (non-linear) trend
    // in residual vs p — unlike a flat B-field scale error.
    const AxisSpec axPMCDbg{nBinsP, 0.f, pTrackMax,
                            "p^{MC} (GeV/c)"};
    // eta axis: absorber path length varies with angle, so material-budget
    // mismatches produce an eta-dependent residual shift.
    const AxisSpec axEtaMCDbg{nBinsEta, etaMin, etaMax,
                              "#eta^{MC}"};

    registry.add("hDebugQOverPtResidual_MuPlus",
                 "Debug (q/p_{T}) residual #mu^{+}",
                 kTH1F, {axQOverPtRes});
    registry.add("hDebugQOverPtResidual_MuMinus",
                 "Debug (q/p_{T}) residual #mu^{-}",
                 kTH1F, {axQOverPtRes});
    registry.add("hDebugQOverPtResoVsPtMC_MuPlus",
                 "Debug (q/p_{T}) residual vs p_{T}^{MC} #mu^{+}",
                 kTH2F, {axPtMCDbg, axQOverPtRes});
    registry.add("hDebugQOverPtResoVsPtMC_MuMinus",
                 "Debug (q/p_{T}) residual vs p_{T}^{MC} #mu^{-}",
                 kTH2F, {axPtMCDbg, axQOverPtRes});
    // [Debug-3] residual vs total momentum p^{MC}
    // Material budget: curved (Bethe-Bloch) dependence at low p.
    // B-field scale error: flat offset independent of p.
    registry.add("hDebugQOverPtResoVsPMC_MuPlus",
                 "Debug (q/p_{T}) residual vs p^{MC} #mu^{+}",
                 kTH2F, {axPMCDbg, axQOverPtRes});
    registry.add("hDebugQOverPtResoVsPMC_MuMinus",
                 "Debug (q/p_{T}) residual vs p^{MC} #mu^{-}",
                 kTH2F, {axPMCDbg, axQOverPtRes});
    // [Debug-4] residual vs eta^{MC}
    // Material budget: absorber thickness increases at large |eta| (more
    // oblique path), giving an eta-dependent shift.
    // B-field scale error: nearly flat in eta.
    registry.add("hDebugQOverPtResoVsEtaMC_MuPlus",
                 "Debug (q/p_{T}) residual vs #eta^{MC} #mu^{+}",
                 kTH2F, {axEtaMCDbg, axQOverPtRes});
    registry.add("hDebugQOverPtResoVsEtaMC_MuMinus",
                 "Debug (q/p_{T}) residual vs #eta^{MC} #mu^{-}",
                 kTH2F, {axEtaMCDbg, axQOverPtRes});

    // --- [Debug-2] Delta(Delta-phi) bimodal ---
    // A global phi rotation by +/-delta shifts phi(mu+) and phi(mu-) in
    // opposite directions, so the difference-of-differences
    //   (phi_{mu+} - phi_{mu-})^{reco} - (phi_{mu+} - phi_{mu-})^{MC}
    // equals +/-2*delta and is independent of pT magnitude.
    // A non-zero alignment offset therefore appears as a clean bimodal at
    // +/-2*delta, completely free from momentum-magnitude smearing.
    // Filled for all unlike-sign MC-muon pairs (before kinematic cuts) to
    // maximise statistics.
    const AxisSpec axDeltaDeltaPhi{400, -0.1f, 0.1f,
                                   "(#phi_{#mu^{+}}-#phi_{#mu^{-}})^{reco}"
                                   "-(#phi_{#mu^{+}}-#phi_{#mu^{-}})^{MC} (rad)"};
    registry.add("hDebugDeltaDeltaPhi",
                 "Debug #Delta(#Delta#phi) bimodal (all #mu^{+}#mu^{-} pairs, pre-kin. cuts)",
                 kTH1F, {axDeltaDeltaPhi});
    const AxisSpec axPairPtMCDbg{nBinsPt, 0.f, ptMax,
                                 "#it{p}_{T,#mu#mu}^{MC} (GeV/#it{c})"};
    registry.add("hDebugDeltaDeltaPhiVsPtMC",
                 "Debug #Delta(#Delta#phi) vs #it{p}_{T,#mu#mu}^{MC} (pre-kin. cuts)",
                 kTH2F, {axPairPtMCDbg, axDeltaDeltaPhi});

    const AxisSpec axSinglePtMCDbg{nBinsPt, 0.f, ptTrackMax, "#it{p}_{T,#mu}^{MC} (GeV/#it{c})"};
    registry.add("hDebugDeltaDeltaPhiVsPtMCPos",
                 "Debug #Delta(#Delta#phi) vs #it{p}_{T,#mu}^{MC} (pre-kin. cuts, positive only)",
                 kTH2F, {axSinglePtMCDbg, axDeltaDeltaPhi});
    registry.add("hDebugDeltaDeltaPhiVsPtMCNeg",
                 "Debug #Delta(#Delta#phi) vs #it{p}_{T,#mu}^{MC} (pre-kin. cuts, negative only)",
                 kTH2F, {axSinglePtMCDbg, axDeltaDeltaPhi});

    // --- [Debug-7] Plain reco #Delta#phi (mu+ - mu-) vs pair-level kinematics ---
    // Unlike hDebugDeltaDeltaPhi (which uses the reco-MC residual to isolate a
    // pure alignment rotation), this looks at the raw reco opening angle in
    // phi directly, to check whether it correlates with the pair pT
    // resolution, the absolute reco pair pT, or the reco pair mass.
    const AxisSpec axDeltaPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(),
                                  "(#phi_{#mu^{+}}-#phi_{#mu^{-}})^{reco} (rad)"};
    const AxisSpec axPairPtDeltaDbg{nBinsPt, -ptMax, ptMax,
                                    "#it{p}_{T,#mu#mu}^{reco}-#it{p}_{T,#mu#mu}^{MC} (GeV/#it{c})"};
    const AxisSpec axPairPtRecoDbg{nBinsPt, 0.f, ptMax,
                                   "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axPairMassRecoDbg{nBinsMass, massAxisMin, massAxisMax,
                                     "#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})"};

    registry.add("hDebugDeltaPhiVsDeltaPairPt",
                 "Debug #Delta#phi vs #it{p}_{T,#mu#mu}^{reco}-#it{p}_{T,#mu#mu}^{MC} (pre-kin. cuts)",
                 kTH2F, {axPairPtDeltaDbg, axDeltaPhiReco});
    registry.add("hDebugDeltaPhiVsPairPtReco",
                 "Debug #Delta#phi vs #it{p}_{T,#mu#mu}^{reco} (pre-kin. cuts)",
                 kTH2F, {axPairPtRecoDbg, axDeltaPhiReco});
    registry.add("hDebugDeltaPhiVsPairMassReco",
                 "Debug #Delta#phi vs #it{M}_{#mu#mu}^{reco} (pre-kin. cuts)",
                 kTH2F, {axPairMassRecoDbg, axDeltaPhiReco});

    // --- [Debug-5] Pair mass / pT resolution vs muon pT asymmetry ---
    // A pT-dependent momentum-scale or material-budget bias affects the two
    // legs differently when their pT split is very unequal; this can show up
    // as a pair mass/pT resolution that shifts or widens with
    // (p_{T,mu+} - p_{T,mu-})^{MC}. Filled post-cut, same population as the
    // standard hPairMassResoVsMassMC / hPairPtResoVsPtMC histograms.
    const AxisSpec axDeltaPtMuonsMC{nBinsDeltaPtMuons, -deltaPtMuonsMax, deltaPtMuonsMax,
                                    "p_{T,#mu^{+}}^{MC} - p_{T,#mu^{-}}^{MC} (GeV/c)"};
    const AxisSpec axPairMassResDbg{200, -0.5f, 0.5f,
                                    "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairPtRelResDbg{200, -1.0f, 1.0f,
                                     "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};

    registry.add("hPairMassResoVsDeltaPtMuonsMC",
                 "Debug pair mass resolution vs muon p_{T} asymmetry (post-cut)",
                 kTH2F, {axDeltaPtMuonsMC, axPairMassResDbg});
    registry.add("hPairPtResoVsDeltaPtMuonsMC",
                 "Debug pair p_{T} resolution vs muon p_{T} asymmetry (post-cut)",
                 kTH2F, {axDeltaPtMuonsMC, axPairPtRelResDbg});

    // --- [Debug-6] Pair pT resolution vs softer-leg pT ---
    // The softer (smaller-pT) leg generally has the worse single-track pT
    // resolution (more affected by multiple scattering at low p), so it is
    // likely to dominate the pair pT resolution. Binning by min(pT1,pT2)^MC
    // tests this directly, independent of which leg carries which charge.
    const AxisSpec axSoftMuonPtMC{nBinsP, 0.f, ptTrackMax, "min(p_{T,#mu1},p_{T,#mu2})^{MC} (GeV/c)"};
    registry.add("hPairPtResoVsSoftMuonPtMC",
                 "Debug pair p_{T} resolution vs softer-leg p_{T}^{MC} (post-cut)",
                 kTH2F, {axSoftMuonPtMC, axPairPtRelResDbg});
  }

  // ---------------------------------------------------------------------------
  // [Debug] Opening-angle / rapidity / helicity-angle diagnostic histograms
  // ---------------------------------------------------------------------------
  void initDebugKinematicHistograms()
  {
    const AxisSpec axPairMassResDbg2{200, -0.5f, 0.5f,
                                     "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairPtRelResDbg2{200, -1.0f, 1.0f,
                                      "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};

    // --- [Debug-8] Opening angle #theta / #DeltaR vs mass resolution ---
    // For small opening angles M^2 ~ p_1 p_2 theta^2, so if the mass deficit
    // grows with the *reconstructed* opening angle (or DeltaR), the two legs'
    // angular separation is being measured too small (or too large); if the
    // deficit is flat in angle, the bias sits in the single-track momentum
    // scale instead. Checked with both the 3D opening angle theta and the
    // (eta,phi)-based DeltaR, and with both MC-truth and reco binning.
    const AxisSpec axOpeningAngleMC{nBinsEta, 0.f, openingAngleMax, "#theta_{open}^{MC} (rad)"};
    const AxisSpec axOpeningAngleReco{nBinsEta, 0.f, openingAngleMax, "#theta_{open}^{reco} (rad)"};
    const AxisSpec axDeltaRMC{nBinsEta, 0.f, deltaRMax, "#DeltaR^{MC}"};
    const AxisSpec axDeltaRReco{nBinsEta, 0.f, deltaRMax, "#DeltaR^{reco}"};

    registry.add("hResponseMatrixOpeningAngle", "Pair opening-angle #theta Response Matrix",
                 kTH2F, {axOpeningAngleMC, axOpeningAngleReco});
    registry.add("hResponseMatrixDeltaR", "Pair #DeltaR Response Matrix",
                 kTH2F, {axDeltaRMC, axDeltaRReco});

    registry.add("hPairMassResoVsOpeningAngleMC",
                 "Debug pair mass resolution vs opening angle #theta^{MC} (post-cut)",
                 kTH2F, {axOpeningAngleMC, axPairMassResDbg2});
    registry.add("hPairMassResoVsOpeningAngleReco",
                 "Debug pair mass resolution vs opening angle #theta^{reco} (post-cut)",
                 kTH2F, {axOpeningAngleReco, axPairMassResDbg2});
    registry.add("hPairMassResoVsDeltaRMC",
                 "Debug pair mass resolution vs #DeltaR^{MC} (post-cut)",
                 kTH2F, {axDeltaRMC, axPairMassResDbg2});
    registry.add("hPairMassResoVsDeltaRReco",
                 "Debug pair mass resolution vs #DeltaR^{reco} (post-cut)",
                 kTH2F, {axDeltaRReco, axPairMassResDbg2});

    // --- [Debug-9] Pair p_{T} / mass resolution vs pair rapidity y^{MC} ---
    const AxisSpec axPairRapMCDbg{nBinsEta, pairRapidityMin, pairRapidityMax, "#it{y}_{#mu#mu}^{MC}"};
    registry.add("hPairPtResoVsRapMC",
                 "Debug pair p_{T} resolution vs pair rapidity #it{y}^{MC} (post-cut)",
                 kTH2F, {axPairRapMCDbg, axPairPtRelResDbg2});
    registry.add("hPairMassResoVsRapMC",
                 "Debug pair mass resolution vs pair rapidity #it{y}^{MC} (post-cut)",
                 kTH2F, {axPairRapMCDbg, axPairMassResDbg2});

    // --- [Debug-10] Pair p_{T} / mass resolution vs helicity decay angle cos(theta_HE) ---
    // cos(theta_HE) is the angle of the mu+ in the dimuon rest frame w.r.t.
    // the dimuon's lab-frame flight direction. cos(theta_HE) ~ 0: the two
    // muons split the pair momentum symmetrically in the lab. cos(theta_HE)
    // ~ +/-1: one muon carries most of the pair momentum while the other is
    // soft, close to the single-muon pT/acceptance cut. Comparing the two
    // regimes separates a generic high-p resolution deficit from an
    // acceptance-edge/low-pT-cut distortion.
    const AxisSpec axCosThetaHEMC{100, -1.f, 1.f, "cos(#theta_{HE})^{MC}"};
    const AxisSpec axCosThetaHEReco{100, -1.f, 1.f, "cos(#theta_{HE})^{reco}"};
    registry.add("hResponseMatrixCosThetaHE", "Pair cos(#theta_{HE}) Response Matrix",
                 kTH2F, {axCosThetaHEMC, axCosThetaHEReco});
    registry.add("hPairPtResoVsCosThetaHEMC",
                 "Debug pair p_{T} resolution vs cos(#theta_{HE})^{MC} (post-cut)",
                 kTH2F, {axCosThetaHEMC, axPairPtRelResDbg2});
    registry.add("hPairMassResoVsCosThetaHEMC",
                 "Debug pair mass resolution vs cos(#theta_{HE})^{MC} (post-cut)",
                 kTH2F, {axCosThetaHEMC, axPairMassResDbg2});
  }

  // ---------------------------------------------------------------------------
  // [Debug] Resolution vs number of MFT clusters
  // ---------------------------------------------------------------------------
  void initDebugMFTClusterHistograms()
  {
    // --- [Debug-11] Single-track / pair resolution vs N MFT clusters ---
    // A track fit constrained by fewer MFT hits is less precisely pointed,
    // so both the single-track (q/pT) residual and the pair mass/pT
    // resolution are expected to degrade as the number of MFT clusters used
    // in the fit drops (MFT has 10 sensor layers, so 0-10 hits per track).
    const AxisSpec axNMFTClusters{11, -0.5f, 10.5f, "N_{MFT clusters}"};
    const AxisSpec axQOverPtResDbg3{400, -0.5f, 0.5f,
                                    "(q/p_{T})^{reco}-(q/p_{T})^{MC} (GeV/c)^{-1}"};
    const AxisSpec axPtRelResDbg3{200, -0.5f, 0.5f,
                                  "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};

    registry.add("hDebugNMFTClusters_MuPlus", "Debug N_{MFT clusters} #mu^{+}", kTH1I, {axNMFTClusters});
    registry.add("hDebugNMFTClusters_MuMinus", "Debug N_{MFT clusters} #mu^{-}", kTH1I, {axNMFTClusters});

    registry.add("hDebugPtResoVsNMFTClusters",
                 "Debug track p_{T} resolution vs N_{MFT clusters}",
                 kTH2F, {axNMFTClusters, axPtRelResDbg3});
    registry.add("hDebugQOverPtResoVsNMFTClusters_MuPlus",
                 "Debug (q/p_{T}) residual vs N_{MFT clusters} #mu^{+}",
                 kTH2F, {axNMFTClusters, axQOverPtResDbg3});
    registry.add("hDebugQOverPtResoVsNMFTClusters_MuMinus",
                 "Debug (q/p_{T}) residual vs N_{MFT clusters} #mu^{-}",
                 kTH2F, {axNMFTClusters, axQOverPtResDbg3});

    const AxisSpec axPairMassResDbg3{200, -0.5f, 0.5f,
                                     "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairPtRelResDbg3{200, -1.0f, 1.0f,
                                      "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};
    // Binned by min(N_{MFT clusters}) of the two legs, mirroring the
    // [Debug-6] softer-leg-pT convention: the worse-constrained leg is
    // expected to dominate the pair resolution.
    registry.add("hPairMassResoVsMinNMFTClusters",
                 "Debug pair mass resolution vs min(N_{MFT clusters}) (post-cut)",
                 kTH2F, {axNMFTClusters, axPairMassResDbg3});
    registry.add("hPairPtResoVsMinNMFTClusters",
                 "Debug pair p_{T} resolution vs min(N_{MFT clusters}) (post-cut)",
                 kTH2F, {axNMFTClusters, axPairPtRelResDbg3});
  }

  // ---------------------------------------------------------------------------
  // [Debug] Pair pT / mass vs number of MCH clusters (UDFwdTracksExtra_001::fNClusters)
  // ---------------------------------------------------------------------------
  void initDebugNClustersHistograms()
  {
    // --- [Debug-12] Pair p_{T} / mass vs min(N_{MCH clusters}) ---
    // tr.nClusters() (UDFwdTracksExtra_001::fNClusters) is directly available
    // on ForwardTracks/CompleteFwdTracks, no extra table join needed. Binned
    // by min(N_{MCH clusters}) of the two legs (same worse-constrained-leg
    // logic as [Debug-6]/[Debug-11]). Both the raw reco pair p_{T}/mass
    // distribution and the p_{T}/mass resolution (reco-MC) are filled, so a
    // shift in the raw distribution vs a widening/shift in the resolution can
    // be told apart.
    const AxisSpec axNClusters{11, -0.5f, 10.5f, "N_{MCH clusters}"};
    const AxisSpec axPairPtRecoDbg3{nBinsPt, 0.f, ptMax, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axPairMassRecoDbg3{nBinsMass, massAxisMin, massAxisMax, "#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassResDbg4{200, -0.5f, 0.5f,
                                     "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairPtRelResDbg4{200, -1.0f, 1.0f,
                                      "(#it{p}_{T}^{reco} - #it{p}_{T}^{MC}) / #it{p}_{T}^{MC}"};

    registry.add("hPairPtRecoVsMinNClusters",
                 "Debug pair p_{T}^{reco} vs min(N_{MCH clusters}) (post-cut)",
                 kTH2F, {axNClusters, axPairPtRecoDbg3});
    registry.add("hPairMassRecoVsMinNClusters",
                 "Debug pair mass^{reco} vs min(N_{MCH clusters}) (post-cut)",
                 kTH2F, {axNClusters, axPairMassRecoDbg3});
    registry.add("hPairPtResoVsMinNClusters",
                 "Debug pair p_{T} resolution vs min(N_{MCH clusters}) (post-cut)",
                 kTH2F, {axNClusters, axPairPtRelResDbg4});
    registry.add("hPairMassResoVsMinNClusters",
                 "Debug pair mass resolution vs min(N_{MCH clusters}) (post-cut)",
                 kTH2F, {axNClusters, axPairMassResDbg4});
  }

  // ---------------------------------------------------------------------------
  // [Debug] Softer-leg total momentum vs pair mass
  // ---------------------------------------------------------------------------
  void initDebugSoftLegMomentumHistograms()
  {
    // --- [Debug-13] Softer-leg total momentum p vs pair mass ---
    // M^2 ~ p_1 p_2 theta^2: [Debug-8] tested the angular half of this
    // relation; this tests the momentum half. The "softer leg" is the muon
    // with the smaller total momentum p^{MC} of the pair (not pT).
    //   hPairMassRecoVsSoftLegPMC: raw reco pair mass vs p^{MC}_{soft leg}.
    //     A mass deficit that only appears at low soft-leg momentum, while
    //     the raw mass otherwise sits at the nominal value, points to an
    //     acceptance/low-p edge effect.
    //   hPairMassResoVsSoftLegPReso: pair mass resolution (reco-MC) vs the
    //     soft leg's own momentum resolution. A direct correlation here means
    //     the soft leg's momentum underestimation is itself what drives the
    //     pair mass deficit (as opposed to an angular mismeasurement).
    const AxisSpec axSoftLegPMC{nBinsP, 0.f, pTrackMax, "p_{soft leg}^{MC} (GeV/#it{c})"};
    const AxisSpec axSoftLegPRelRes{200, -0.5f, 0.5f,
                                    "(p_{soft leg}^{reco} - p_{soft leg}^{MC}) / p_{soft leg}^{MC}"};
    const AxisSpec axPairMassRecoDbg5{nBinsMass, massAxisMin, massAxisMax, "#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassResDbg5{200, -0.5f, 0.5f,
                                     "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};

    registry.add("hPairMassRecoVsSoftLegPMC",
                 "Debug pair mass^{reco} vs softer-leg p^{MC} (post-cut)",
                 kTH2F, {axSoftLegPMC, axPairMassRecoDbg5});
    registry.add("hPairMassResoVsSoftLegPReso",
                 "Debug pair mass resolution vs softer-leg p resolution (post-cut)",
                 kTH2F, {axSoftLegPRelRes, axPairMassResDbg5});
  }

  // ---------------------------------------------------------------------------
  // [Debug-14] MC-truth direct comparison: opening-angle / single-muon-pT
  // residuals vs pair pT, and mass-deficit-tail gating
  // ---------------------------------------------------------------------------
  void initDebugMcTruthDirectComparisonHistograms()
  {
    // --- (1) Delta-theta = theta^{reco} - theta^{MC} vs pair p_{T}^{reco} ---
    // Filled pre-kinematic-cuts (same convention as [Debug-7]) so the full
    // p_{T,mu mu}^{reco} range is visible, not just the <pairPtMax signal
    // region. If an underestimated opening angle (Delta-theta<0) is what
    // inflates the reco pair pT, Delta-theta<0 events should extend to
    // higher pair pT than Delta-theta~0 events.
    const AxisSpec axDeltaOpeningAngle{400, -0.02f, 0.02f, "#theta^{reco}-#theta^{MC} (rad)"};
    const AxisSpec axPairPtRecoDbg6{nBinsPt, 0.f, 2.0, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};
    registry.add("hDebugDeltaOpeningAngleVsPairPtReco",
                 "Debug #theta^{reco}-#theta^{MC} vs #it{p}_{T,#mu#mu}^{reco} (pre-kin. cuts)",
                 kTH2F, {axPairPtRecoDbg6, axDeltaOpeningAngle});

    // --- (2) Delta-pT_mu = pT_mu^{reco} - pT_mu^{MC} (either leg) vs pair p_{T}^{reco} ---
    // Filled once per leg per pair (pre-kinematic-cuts). Separates scenario A
    // (one leg's pT is overestimated -> pair pT rises with that leg's
    // Delta-pT) from scenario B (both legs' Delta-pT stay small while the
    // pair pT still rises, pointing back at the opening-angle mismeasurement
    // in (1) instead).
    const AxisSpec axDeltaPtMu{400, -2.0f, 2.0f, "#it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} (GeV/#it{c})"};
    registry.add("hDebugDeltaPtMuVsPairPtReco",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs #it{p}_{T,#mu#mu}^{reco} (pre-kin. cuts, both legs)",
                 kTH2F, {axPairPtRecoDbg6, axDeltaPtMu});

    // --- (3) Mass-deficit-tail gate ---
    // Post-cut (same accepted sample as [Debug-8] etc.): compare the
    // Delta-theta and Delta-pT_mu distributions for pairs whose reco mass
    // falls in the low-mass tail [massTailMin,massTailMax] against the
    // inclusive (all accepted pairs) distributions, to see whether the
    // mass-deficit tail is populated by a distinct (angle- or
    // momentum-mismeasured) sub-population.
    registry.add("hDebugDeltaOpeningAngle_Inclusive",
                 "Debug #theta^{reco}-#theta^{MC}, all accepted pairs (post-cut)",
                 kTH1F, {axDeltaOpeningAngle});
    registry.add("hDebugDeltaOpeningAngle_MassTail",
                 "Debug #theta^{reco}-#theta^{MC}, mass-deficit tail only (post-cut)",
                 kTH1F, {axDeltaOpeningAngle});
    registry.add("hDebugDeltaPtMu_Inclusive",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axDeltaPtMu});
    registry.add("hDebugDeltaPtMu_MassTail",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axDeltaPtMu});

    // --- (4) Single-muon kinematics + R at absorber end: mass-deficit-tail
    // gate ---
    // Same Inclusive/MassTail comparison as (3), but for the raw single-muon
    // reco pT/eta/phi and R_abs (both legs), to see whether the tail muons
    // sit in a distinct region of single-track phase space (e.g. low pT,
    // eta/phi acceptance edge) or a distinct R_abs range (a cheap proxy for
    // the track's (x,y) position at the absorber end, since UDFwdTracksProp
    // is not joined here).
    const AxisSpec axSingleMuPtReco{nBinsPt, 0.f, ptTrackMax, "#it{p}_{T,#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axSingleMuEtaReco{nBinsEta, etaMin, etaMax, "#eta_{#mu}^{reco}"};
    const AxisSpec axSingleMuPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi_{#mu}^{reco} (rad)"};
    const AxisSpec axSingleMuRAbs{200, rAbsMin, rAbsMax, "R_{abs,#mu}^{reco} (cm)"};

    registry.add("hDebugPtMu_Inclusive",
                 "Debug #it{p}_{T,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuPtReco});
    registry.add("hDebugPtMu_MassTail",
                 "Debug #it{p}_{T,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuPtReco});

    // MC-truth comparison for hDebugPtMu, plus the reco-vs-MC response
    // itself (not just the resolution difference already covered by
    // hDebugDeltaPtMu*), so the reco and MC-truth pT shapes -- and how
    // reco tracks pT^{MC} across the full range -- can be read directly.
    const AxisSpec axSingleMuPtMC{nBinsPt, 0.f, ptTrackMax, "#it{p}_{T,#mu}^{MC} (GeV/#it{c})"};

    registry.add("hDebugPtMuMC_Inclusive",
                 "Debug #it{p}_{T,#mu}^{MC}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuPtMC});
    registry.add("hDebugPtMuMC_MassTail",
                 "Debug #it{p}_{T,#mu}^{MC}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuPtMC});

    registry.add("hDebugPtMuRecoVsMC_Inclusive",
                 "Debug #it{p}_{T,#mu}^{reco} vs #it{p}_{T,#mu}^{MC}, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axSingleMuPtMC, axSingleMuPtReco});
    registry.add("hDebugPtMuRecoVsMC_MassTail",
                 "Debug #it{p}_{T,#mu}^{reco} vs #it{p}_{T,#mu}^{MC}, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axSingleMuPtMC, axSingleMuPtReco});

    registry.add("hDebugEtaMu_Inclusive",
                 "Debug #eta_{#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuEtaReco});
    registry.add("hDebugEtaMu_MassTail",
                 "Debug #eta_{#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuEtaReco});

    registry.add("hDebugPhiMu_Inclusive",
                 "Debug #phi_{#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuPhiReco});
    registry.add("hDebugPhiMu_MassTail",
                 "Debug #phi_{#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuPhiReco});

    registry.add("hDebugRAbsMu_Inclusive",
                 "Debug R_{abs,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuRAbs});
    registry.add("hDebugRAbsMu_MassTail",
                 "Debug R_{abs,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuRAbs});

    // --- (5) Track fit quality: global chi2 and MCH-MFT match chi2 ---
    // Same Inclusive/MassTail gate as (3)/(4). pT already showed a
    // correlation with the tail; chi2/chi2MatchMCHMFT are checked next since
    // a track fit that is a worse match to its cluster hits is a direct
    // mismeasurement handle, and worth checking against pT as a possible
    // confound (soft tracks tend to fit worse anyway).
    const AxisSpec axSingleMuChi2{200, 0.f, maxChi2, "#chi^{2}_{#mu}^{reco}"};
    const AxisSpec axSingleMuChi2MatchMCHMFT{200, 0.f, maxChi2MatchMCHMFT, "#chi^{2}_{MCH-MFT match,#mu}^{reco}"};

    registry.add("hDebugChi2Mu_Inclusive",
                 "Debug #chi^{2}_{#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuChi2});
    registry.add("hDebugChi2Mu_MassTail",
                 "Debug #chi^{2}_{#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuChi2});

    registry.add("hDebugChi2MatchMCHMFTMu_Inclusive",
                 "Debug #chi^{2}_{MCH-MFT match,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuChi2MatchMCHMFT});
    registry.add("hDebugChi2MatchMCHMFTMu_MassTail",
                 "Debug #chi^{2}_{MCH-MFT match,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuChi2MatchMCHMFT});

    // --- (6) MCH cluster count and pDCA: mass-deficit-tail gate ---
    // Same Inclusive/MassTail gate as (3)-(5). N_{MCH clusters} is already
    // used per-leg in [Debug-12] (as the min over both legs); pDca is the
    // p*DCA quality cut applied in passRecoTrackCuts (threshold depends on
    // rAbs via the 26.5 cm boundary, so it is not a pure pT proxy).
    const AxisSpec axSingleMuNClusters{11, -0.5f, 20.5f, "N_{MCH clusters,#mu}^{reco}"};
    const AxisSpec axSingleMuPDca{200, 0.f, 350.f, "p#times DCA_{#mu}^{reco} (cm#upointGeV/#it{c})"};

    registry.add("hDebugNClustersMu_Inclusive",
                 "Debug N_{MCH clusters,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1I, {axSingleMuNClusters});
    registry.add("hDebugNClustersMu_MassTail",
                 "Debug N_{MCH clusters,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1I, {axSingleMuNClusters});

    registry.add("hDebugPDcaMu_Inclusive",
                 "Debug p#times DCA_{#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuPDca});
    registry.add("hDebugPDcaMu_MassTail",
                 "Debug p#times DCA_{#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuPDca});

    // --- (7) N_{MFT clusters}: mass-deficit-tail gate ---
    // Complements chi2MatchMCHMFT in (5): a low match chi2 can still ride on
    // very few MFT clusters (few points to disagree with), so cluster count
    // and match chi2 need to be read together to judge whether the MFT-MCH
    // matching itself is degraded for tail muons. nMFTClusters1/2 are already
    // passed into this function (used combined as minNMFTClusters in
    // [Debug-11]); here each leg is filled separately and gated as in (3)-(6).
    const AxisSpec axSingleMuNMFTClusters{11, -0.5f, 10.5f, "N_{MFT clusters,#mu}^{reco}"};

    registry.add("hDebugNMFTClustersMu_Inclusive",
                 "Debug N_{MFT clusters,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH1I, {axSingleMuNMFTClusters});
    registry.add("hDebugNMFTClustersMu_MassTail",
                 "Debug N_{MFT clusters,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH1I, {axSingleMuNMFTClusters});

    // --- (8) MC-truth MCH-MFT match correctness (mcMask fake bit) ---
    // Unlike (5)-(7), which only check the *quality* of the match (chi2,
    // cluster count), this checks *correctness* using MC truth: the AOD
    // producer sets bit 7 (0x80) of UDMcFwdTrackLabels::mcMask when the
    // MFTMCH-source global track's label is fake, i.e. when the MCH leg and
    // the MFT leg of the geometric/kinematic match do not trace back to the
    // same true MC particle (a wrong combinatorial pairing masquerading as
    // one muon). Bit 6 (0x40) similarly flags a match to pure noise. Both are
    // available directly on tr1/tr2 (CompleteFwdTracks already joins
    // UDMcFwdTrackLabels), no extra table join needed. Same Inclusive/
    // MassTail gate as (3)-(7).
    const AxisSpec axSingleMuMcMaskFlag{2, -0.5f, 1.5f, "flag (0/1)"};

    registry.add("hDebugMcMaskFakeMu_Inclusive",
                 "Debug MCH-MFT global label isFake (mcMask bit 7), all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuMcMaskFlag});
    registry.add("hDebugMcMaskFakeMu_MassTail",
                 "Debug MCH-MFT global label isFake (mcMask bit 7), mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuMcMaskFlag});

    registry.add("hDebugMcMaskNoiseMu_Inclusive",
                 "Debug MCH-MFT global label isNoise (mcMask bit 6), all accepted pairs (post-cut, both legs)",
                 kTH1F, {axSingleMuMcMaskFlag});
    registry.add("hDebugMcMaskNoiseMu_MassTail",
                 "Debug MCH-MFT global label isNoise (mcMask bit 6), mass-deficit tail only (post-cut, both legs)",
                 kTH1F, {axSingleMuMcMaskFlag});

    // --- (9) Confound check: chi2 / chi2MatchMCHMFT vs pT ---
    // (5) showed chi2 and chi2MatchMCHMFT correlate with the mass-deficit
    // tail, same as pT did earlier. Softer tracks generically fit worse, so
    // that correlation could just be riding on the pT one instead of being
    // an independent effect. These 2D histograms (both legs, reusing the
    // pT/chi2 axes from (4)/(5)) let chi2 vs pT be read directly within the
    // MassTail-gated sample and compared to the Inclusive sample's own
    // chi2-pT correlation shape.
    registry.add("hDebugChi2VsPtMu_Inclusive",
                 "Debug #chi^{2}_{#mu}^{reco} vs #it{p}_{T,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axSingleMuPtReco, axSingleMuChi2});
    registry.add("hDebugChi2VsPtMu_MassTail",
                 "Debug #chi^{2}_{#mu}^{reco} vs #it{p}_{T,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axSingleMuPtReco, axSingleMuChi2});

    registry.add("hDebugChi2MatchMCHMFTVsPtMu_Inclusive",
                 "Debug #chi^{2}_{MCH-MFT match,#mu}^{reco} vs #it{p}_{T,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axSingleMuPtReco, axSingleMuChi2MatchMCHMFT});
    registry.add("hDebugChi2MatchMCHMFTVsPtMu_MassTail",
                 "Debug #chi^{2}_{MCH-MFT match,#mu}^{reco} vs #it{p}_{T,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axSingleMuPtReco, axSingleMuChi2MatchMCHMFT});
  }

  // ---------------------------------------------------------------------------
  // [Debug-15] Verification method B: mass resolution normalized by true pair pT
  // ---------------------------------------------------------------------------
  void initDebugNormalizedMassResoHistograms()
  {
    // The pT^2-resolution cross term 2*p_{T,true}*Sigma(eps_x) scales with
    // the event's true pair pT, so it dominates the raw (M_reco-M_MC) spread
    // for incoherent events with a wide p_{T,true} spread. Dividing the mass
    // residual by p_{T,mu mu}^{true} removes this scale dependence; if the
    // cross term is indeed the driver, the normalized residual should be
    // flat vs p_{T,true} where the raw one is not.
    const AxisSpec axPairPtMCDbg7{nBinsPt, 0.f, 2.0f, "p_{T,#mu#mu}^{true} (GeV/#it{c})"};
    const AxisSpec axPairMassResRaw{200, -0.5f, 0.5f,
                                    "#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC} (GeV/#it{c}^{2})"};
    const AxisSpec axPairMassResNorm{400, -20.0f, 20.0f,
                                     "(#it{M}_{#mu#mu}^{reco} - #it{M}_{#mu#mu}^{MC}) / #it{p}_{T,#mu#mu}^{true}"};

    registry.add("hPairMassResoVsPtMC",
                 "Debug pair mass resolution (raw) vs #it{p}_{T,#mu#mu}^{true} (post-cut)",
                 kTH2F, {axPairPtMCDbg7, axPairMassResRaw});
    registry.add("hPairMassResoNormVsPtMC",
                 "Debug pair mass resolution / #it{p}_{T,#mu#mu}^{true} (Verification B) (post-cut)",
                 kTH2F, {axPairPtMCDbg7, axPairMassResNorm});
  }

  // ---------------------------------------------------------------------------
  // [Debug-16] Delta-pT_mu vs single-muon reco pT (Inclusive / MassTail)
  // ---------------------------------------------------------------------------
  void initDebugDeltaPtMuVsRecoPtHistograms()
  {
    // hDebugDeltaPtMu_Inclusive/_MassTail ([Debug-14](3)) are 1D and only
    // show the overall shape of the momentum-underestimated sub-population
    // in the mass-deficit tail. This adds the same split but 2D against the
    // muon's own reco pT, to localise *where* in reco-pT space that
    // sub-population sits (e.g. concentrated near the pTmin acceptance edge
    // vs spread across the full range).
    const AxisSpec axMuPtReco{nBinsP, 0.f, ptTrackMax, "p_{T,#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axDeltaPtMuDbg8{400, -2.0f, 2.0f, "p_{T,#mu}^{reco}-p_{T,#mu}^{MC} (GeV/#it{c})"};

    registry.add("hDebugDeltaPtMuVsRecoPt_Inclusive",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs #it{p}_{T,#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axMuPtReco, axDeltaPtMuDbg8});
    registry.add("hDebugDeltaPtMuVsRecoPt_MassTail",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs #it{p}_{T,#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axMuPtReco, axDeltaPtMuDbg8});
  }

  // ---------------------------------------------------------------------------
  // [Debug-17] Delta-pT_mu vs single-muon phi (reco), and split by charge
  // ---------------------------------------------------------------------------
  void initDebugDeltaPtMuVsPhiAndSignHistograms()
  {
    // (1) vs phi^{reco}: looks for a detector-sector dependence (MFT disk
    // boundary, MCH quadrant) in the momentum-underestimated sub-population,
    // same Inclusive/MassTail gate as [Debug-14](3)/[Debug-16].
    const AxisSpec axMuPhiReco{nBinsPhi, -TMath::Pi(), TMath::Pi(), "#phi_{#mu}^{reco} (rad)"};
    const AxisSpec axDeltaPtMuDbg9{400, -2.0f, 2.0f, "p_{T,#mu}^{reco}-p_{T,#mu}^{MC} (GeV/#it{c})"};

    registry.add("hDebugDeltaPtMuVsPhiMu_Inclusive",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs #phi_{#mu}^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axMuPhiReco, axDeltaPtMuDbg9});
    registry.add("hDebugDeltaPtMuVsPhiMu_MassTail",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs #phi_{#mu}^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axMuPhiReco, axDeltaPtMuDbg9});

    // (2) split by charge: if the momentum-underestimated sub-population is
    // concentrated in mu+ or mu- specifically, that points back at a
    // charge-dependent bias (e.g. an MFT/MCH rotation/twist, as in [Debug-1]).
    registry.add("hDebugDeltaPtMu_Inclusive_MuPlus",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC}, #mu^{+}, all accepted pairs (post-cut)",
                 kTH1F, {axDeltaPtMuDbg9});
    registry.add("hDebugDeltaPtMu_Inclusive_MuMinus",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC}, #mu^{-}, all accepted pairs (post-cut)",
                 kTH1F, {axDeltaPtMuDbg9});
    registry.add("hDebugDeltaPtMu_MassTail_MuPlus",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC}, #mu^{+}, mass-deficit tail only (post-cut)",
                 kTH1F, {axDeltaPtMuDbg9});
    registry.add("hDebugDeltaPtMu_MassTail_MuMinus",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC}, #mu^{-}, mass-deficit tail only (post-cut)",
                 kTH1F, {axDeltaPtMuDbg9});
  }

  // ---------------------------------------------------------------------------
  // [Debug-18] Approximate single-muon (x,y) position at the absorber end
  // ---------------------------------------------------------------------------
  void initDebugXYPositionHistograms()
  {
    // True fitted (x,y) position isn't available in the current UD skimmed
    // tables (UDFwdTracksProp/UDFwdTracksCovProp are not produced by
    // upcCandProducerGlobalMuon.cxx). Approximate the transverse position at
    // the absorber end from already-available quantities:
    //   x ~ R_abs * cos(phi^{reco}), y ~ R_abs * sin(phi^{reco})
    // using the reco momentum's azimuthal angle as a stand-in for the
    // position's azimuthal angle (exact for a straight line, off by the
    // magnetic bending angle otherwise -- fine for a coarse spatial map).
    const AxisSpec axXApprox{200, -100.f, 100.f, "x_{abs}^{approx} = R_{abs}cos(#phi^{reco}) (cm)"};
    const AxisSpec axYApprox{200, -100.f, 100.f, "y_{abs}^{approx} = R_{abs}sin(#phi^{reco}) (cm)"};

    registry.add("hDebugXYMu_Inclusive",
                 "Debug approx (x,y) at absorber end, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axXApprox, axYApprox});
    registry.add("hDebugXYMu_MassTail",
                 "Debug approx (x,y) at absorber end, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axXApprox, axYApprox});
    registry.add("hDebugXYMu_AnomalousPt",
                 "Debug approx (x,y) at absorber end, legs with anomalously underestimated p_{T} only (post-cut)",
                 kTH2F, {axXApprox, axYApprox});
  }

  // ---------------------------------------------------------------------------
  // [Debug-19] Delta-pT_mu vs pDca, and event vertex position
  // ---------------------------------------------------------------------------
  void initDebugDcaVertexHistograms()
  {
    // (1) vs pDca^{reco}: hDebugPDcaMu_Inclusive/_MassTail ([Debug-14](6))
    // are 1D; this correlates the momentum residual directly with pDca to
    // see whether the momentum-underestimated sub-population also carries
    // an anomalous DCA (both are outputs of the same track fit).
    const AxisSpec axPDcaMu{200, 0.f, 1000.f, "pDCA_{#mu}^{reco} (cm GeV/c)"};
    const AxisSpec axDeltaPtMuDbg10{400, -2.0f, 2.0f, "p_{T,#mu}^{reco}-p_{T,#mu}^{MC} (GeV/#it{c})"};

    registry.add("hDebugDeltaPtMuVsPDca_Inclusive",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs pDCA^{reco}, all accepted pairs (post-cut, both legs)",
                 kTH2F, {axPDcaMu, axDeltaPtMuDbg10});
    registry.add("hDebugDeltaPtMuVsPDca_MassTail",
                 "Debug #it{p}_{T,#mu}^{reco}-#it{p}_{T,#mu}^{MC} vs pDCA^{reco}, mass-deficit tail only (post-cut, both legs)",
                 kTH2F, {axPDcaMu, axDeltaPtMuDbg10});

    // (2) Event (collision) vertex position: already available via
    // UDCollisions (posX/posY/posZ), joined as CandidatesFwd but previously
    // unused. A poorly-reconstructed vertex feeds directly into the MCH
    // track fit's DCA/pDca calculation, so a vertex-position anomaly could
    // explain a pDca/pT anomaly. Same Inclusive/MassTail/AnomalousPt split
    // as [Debug-18] (one fill per pair for Inclusive/MassTail, one fill per
    // anomalous leg found for AnomalousPt).
    const AxisSpec axVtxXY{200, -1.f, 1.f, "vtx (cm)"};
    const AxisSpec axVtxZ{200, -20.f, 20.f, "vtx_{z} (cm)"};

    registry.add("hDebugVertexXY_Inclusive",
                 "Debug event vertex (x,y), all accepted pairs (post-cut)",
                 kTH2F, {axVtxXY, axVtxXY});
    registry.add("hDebugVertexXY_MassTail",
                 "Debug event vertex (x,y), mass-deficit tail only (post-cut)",
                 kTH2F, {axVtxXY, axVtxXY});
    registry.add("hDebugVertexXY_AnomalousPt",
                 "Debug event vertex (x,y), pairs with an anomalously underestimated leg (post-cut)",
                 kTH2F, {axVtxXY, axVtxXY});

    registry.add("hDebugVertexZ_Inclusive",
                 "Debug event vertex z, all accepted pairs (post-cut)",
                 kTH1F, {axVtxZ});
    registry.add("hDebugVertexZ_MassTail",
                 "Debug event vertex z, mass-deficit tail only (post-cut)",
                 kTH1F, {axVtxZ});
    registry.add("hDebugVertexZ_AnomalousPt",
                 "Debug event vertex z, pairs with an anomalously underestimated leg (post-cut)",
                 kTH1F, {axVtxZ});
  }

  // ---------------------------------------------------------------------------
  // [Debug-20] Muon pT asymmetry (reco) vs pair pT (reco)
  // ---------------------------------------------------------------------------
  void initDebugDeltaPtMuonsRecoHistograms()
  {
    // Exclusive UPC production predicts p_{T,mu1}^{true} ~ -p_{T,mu2}^{true}
    // (near back-to-back in the transverse plane), so a pair with strongly
    // asymmetric leg pT's is expected to show a larger vector-sum (pair) pT
    // than a balanced pair. Tests that kinematic relationship directly with
    // reco quantities only (no MC truth), same Inclusive/MassTail gate as
    // [Debug-16]/[Debug-17].
    const AxisSpec axPairPtRecoDbg8{nBinsPt, 0.f, 2.0f, "#it{p}_{T,#mu#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axDeltaPtMuonsReco{200, 0.f, 5.f, "|p_{T,#mu1}^{reco}-p_{T,#mu2}^{reco}| (GeV/#it{c})"};

    registry.add("hDeltaPtMuonsRecoVsPairPtReco_Inclusive",
                 "Debug |p_{T,#mu1}^{reco}-p_{T,#mu2}^{reco}| vs #it{p}_{T,#mu#mu}^{reco}, all accepted pairs (post-cut)",
                 kTH2F, {axPairPtRecoDbg8, axDeltaPtMuonsReco});
    registry.add("hDeltaPtMuonsRecoVsPairPtReco_MassTail",
                 "Debug |p_{T,#mu1}^{reco}-p_{T,#mu2}^{reco}| vs #it{p}_{T,#mu#mu}^{reco}, mass-deficit tail only (post-cut)",
                 kTH2F, {axPairPtRecoDbg8, axDeltaPtMuonsReco});
  }

  // ---------------------------------------------------------------------------
  // [Debug-21] Muon pair track-time difference
  // ---------------------------------------------------------------------------
  void initDebugDeltaTrackTimeHistograms()
  {
    // trackTime()/trackTimeRes() are already on UDFwdTracks (no extra table
    // join needed) but were unused until now. If the two legs of a "pair"
    // are actually mismatched to different bunch crossings (e.g. one leg is
    // a pileup/neighbouring-BC track wrongly combined with the other), that
    // should show up as an anomalously large |trackTime_mu1-trackTime_mu2|,
    // which would also explain a wrong vertex/cluster assignment feeding
    // into the momentum anomaly investigated in [Debug-14]-[Debug-20]. Step
    // (1): plain Inclusive/MassTail comparison; step (2) (vs Delta-pT_mu /
    // pair pT) to follow once a trend is seen here.
    const AxisSpec axDeltaTrackTime{200, 0.f, 50.f, "|trackTime_{#mu1}-trackTime_{#mu2}| (ns)"};

    registry.add("hDebugDeltaTrackTime_Inclusive",
                 "Debug |trackTime_{#mu1}-trackTime_{#mu2}|, all accepted pairs (post-cut)",
                 kTH1F, {axDeltaTrackTime});
    registry.add("hDebugDeltaTrackTime_MassTail",
                 "Debug |trackTime_{#mu1}-trackTime_{#mu2}|, mass-deficit tail only (post-cut)",
                 kTH1F, {axDeltaTrackTime});
  }

  // ---------------------------------------------------------------------------
  // [Debug-22] Pair mass vs single-muon pT / eta (continuous, both legs)
  // ---------------------------------------------------------------------------
  void initDebugMassVsSingleMuonKinematicsHistograms()
  {
    // Unlike the massTailMin/massTailMax binary gate used throughout
    // [Debug-14] etc., this correlates the reco pair mass directly and
    // continuously with each leg's own reco pT/eta, to see whether a mass
    // distortion (not just the low-mass tail specifically) tracks a
    // particular single-muon pT or eta region (e.g. low pT, or the eta
    // acceptance edge) across the whole mass range.
    const AxisSpec axSingleMuPtRecoDbg{nBinsPt, 0.f, ptTrackMax, "#it{p}_{T,#mu}^{reco} (GeV/#it{c})"};
    const AxisSpec axSingleMuEtaRecoDbg{nBinsEta, etaMin, etaMax, "#eta_{#mu}^{reco}"};
    const AxisSpec axPairMassRecoDbg6{nBinsMass, massAxisMin, massAxisMax, "#it{M}_{#mu#mu}^{reco} (GeV/#it{c}^{2})"};

    registry.add("hDebugPairMassRecoVsPtMu",
                 "Debug pair mass^{reco} vs #it{p}_{T,#mu}^{reco} (post-cut, both legs)",
                 kTH2F, {axSingleMuPtRecoDbg, axPairMassRecoDbg6});
    registry.add("hDebugPairMassRecoVsEtaMu",
                 "Debug pair mass^{reco} vs #eta_{#mu}^{reco} (post-cut, both legs)",
                 kTH2F, {axSingleMuEtaRecoDbg, axPairMassRecoDbg6});
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
  void fillSingleTrackResolution(const TTrack& tr, const TMcParticle& mc, int nMFTClusters)
  {
    TLorentzVector vReco, vMC;
    vReco.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    vMC.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);

    registry.fill(HIST("hResoTrkPhiMC"), vMC.Phi());
    registry.fill(HIST("hResoTrkPhiReco"), vReco.Phi());
    registry.fill(HIST("hResoTrkEtaMC"), vMC.Eta());
    registry.fill(HIST("hResoTrkEtaReco"), vReco.Eta());

    registry.fill(HIST("hResoTrkPMC"), vMC.P());
    registry.fill(HIST("hResoTrkPReco"), vReco.P());
    registry.fill(HIST("hResoTrkPtMC"), vMC.Pt());
    registry.fill(HIST("hResoTrkPtReco"), vReco.Pt());
    registry.fill(HIST("hResoTrkPxMC"), vMC.Px());
    registry.fill(HIST("hResoTrkPxReco"), vReco.Px());
    registry.fill(HIST("hResoTrkPyMC"), vMC.Py());
    registry.fill(HIST("hResoTrkPyReco"), vReco.Py());
    registry.fill(HIST("hResoTrkPzMC"), vMC.Pz());
    registry.fill(HIST("hResoTrkPzReco"), vReco.Pz());

    const float phiRes = TVector2::Phi_mpi_pi(vReco.Phi() - vMC.Phi());
    registry.fill(HIST("hResoPhi"), vMC.Phi(), phiRes);
    if (tr.sign() > 0)
      registry.fill(HIST("hResoPhi_MuPlus"), vMC.Phi(), phiRes);
    else
      registry.fill(HIST("hResoPhi_MuMinus"), vMC.Phi(), phiRes);
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

    // [Debug-11] Single-track pT resolution vs N MFT clusters used in the fit
    if (vMC.Pt() > 0.f)
      registry.fill(HIST("hDebugPtResoVsNMFTClusters"), nMFTClusters, (vReco.Pt() - vMC.Pt()) / vMC.Pt());
  }

  // ---------------------------------------------------------------------------
  // [Debug-1] (q/pT) residual per charge
  // ---------------------------------------------------------------------------
  template <typename TTrack, typename TMcParticle>
  void fillDebugQOverPtResidual(const TTrack& tr, const TMcParticle& mc, int nMFTClusters)
  {
    TLorentzVector vReco, vMC;
    vReco.SetXYZM(tr.px(), tr.py(), tr.pz(), mMu);
    vMC.SetXYZM(mc.px(), mc.py(), mc.pz(), mMu);

    if (vReco.Pt() <= 0.f || vMC.Pt() <= 0.f)
      return;

    // q/pT is the track fitter's linear curvature variable.
    // A MFT/MCH rotation/twist adds a spurious bending:
    //   mu+ (q=+1): curvature over-estimated  -> residual > 0
    //   mu- (q=-1): curvature under-estimated -> residual < 0
    const float q        = static_cast<float>(tr.sign());
    const float residual = q / vReco.Pt() - q / vMC.Pt();

    if (tr.sign() > 0) {
      registry.fill(HIST("hDebugQOverPtResidual_MuPlus"), residual);
      registry.fill(HIST("hDebugQOverPtResoVsPtMC_MuPlus"), vMC.Pt(), residual);
      registry.fill(HIST("hDebugQOverPtResoVsPMC_MuPlus"), vMC.P(), residual);
      registry.fill(HIST("hDebugQOverPtResoVsEtaMC_MuPlus"), vMC.Eta(), residual);
      // [Debug-11] (q/pT) residual vs N MFT clusters used in the fit
      registry.fill(HIST("hDebugNMFTClusters_MuPlus"), nMFTClusters);
      registry.fill(HIST("hDebugQOverPtResoVsNMFTClusters_MuPlus"), nMFTClusters, residual);
    } else {
      registry.fill(HIST("hDebugQOverPtResidual_MuMinus"), residual);
      registry.fill(HIST("hDebugQOverPtResoVsPtMC_MuMinus"), vMC.Pt(), residual);
      registry.fill(HIST("hDebugQOverPtResoVsPMC_MuMinus"), vMC.P(), residual);
      registry.fill(HIST("hDebugQOverPtResoVsEtaMC_MuMinus"), vMC.Eta(), residual);
      // [Debug-11] (q/pT) residual vs N MFT clusters used in the fit
      registry.fill(HIST("hDebugNMFTClusters_MuMinus"), nMFTClusters);
      registry.fill(HIST("hDebugQOverPtResoVsNMFTClusters_MuMinus"), nMFTClusters, residual);
    }
  }

  template <typename TTrack, typename TMcParticle>
  void fillPairAnalysis(const TTrack& tr1, const TMcParticle& mc1,
                        const TTrack& tr2, const TMcParticle& mc2,
                        int nMFTClusters1, int nMFTClusters2,
                        float vtxX, float vtxY, float vtxZ)
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

    // mu+ / mu- identified once, reused by the Delta(Delta-phi), opening-angle
    // and helicity-angle diagnostics below.
    const TLorentzVector& recoVecPos = (tr1.sign() > 0) ? recoVec1 : recoVec2;
    const TLorentzVector& recoVecNeg = (tr1.sign() > 0) ? recoVec2 : recoVec1;
    const TLorentzVector& mcVecPos   = (tr1.sign() > 0) ? mcVec1   : mcVec2;
    const TLorentzVector& mcVecNeg   = (tr1.sign() > 0) ? mcVec2   : mcVec1;

    // Opening angle theta (MC and reco), hoisted here so it is available to
    // both the pre-cut [Debug-14] check and the post-cut [Debug-8] check.
    const float openingAngleMC = mcVec1.Vect().Angle(mcVec2.Vect());
    const float openingAngleReco = recoVec1.Vect().Angle(recoVec2.Vect());

    const float PtMCPos = mcVecPos.Pt();
    const float PtMCNeg = mcVecNeg.Pt();
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

    // --- [Debug-2] Delta(Delta-phi) bimodal (filled here, before kinematic
    //     cuts, to maximise statistics for the alignment diagnostic) ---
    // Identify mu+ (sign>0) and mu- (sign<0) among the two tracks.
    // A global phi rotation by delta shifts phi(mu+) by +delta and
    // phi(mu-) by -delta (or vice versa), so the difference-of-differences
    //   (phi_{mu+}-phi_{mu-})^{reco} - (phi_{mu+}-phi_{mu-})^{MC}
    // equals +/-2*delta and forms a clean bimodal, free from pT smearing.
    {
      const float deltaPhiReco = TVector2::Phi_mpi_pi(recoVecPos.Phi() - recoVecNeg.Phi());
      const float deltaPhiMC   = TVector2::Phi_mpi_pi(mcVecPos.Phi()   - mcVecNeg.Phi());
      registry.fill(HIST("hDebugDeltaDeltaPhi"),
                    TVector2::Phi_mpi_pi(deltaPhiReco - deltaPhiMC));
      registry.fill(HIST("hDebugDeltaDeltaPhiVsPtMC"), pairPtMC,
                    TVector2::Phi_mpi_pi(deltaPhiReco - deltaPhiMC));

      registry.fill(HIST("hDebugDeltaDeltaPhiVsPtMCPos"), PtMCPos,
                    TVector2::Phi_mpi_pi(deltaPhiReco - deltaPhiMC));
      registry.fill(HIST("hDebugDeltaDeltaPhiVsPtMCNeg"), PtMCNeg,
                    TVector2::Phi_mpi_pi(deltaPhiReco - deltaPhiMC));

      // [Debug-7] plain reco Delta-phi vs pair-level kinematics
      registry.fill(HIST("hDebugDeltaPhiVsDeltaPairPt"), pairPtReco - pairPtMC, deltaPhiReco);
      registry.fill(HIST("hDebugDeltaPhiVsPairPtReco"), pairPtReco, deltaPhiReco);
      registry.fill(HIST("hDebugDeltaPhiVsPairMassReco"), pairMassReco, deltaPhiReco);

      // [Debug-14] (1) Delta-theta = theta^{reco}-theta^{MC} vs pair pT^{reco}
      registry.fill(HIST("hDebugDeltaOpeningAngleVsPairPtReco"), pairPtReco, openingAngleReco - openingAngleMC);

      // [Debug-14] (2) Delta-pT_mu = pT_mu^{reco}-pT_mu^{MC} vs pair pT^{reco},
      // filled once per leg so either scenario A (one leg dominates) or
      // scenario B (both legs' Delta-pT stay small) is directly visible.
      registry.fill(HIST("hDebugDeltaPtMuVsPairPtReco"), pairPtReco, recoVec1.Pt() - mcVec1.Pt());
      registry.fill(HIST("hDebugDeltaPtMuVsPairPtReco"), pairPtReco, recoVec2.Pt() - mcVec2.Pt());
    }

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

    // [Debug-5] Pair mass / pT resolution vs muon pT asymmetry (PtMCPos/PtMCNeg
    // already identified above for the Delta(Delta-phi) check)
    const float deltaPtMuonsMC = PtMCPos - PtMCNeg;
    registry.fill(HIST("hPairMassResoVsDeltaPtMuonsMC"), deltaPtMuonsMC, pairMassReco - pairMassMC);
    if (pairPtMC > 0.f)
      registry.fill(HIST("hPairPtResoVsDeltaPtMuonsMC"), deltaPtMuonsMC, (pairPtReco - pairPtMC) / pairPtMC);

    // [Debug-6] Pair pT resolution vs softer-leg pT
    const float softMuonPtMC = (mcVec1.Pt() < mcVec2.Pt()) ? mcVec1.Pt() : mcVec2.Pt();
    if (pairPtMC > 0.f)
      registry.fill(HIST("hPairPtResoVsSoftMuonPtMC"), softMuonPtMC, (pairPtReco - pairPtMC) / pairPtMC);

    // [Debug-8] Opening angle theta / DeltaR vs mass resolution.
    // M^2 ~ p_1 p_2 theta^2 for small theta: if the mass deficit tracks the
    // *reco* angle, the angular separation is mismeasured; if it is flat in
    // angle, the bias is in the single-track momentum scale instead.
    // (openingAngleMC/openingAngleReco computed earlier, before the kinematic cuts.)
    const float trkDEtaMC = mcVec1.Eta() - mcVec2.Eta();
    const float trkDPhiMC = TVector2::Phi_mpi_pi(mcVec1.Phi() - mcVec2.Phi());
    const float deltaRMC = std::sqrt(trkDEtaMC * trkDEtaMC + trkDPhiMC * trkDPhiMC);
    const float trkDEtaReco = recoVec1.Eta() - recoVec2.Eta();
    const float trkDPhiReco = TVector2::Phi_mpi_pi(recoVec1.Phi() - recoVec2.Phi());
    const float deltaRReco = std::sqrt(trkDEtaReco * trkDEtaReco + trkDPhiReco * trkDPhiReco);

    registry.fill(HIST("hResponseMatrixOpeningAngle"), openingAngleMC, openingAngleReco);
    registry.fill(HIST("hResponseMatrixDeltaR"), deltaRMC, deltaRReco);
    registry.fill(HIST("hPairMassResoVsOpeningAngleMC"), openingAngleMC, pairMassReco - pairMassMC);
    registry.fill(HIST("hPairMassResoVsOpeningAngleReco"), openingAngleReco, pairMassReco - pairMassMC);
    registry.fill(HIST("hPairMassResoVsDeltaRMC"), deltaRMC, pairMassReco - pairMassMC);
    registry.fill(HIST("hPairMassResoVsDeltaRReco"), deltaRReco, pairMassReco - pairMassMC);

    // [Debug-9] Pair pT / mass resolution vs pair rapidity y^{MC}
    registry.fill(HIST("hPairMassResoVsRapMC"), pairRapMC, pairMassReco - pairMassMC);
    if (pairPtMC > 0.f)
      registry.fill(HIST("hPairPtResoVsRapMC"), pairRapMC, (pairPtReco - pairPtMC) / pairPtMC);

    // [Debug-15] Verification method B: mass resolution normalized by true pair pT
    if (pairPtMC > 0.f) {
      registry.fill(HIST("hPairMassResoVsPtMC"), pairPtMC, pairMassReco - pairMassMC);
      registry.fill(HIST("hPairMassResoNormVsPtMC"), pairPtMC, (pairMassReco - pairMassMC) / pairPtMC);
    }

    // [Debug-10] Pair pT / mass resolution vs helicity decay angle cos(theta_HE).
    // cos(theta_HE) is the angle of mu+ in the dimuon rest frame w.r.t. the
    // dimuon's own flight direction in the lab frame (standard helicity-frame
    // z-axis for an exclusively produced resonance).
    if (pairMC.Vect().Mag2() > 0.f && pairReco.Vect().Mag2() > 0.f) {
      const TVector3 zAxisMC = pairMC.Vect().Unit();
      TLorentzVector muPlusHEMC = mcVecPos;
      muPlusHEMC.Boost(-pairMC.BoostVector());

      const TVector3 zAxisReco = pairReco.Vect().Unit();
      TLorentzVector muPlusHEReco = recoVecPos;
      muPlusHEReco.Boost(-pairReco.BoostVector());

      if (muPlusHEMC.Vect().Mag2() > 0.f && muPlusHEReco.Vect().Mag2() > 0.f) {
        const float cosThetaHEMC = muPlusHEMC.Vect().Unit().Dot(zAxisMC);
        const float cosThetaHEReco = muPlusHEReco.Vect().Unit().Dot(zAxisReco);

        registry.fill(HIST("hResponseMatrixCosThetaHE"), cosThetaHEMC, cosThetaHEReco);
        registry.fill(HIST("hPairMassResoVsCosThetaHEMC"), cosThetaHEMC, pairMassReco - pairMassMC);
        if (pairPtMC > 0.f)
          registry.fill(HIST("hPairPtResoVsCosThetaHEMC"), cosThetaHEMC, (pairPtReco - pairPtMC) / pairPtMC);
      }
    }

    // [Debug-11] Pair mass / pT resolution vs min(N_{MFT clusters}) of the two legs
    const int minNMFTClusters = std::min(nMFTClusters1, nMFTClusters2);
    registry.fill(HIST("hPairMassResoVsMinNMFTClusters"), minNMFTClusters, pairMassReco - pairMassMC);
    if (pairPtMC > 0.f)
      registry.fill(HIST("hPairPtResoVsMinNMFTClusters"), minNMFTClusters, (pairPtReco - pairPtMC) / pairPtMC);

    // [Debug-12] Pair p_{T} / mass (raw reco distribution and resolution) vs
    // min(N_{MCH clusters}) of the two legs (UDFwdTracksExtra_001::fNClusters,
    // available directly on tr1/tr2, no extra table join needed).
    const int minNClusters = std::min(static_cast<int>(tr1.nClusters()), static_cast<int>(tr2.nClusters()));
    registry.fill(HIST("hPairPtRecoVsMinNClusters"), minNClusters, pairPtReco);
    registry.fill(HIST("hPairMassRecoVsMinNClusters"), minNClusters, pairMassReco);
    registry.fill(HIST("hPairMassResoVsMinNClusters"), minNClusters, pairMassReco - pairMassMC);
    if (pairPtMC > 0.f)
      registry.fill(HIST("hPairPtResoVsMinNClusters"), minNClusters, (pairPtReco - pairPtMC) / pairPtMC);

    // [Debug-13] Softer-leg (by total momentum) p vs pair mass. The soft leg
    // is picked by p^{MC} here, independent of the pT-based softness used in
    // [Debug-6]/[Debug-11]/[Debug-12].
    const bool track1IsSoftByP = mcVec1.P() < mcVec2.P();
    const float softLegPMC = track1IsSoftByP ? mcVec1.P() : mcVec2.P();
    const float softLegPReco = track1IsSoftByP ? recoVec1.P() : recoVec2.P();
    registry.fill(HIST("hPairMassRecoVsSoftLegPMC"), softLegPMC, pairMassReco);
    if (softLegPMC > 0.f)
      registry.fill(HIST("hPairMassResoVsSoftLegPReso"), (softLegPReco - softLegPMC) / softLegPMC, pairMassReco - pairMassMC);

    // [Debug-14] (3) Mass-deficit-tail gate: compare Delta-theta / Delta-pT_mu
    // for accepted pairs whose reco mass falls in the tail against the
    // inclusive (all accepted pairs) distributions.
    const float deltaOpeningAngle = openingAngleReco - openingAngleMC;
    registry.fill(HIST("hDebugDeltaOpeningAngle_Inclusive"), deltaOpeningAngle);
    registry.fill(HIST("hDebugDeltaPtMu_Inclusive"), recoVec1.Pt() - mcVec1.Pt());
    registry.fill(HIST("hDebugDeltaPtMu_Inclusive"), recoVec2.Pt() - mcVec2.Pt());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugDeltaOpeningAngle_MassTail"), deltaOpeningAngle);
      registry.fill(HIST("hDebugDeltaPtMu_MassTail"), recoVec1.Pt() - mcVec1.Pt());
      registry.fill(HIST("hDebugDeltaPtMu_MassTail"), recoVec2.Pt() - mcVec2.Pt());
    }

    // [Debug-14] (4) Single-muon kinematics + R_abs, same Inclusive/MassTail
    // gate as (3), to look for a distinct sub-population in single-track
    // phase space rather than in the reco-MC residuals.
    registry.fill(HIST("hDebugPtMu_Inclusive"), recoVec1.Pt());
    registry.fill(HIST("hDebugPtMu_Inclusive"), recoVec2.Pt());
    registry.fill(HIST("hDebugPtMuMC_Inclusive"), mcVec1.Pt());
    registry.fill(HIST("hDebugPtMuMC_Inclusive"), mcVec2.Pt());
    registry.fill(HIST("hDebugPtMuRecoVsMC_Inclusive"), mcVec1.Pt(), recoVec1.Pt());
    registry.fill(HIST("hDebugPtMuRecoVsMC_Inclusive"), mcVec2.Pt(), recoVec2.Pt());
    registry.fill(HIST("hDebugEtaMu_Inclusive"), recoVec1.Eta());
    registry.fill(HIST("hDebugEtaMu_Inclusive"), recoVec2.Eta());
    registry.fill(HIST("hDebugPhiMu_Inclusive"), recoVec1.Phi());
    registry.fill(HIST("hDebugPhiMu_Inclusive"), recoVec2.Phi());
    registry.fill(HIST("hDebugRAbsMu_Inclusive"), tr1.rAtAbsorberEnd());
    registry.fill(HIST("hDebugRAbsMu_Inclusive"), tr2.rAtAbsorberEnd());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugPtMu_MassTail"), recoVec1.Pt());
      registry.fill(HIST("hDebugPtMu_MassTail"), recoVec2.Pt());
      registry.fill(HIST("hDebugPtMuMC_MassTail"), mcVec1.Pt());
      registry.fill(HIST("hDebugPtMuMC_MassTail"), mcVec2.Pt());
      registry.fill(HIST("hDebugPtMuRecoVsMC_MassTail"), mcVec1.Pt(), recoVec1.Pt());
      registry.fill(HIST("hDebugPtMuRecoVsMC_MassTail"), mcVec2.Pt(), recoVec2.Pt());
      registry.fill(HIST("hDebugEtaMu_MassTail"), recoVec1.Eta());
      registry.fill(HIST("hDebugEtaMu_MassTail"), recoVec2.Eta());
      registry.fill(HIST("hDebugPhiMu_MassTail"), recoVec1.Phi());
      registry.fill(HIST("hDebugPhiMu_MassTail"), recoVec2.Phi());
      registry.fill(HIST("hDebugRAbsMu_MassTail"), tr1.rAtAbsorberEnd());
      registry.fill(HIST("hDebugRAbsMu_MassTail"), tr2.rAtAbsorberEnd());
    }

    // [Debug-14] (5) Track fit quality (chi2, chi2MatchMCHMFT), same
    // Inclusive/MassTail gate.
    registry.fill(HIST("hDebugChi2Mu_Inclusive"), tr1.chi2());
    registry.fill(HIST("hDebugChi2Mu_Inclusive"), tr2.chi2());
    registry.fill(HIST("hDebugChi2MatchMCHMFTMu_Inclusive"), tr1.chi2MatchMCHMFT());
    registry.fill(HIST("hDebugChi2MatchMCHMFTMu_Inclusive"), tr2.chi2MatchMCHMFT());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugChi2Mu_MassTail"), tr1.chi2());
      registry.fill(HIST("hDebugChi2Mu_MassTail"), tr2.chi2());
      registry.fill(HIST("hDebugChi2MatchMCHMFTMu_MassTail"), tr1.chi2MatchMCHMFT());
      registry.fill(HIST("hDebugChi2MatchMCHMFTMu_MassTail"), tr2.chi2MatchMCHMFT());
    }

    // [Debug-14] (6) N_{MCH clusters} / pDCA, same Inclusive/MassTail gate.
    registry.fill(HIST("hDebugNClustersMu_Inclusive"), tr1.nClusters());
    registry.fill(HIST("hDebugNClustersMu_Inclusive"), tr2.nClusters());
    registry.fill(HIST("hDebugPDcaMu_Inclusive"), tr1.pDca());
    registry.fill(HIST("hDebugPDcaMu_Inclusive"), tr2.pDca());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugNClustersMu_MassTail"), tr1.nClusters());
      registry.fill(HIST("hDebugNClustersMu_MassTail"), tr2.nClusters());
      registry.fill(HIST("hDebugPDcaMu_MassTail"), tr1.pDca());
      registry.fill(HIST("hDebugPDcaMu_MassTail"), tr2.pDca());
    }

    // [Debug-14] (7) N_{MFT clusters}, same Inclusive/MassTail gate.
    registry.fill(HIST("hDebugNMFTClustersMu_Inclusive"), nMFTClusters1);
    registry.fill(HIST("hDebugNMFTClustersMu_Inclusive"), nMFTClusters2);
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugNMFTClustersMu_MassTail"), nMFTClusters1);
      registry.fill(HIST("hDebugNMFTClustersMu_MassTail"), nMFTClusters2);
    }

    // [Debug-14] (8) MC-truth MCH-MFT match correctness (mcMask fake/noise
    // bits), same Inclusive/MassTail gate. bit 7 (0x80) = isFake: the MCH leg
    // and MFT leg of this global track do not trace back to the same true MC
    // particle. bit 6 (0x40) = isNoise: matched to combinatorial noise.
    static constexpr uint8_t kMcMaskFakeBit = 0x1 << 7;
    static constexpr uint8_t kMcMaskNoiseBit = 0x1 << 6;
    registry.fill(HIST("hDebugMcMaskFakeMu_Inclusive"), static_cast<float>((tr1.mcMask() & kMcMaskFakeBit) != 0));
    registry.fill(HIST("hDebugMcMaskFakeMu_Inclusive"), static_cast<float>((tr2.mcMask() & kMcMaskFakeBit) != 0));
    registry.fill(HIST("hDebugMcMaskNoiseMu_Inclusive"), static_cast<float>((tr1.mcMask() & kMcMaskNoiseBit) != 0));
    registry.fill(HIST("hDebugMcMaskNoiseMu_Inclusive"), static_cast<float>((tr2.mcMask() & kMcMaskNoiseBit) != 0));
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugMcMaskFakeMu_MassTail"), static_cast<float>((tr1.mcMask() & kMcMaskFakeBit) != 0));
      registry.fill(HIST("hDebugMcMaskFakeMu_MassTail"), static_cast<float>((tr2.mcMask() & kMcMaskFakeBit) != 0));
      registry.fill(HIST("hDebugMcMaskNoiseMu_MassTail"), static_cast<float>((tr1.mcMask() & kMcMaskNoiseBit) != 0));
      registry.fill(HIST("hDebugMcMaskNoiseMu_MassTail"), static_cast<float>((tr2.mcMask() & kMcMaskNoiseBit) != 0));
    }

    // [Debug-14] (9) Confound check: chi2 / chi2MatchMCHMFT vs pT, same
    // Inclusive/MassTail gate. Lets the MassTail chi2-pT correlation shape be
    // compared directly against the Inclusive one, to see whether (5)'s chi2
    // correlation is independent of the pT correlation already found, or
    // just rides on it (soft tracks generically fit worse).
    registry.fill(HIST("hDebugChi2VsPtMu_Inclusive"), recoVec1.Pt(), tr1.chi2());
    registry.fill(HIST("hDebugChi2VsPtMu_Inclusive"), recoVec2.Pt(), tr2.chi2());
    registry.fill(HIST("hDebugChi2MatchMCHMFTVsPtMu_Inclusive"), recoVec1.Pt(), tr1.chi2MatchMCHMFT());
    registry.fill(HIST("hDebugChi2MatchMCHMFTVsPtMu_Inclusive"), recoVec2.Pt(), tr2.chi2MatchMCHMFT());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugChi2VsPtMu_MassTail"), recoVec1.Pt(), tr1.chi2());
      registry.fill(HIST("hDebugChi2VsPtMu_MassTail"), recoVec2.Pt(), tr2.chi2());
      registry.fill(HIST("hDebugChi2MatchMCHMFTVsPtMu_MassTail"), recoVec1.Pt(), tr1.chi2MatchMCHMFT());
      registry.fill(HIST("hDebugChi2MatchMCHMFTVsPtMu_MassTail"), recoVec2.Pt(), tr2.chi2MatchMCHMFT());
    }

    // [Debug-16] Delta-pT_mu vs single-muon reco pT, same Inclusive/MassTail
    // gate: localises where in reco-pT space the MassTail's momentum-
    // underestimated sub-population (seen in hDebugDeltaPtMu_MassTail) sits.
    registry.fill(HIST("hDebugDeltaPtMuVsRecoPt_Inclusive"), recoVec1.Pt(), recoVec1.Pt() - mcVec1.Pt());
    registry.fill(HIST("hDebugDeltaPtMuVsRecoPt_Inclusive"), recoVec2.Pt(), recoVec2.Pt() - mcVec2.Pt());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugDeltaPtMuVsRecoPt_MassTail"), recoVec1.Pt(), recoVec1.Pt() - mcVec1.Pt());
      registry.fill(HIST("hDebugDeltaPtMuVsRecoPt_MassTail"), recoVec2.Pt(), recoVec2.Pt() - mcVec2.Pt());
    }

    // [Debug-17] (1) Delta-pT_mu vs single-muon phi (reco): detector-sector
    // dependence check, same Inclusive/MassTail gate.
    registry.fill(HIST("hDebugDeltaPtMuVsPhiMu_Inclusive"), recoVec1.Phi(), recoVec1.Pt() - mcVec1.Pt());
    registry.fill(HIST("hDebugDeltaPtMuVsPhiMu_Inclusive"), recoVec2.Phi(), recoVec2.Pt() - mcVec2.Pt());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugDeltaPtMuVsPhiMu_MassTail"), recoVec1.Phi(), recoVec1.Pt() - mcVec1.Pt());
      registry.fill(HIST("hDebugDeltaPtMuVsPhiMu_MassTail"), recoVec2.Phi(), recoVec2.Pt() - mcVec2.Pt());
    }

    // [Debug-17] (2) Delta-pT_mu split by charge: if the momentum-
    // underestimated sub-population is concentrated in mu+ or mu-, that
    // points back at a charge-dependent bias (e.g. an MFT/MCH rotation/
    // twist, as in [Debug-1]).
    if (tr1.sign() > 0)
      registry.fill(HIST("hDebugDeltaPtMu_Inclusive_MuPlus"), recoVec1.Pt() - mcVec1.Pt());
    else
      registry.fill(HIST("hDebugDeltaPtMu_Inclusive_MuMinus"), recoVec1.Pt() - mcVec1.Pt());
    if (tr2.sign() > 0)
      registry.fill(HIST("hDebugDeltaPtMu_Inclusive_MuPlus"), recoVec2.Pt() - mcVec2.Pt());
    else
      registry.fill(HIST("hDebugDeltaPtMu_Inclusive_MuMinus"), recoVec2.Pt() - mcVec2.Pt());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      if (tr1.sign() > 0)
        registry.fill(HIST("hDebugDeltaPtMu_MassTail_MuPlus"), recoVec1.Pt() - mcVec1.Pt());
      else
        registry.fill(HIST("hDebugDeltaPtMu_MassTail_MuMinus"), recoVec1.Pt() - mcVec1.Pt());
      if (tr2.sign() > 0)
        registry.fill(HIST("hDebugDeltaPtMu_MassTail_MuPlus"), recoVec2.Pt() - mcVec2.Pt());
      else
        registry.fill(HIST("hDebugDeltaPtMu_MassTail_MuMinus"), recoVec2.Pt() - mcVec2.Pt());
    }

    // [Debug-18] Approximate (x,y) position at the absorber end
    // (x ~ R_abs*cos(phi^{reco}), y ~ R_abs*sin(phi^{reco})): same
    // Inclusive/MassTail gate, plus a third map restricted to legs whose own
    // pT is anomalously underestimated (< deltaPtMuAnomalousMax), to
    // directly localise the momentum-underestimated sub-population.
    const float xApprox1 = tr1.rAtAbsorberEnd() * std::cos(recoVec1.Phi());
    const float yApprox1 = tr1.rAtAbsorberEnd() * std::sin(recoVec1.Phi());
    const float xApprox2 = tr2.rAtAbsorberEnd() * std::cos(recoVec2.Phi());
    const float yApprox2 = tr2.rAtAbsorberEnd() * std::sin(recoVec2.Phi());
    registry.fill(HIST("hDebugXYMu_Inclusive"), xApprox1, yApprox1);
    registry.fill(HIST("hDebugXYMu_Inclusive"), xApprox2, yApprox2);
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugXYMu_MassTail"), xApprox1, yApprox1);
      registry.fill(HIST("hDebugXYMu_MassTail"), xApprox2, yApprox2);
    }
    if (recoVec1.Pt() - mcVec1.Pt() < deltaPtMuAnomalousMax)
      registry.fill(HIST("hDebugXYMu_AnomalousPt"), xApprox1, yApprox1);
    if (recoVec2.Pt() - mcVec2.Pt() < deltaPtMuAnomalousMax)
      registry.fill(HIST("hDebugXYMu_AnomalousPt"), xApprox2, yApprox2);

    // [Debug-19] (1) Delta-pT_mu vs pDca^{reco}, same Inclusive/MassTail gate
    registry.fill(HIST("hDebugDeltaPtMuVsPDca_Inclusive"), tr1.pDca(), recoVec1.Pt() - mcVec1.Pt());
    registry.fill(HIST("hDebugDeltaPtMuVsPDca_Inclusive"), tr2.pDca(), recoVec2.Pt() - mcVec2.Pt());
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugDeltaPtMuVsPDca_MassTail"), tr1.pDca(), recoVec1.Pt() - mcVec1.Pt());
      registry.fill(HIST("hDebugDeltaPtMuVsPDca_MassTail"), tr2.pDca(), recoVec2.Pt() - mcVec2.Pt());
    }

    // [Debug-19] (2) Event vertex position: Inclusive/MassTail (one fill per
    // pair) plus AnomalousPt (one fill per anomalous leg found in the pair).
    registry.fill(HIST("hDebugVertexXY_Inclusive"), vtxX, vtxY);
    registry.fill(HIST("hDebugVertexZ_Inclusive"), vtxZ);
    if (pairMassReco > massTailMin && pairMassReco < massTailMax) {
      registry.fill(HIST("hDebugVertexXY_MassTail"), vtxX, vtxY);
      registry.fill(HIST("hDebugVertexZ_MassTail"), vtxZ);
    }
    if (recoVec1.Pt() - mcVec1.Pt() < deltaPtMuAnomalousMax || recoVec2.Pt() - mcVec2.Pt() < deltaPtMuAnomalousMax) {
      registry.fill(HIST("hDebugVertexXY_AnomalousPt"), vtxX, vtxY);
      registry.fill(HIST("hDebugVertexZ_AnomalousPt"), vtxZ);
    }

    // [Debug-20] Muon pT asymmetry (reco) vs pair pT (reco): tests whether
    // an asymmetric pair (one leg much softer/harder than the other) shows
    // up as a larger reco pair pT, as expected if the true pair pT is small
    // (exclusive UPC production).
    const float deltaPtMuonsReco = std::abs(recoVec1.Pt() - recoVec2.Pt());
    registry.fill(HIST("hDeltaPtMuonsRecoVsPairPtReco_Inclusive"), pairPtReco, deltaPtMuonsReco);
    if (pairMassReco > massTailMin && pairMassReco < massTailMax)
      registry.fill(HIST("hDeltaPtMuonsRecoVsPairPtReco_MassTail"), pairPtReco, deltaPtMuonsReco);

    // [Debug-21] Muon pair track-time difference, same Inclusive/MassTail gate
    const float deltaTrackTime = std::abs(tr1.trackTime() - tr2.trackTime());
    registry.fill(HIST("hDebugDeltaTrackTime_Inclusive"), deltaTrackTime);
    if (pairMassReco > massTailMin && pairMassReco < massTailMax)
      registry.fill(HIST("hDebugDeltaTrackTime_MassTail"), deltaTrackTime);

    // [Debug-22] Pair mass vs single-muon pT/eta (continuous, both legs)
    registry.fill(HIST("hDebugPairMassRecoVsPtMu"), recoVec1.Pt(), pairMassReco);
    registry.fill(HIST("hDebugPairMassRecoVsPtMu"), recoVec2.Pt(), pairMassReco);
    registry.fill(HIST("hDebugPairMassRecoVsEtaMu"), recoVec1.Eta(), pairMassReco);
    registry.fill(HIST("hDebugPairMassRecoVsEtaMu"), recoVec2.Eta(), pairMassReco);

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

      if (v.Eta() <= etaMin || v.Eta() >= etaMax)
        continue;
      registry.fill(HIST("hCutFlowMC"), 9);

      if (v.Pt() < pTmin)
        continue;
      registry.fill(HIST("hCutFlowMC"), 10);

      registry.fill(HIST("hEffTrkPhiMC"), v.Phi());
      registry.fill(HIST("hEffTrkEtaMC"), v.Eta());
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

    (void)mcCollisions;
  }

  // ===========================================================================
  // Process: MC reco
  // ===========================================================================
  void processMcReco(CandidatesFwd const& eventCandidates,
                     CompleteFwdTracks const& fwdTracks,
                     o2::aod::UDFwdTracksCls const& fwdTrackClusters,
                     o2::aod::UDMcCollisions const&,
                     o2::aod::UDMcParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);

    // [Debug-11] N MFT clusters per track, from the per-cluster (x,y,z) table
    // (only populated when the UD table producer ran with fSaveMFTClusters).
    std::unordered_map<int32_t, int> nMFTClustersPerTrack;
    for (const auto& cls : fwdTrackClusters) {
      ++nMFTClustersPerTrack[cls.udFwdTrackId()];
    }

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
        registry.fill(HIST("hEffTrkPtReco"), vReco.Pt());

        // [Debug-11] N MFT clusters used in this track's fit
        auto nClsIt = nMFTClustersPerTrack.find(static_cast<int32_t>(tr.globalIndex()));
        const int nCls = (nClsIt != nMFTClustersPerTrack.end()) ? nClsIt->second : 0;

        // Resolution: single-track histograms
        fillSingleTrackResolution(tr, mc, nCls);

        // [Debug-1] (q/pT) residual per charge
        fillDebugQOverPtResidual(tr, mc, nCls);
      }

      registry.fill(HIST("hMuonMultReco"), static_cast<float>(goodTrkIds.size()));

      // --- 2-track requirement for cutflow step 2 ---
      if (candidateTrkIds.size() == static_cast<size_t>(reqMatchMFT))
        registry.fill(HIST("hCutFlowReco"), 2); // 2: Pass Exact 2 Tracks (candidate level)

      // --- Pair loop: only when exactly reqMatchMFT good muon tracks found ---
      if (goodTrkIds.size() != static_cast<size_t>(reqMatchMFT))
        continue;
      registry.fill(HIST("hCutFlowReco"), 3); // 3: Pass Exact MatchMFT

      // [Debug-19] Event (collision) vertex position, same for every pair in
      // this candidate.
      auto cand = eventCandidates.iteratorAt(candId);
      const float vtxX = cand.posX();
      const float vtxY = cand.posY();
      const float vtxZ = cand.posZ();

      for (size_t i = 0; i < goodTrkIds.size(); ++i) {
        auto tr1 = fwdTracks.iteratorAt(goodTrkIds[i]);
        const auto& mc1 = tr1.udMcParticle();
        auto nCls1It = nMFTClustersPerTrack.find(static_cast<int32_t>(tr1.globalIndex()));
        const int nCls1 = (nCls1It != nMFTClustersPerTrack.end()) ? nCls1It->second : 0;

        for (size_t j = i + 1; j < goodTrkIds.size(); ++j) {
          auto tr2 = fwdTracks.iteratorAt(goodTrkIds[j]);
          const auto& mc2 = tr2.udMcParticle();
          auto nCls2It = nMFTClustersPerTrack.find(static_cast<int32_t>(tr2.globalIndex()));
          const int nCls2 = (nCls2It != nMFTClustersPerTrack.end()) ? nCls2It->second : 0;

          fillPairAnalysis(tr1, mc1, tr2, mc2, nCls1, nCls2, vtxX, vtxY, vtxZ);
        }
      }
    }
  }

  PROCESS_SWITCH(UPCMuonAnalysisAlignmentDebug, processMCTrue,
                 "Fill MC-truth single-track and dimuon efficiency histograms", true);
  PROCESS_SWITCH(UPCMuonAnalysisAlignmentDebug, processMcReco,
                 "Fill reco single-track and dimuon resolution + efficiency histograms", true);
};

// ============================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCMuonAnalysisAlignmentDebug>(cfgc, TaskName{"upc-muon-alignment-debug"}),
  };
}
