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
/// \brief Custom task to compute invariant mass of dimuon at forward rapidity for UPC events.
/// \author Takuma
/// \date 2026

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// O2Physics headers
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"

// ROOT headers
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct MyUPCTask {

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 200, "N bins in InvMass histo (50 MeV/c^2 per bin)"};
  Configurable<float> massMin{"massMin", 2.8f, "Min mass for single muon pT region"};
  Configurable<float> massMax{"massMax", 3.4f, "Max mass for single muon pT region"};
  
  // ZDC variables
  Configurable<float> cutZNAEnergy{"cutZNAEnergy", 1.0f, "Minimum ZNA energy (ADC) to be considered active (neutron)"};
  Configurable<float> cutZNCEnergy{"cutZNCEnergy", 1.0f, "Minimum ZNC energy (ADC) to be considered active (neutron)"};

  using UDCollisionsFwd = o2::aod::UDCollisions; // Option 2: Removed UDCollisionsSels and SelsFwd to run without helper task
  using fwdtraks = soa::Join<aod::UDFwdTracks, aod::UDFwdTracksExtra>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{80, -5.0, -2.0, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0.0, 10.0, "p_{T}"};
    const AxisSpec axisM{nBinsMass, 0.0, 10.0, "M_{mumu}"};
    const AxisSpec axisPhi{120, -TMath::Pi(), -TMath::Pi(), "#phi"};
    const AxisSpec axisRapidity{80, -5.0, -2.0, "#it{y}"};

    // create histograms
    // Fill counter to see effect of each selection criteria
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{15, 0., 15.}});
    TString SelectionCuts[13] = {"NoSelection", "V0A", "rAbs1", "rAbs2", "pDCA", "trackmatch", "eta1", "eta2", "pair pt", "pair rapidity", "mass_cut", "unlikesign", "likesign"};

    for (int i = 0; i < 13; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("hnumContrib", "hnumContrib;;#counts", kTH1D, {{100, -0.5, 99.5}});
    registry.add("hTracks1", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracks2", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksMuons", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});

    registry.add("etaMuon1", "etaMuon1", kTH1F, {axisEta});
    registry.add("ptMuon1", "ptMuon1", kTH1F, {axisPt});
    registry.add("phiMuon1", "phiMuon1", kTH1F, {axisPhi});

    registry.add("etaMuon2", "etaMuon2", kTH1F, {axisEta});
    registry.add("ptMuon2", "ptMuon2", kTH1F, {axisPt});
    registry.add("phiMuon2", "phiMuon2", kTH1F, {axisPhi});

    registry.add("YJpsi", "YJpsi", kTH1F, {axisRapidity});
    registry.add("PhiJpsi", "PhiJpsi", kTH1F, {axisPhi});

    registry.add("MMuonUnlike", "MMuonUnlike", kTH1F, {axisM});
    registry.add("PtMuonUnlike", "PtMuonUnlike", kTH1F, {axisPt});
    registry.add("MMuonLike", "MMuonLike", kTH1F, {axisM});
    registry.add("PtMuonLike", "PtMuonLike", kTH1F, {axisPt});

    registry.add("PtEachMuonUnlikeMassRegion", "Pt of each Muon (Unlike) in Mass Region", kTH1F, {axisPt});
    registry.add("PtEachMuonLikeMassRegion", "Pt of each Muon (Like) in Mass Region", kTH1F, {axisPt});

    // 2D Histograms for Mass vs Pt  
    registry.add("MassVsPtUnlike", "Mass vs Pt (Unlike);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)", kTH2D, {axisM, axisPt});
    registry.add("MassVsPtLike", "Mass vs Pt (Like);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)", kTH2D, {axisM, axisPt});

    // ZDC Histograms
    const AxisSpec axisZDC{1000, -2.5, 199.5, "ZDC Energy"};
    const AxisSpec axisTopology{4, -0.5, 3.5, "0n0n, Xn0n, 0nXn, XnXn"};
    registry.add("EnergyZNA", "Energy ZNA", kTH1F, {axisZDC});
    registry.add("EnergyZNC", "Energy ZNC", kTH1F, {axisZDC});
    registry.add("EnergyZNAvsZNC", "Energy ZNA vs ZNC; ZDC-C; ZDC-A", kTH2D, {axisZDC, axisZDC});
    auto hTopologyCounter = registry.add<TH1>("hTopologyCounter", "Neutron Topology;;Events", HistType::kTH1I, {{4, 0., 4.}});
    TString TopologyClasses[4] = {"0n0n", "Xn0n", "0nXn", "XnXn"};
    for (int i = 0; i < 4; i++) {
        hTopologyCounter->GetXaxis()->SetBinLabel(i + 1, TopologyClasses[i].Data());
    }

    // Topology specific Mass vs Pt (Unlike)
    registry.add("MassVsPtUnlike_0n0n", "Mass vs Pt (Unlike, 0n0n);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)", kTH2D, {axisM, axisPt});
    registry.add("MassVsPtUnlike_Xn0n", "Mass vs Pt (Unlike, Xn0n);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)", kTH2D, {axisM, axisPt});
    registry.add("MassVsPtUnlike_0nXn", "Mass vs Pt (Unlike, 0nXn);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)", kTH2D, {axisM, axisPt});
    registry.add("MassVsPtUnlike_XnXn", "Mass vs Pt (Unlike, XnXn);M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c)", kTH2D, {axisM, axisPt});

  }

  //____________________________________________________________________________________________

  template <typename TTrack1, typename TTrack2>
  void processCandidate(UDCollisionsFwd::iterator const& collision, TTrack1& tr1, TTrack2& tr2, float E_ZNA, float E_ZNC)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    registry.fill(HIST("hnumContrib"), collision.numContrib());
    registry.fill(HIST("hTracks1"), tr1.size());
    registry.fill(HIST("hTracks2"), tr2.size());
    registry.fill(HIST("hSelectionCounter"), 0);

    // V0A selection (Commented out for Option 2 since it requires UDCollisionsSelsFwd)
    /*
    const auto& ampsV0A = collision.amplitudesV0A();
    const auto& ampsRelBCsV0A = collision.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1) {
        if (ampsV0A[i] > 100.)
          return;
      }
    }
    */

    // T0A selection (Commented out for Option 2 since it requires UDCollisionsSelsFwd)
    /*
    const auto& ampsT0A = collision.amplitudesT0A();
    const auto& ampsRelBCsT0A = collision.ampRelBCsT0A();
    for (unsigned int i = 0; i < ampsT0A.size(); ++i) {
      if (std::abs(ampsRelBCsT0A[i]) <= 1) {
        if (ampsT0A[i] > 100.)
          return;
      }
    }
    */

    registry.fill(HIST("hSelectionCounter"), 1);

    // Get ZDC energies safely
    registry.fill(HIST("EnergyZNA"), E_ZNA);
    registry.fill(HIST("EnergyZNC"), E_ZNC);
    registry.fill(HIST("EnergyZNAvsZNC"), E_ZNC, E_ZNA);

    // Determine Topology
    bool isZNA = (E_ZNA > cutZNAEnergy);
    bool isZNC = (E_ZNC > cutZNCEnergy);
    int topology = 0; // 0: 0n0n
    if (isZNA && !isZNC) topology = 1; // 1: Xn0n (A-side only)
    else if (!isZNA && isZNC) topology = 2; // 2: 0nXn (C-side only)
    else if (isZNA && isZNC) topology = 3; // 3: XnXn (both)
    registry.fill(HIST("hTopologyCounter"), topology);

    // absorber end selection
    if (tr1.rAtAbsorberEnd() < 17.6 || tr1.rAtAbsorberEnd() > 89.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 2);
    if (tr2.rAtAbsorberEnd() < 17.6 || tr2.rAtAbsorberEnd() > 89.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 3);

    // pDCA cut
    auto checkPDCA = [](auto& tr) {
      double r = tr.rAtAbsorberEnd();
      double pdca = tr.pDca();
      if (r < 26.5) return pdca < 350.0;
      else return pdca < 200.0;
    };
    if (!checkPDCA(tr1) || !checkPDCA(tr2))
      return;
    registry.fill(HIST("hSelectionCounter"), 4);

    // MCH-MID match selection
    if ((tr1.chi2MatchMCHMID() < 0) || (tr2.chi2MatchMCHMID() < 0))
      return;
    registry.fill(HIST("hSelectionCounter"), 5);

    bool isUnlikeSign = (tr1.sign() + tr2.sign()) == 0;
    bool isLikeSign = ((tr1.sign() + tr2.sign()) != 0);

    // track selection
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassMuon);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassMuon);
    TLorentzVector p = p1 + p2;

    // eta cut on each track
    if (p1.Eta() <= -4.0 || p1.Eta() >= -2.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 6);
    if (p2.Eta() <= -4.0 || p2.Eta() >= -2.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 7);

    // pair pt cut (pT < 0.25 GeV/c)
    if (p.Pt() >= 0.25)
      return;
    registry.fill(HIST("hSelectionCounter"), 8);

    // pair rapidity cut (-4.0 < y < -2.5)
    if ((p.Rapidity() <= -4.0) || (p.Rapidity() >= -2.5))
      return;
    registry.fill(HIST("hSelectionCounter"), 9);

    // cuts on pair kinematics (modify this range depending on J/psi or generic mu-mu)
    // For general mu-mu, we might want to relax this cut or keep it wide
    if (!(p.M() > 1.0 && p.M() < 10.0))
      return;
    registry.fill(HIST("hSelectionCounter"), 10);

    registry.fill(HIST("hTracksMuons"), collision.numContrib());
    registry.fill(HIST("ptMuon1"), p1.Pt());
    registry.fill(HIST("ptMuon2"), p2.Pt());
    registry.fill(HIST("etaMuon1"), p1.Eta());
    registry.fill(HIST("etaMuon2"), p2.Eta());
    registry.fill(HIST("phiMuon1"), p1.Phi());
    registry.fill(HIST("phiMuon2"), p2.Phi());
    registry.fill(HIST("YJpsi"), p.Rapidity());
    registry.fill(HIST("PhiJpsi"), p.Phi());

    if (isUnlikeSign) {
      registry.fill(HIST("hSelectionCounter"), 11);
      registry.fill(HIST("MMuonUnlike"), p.M());
      registry.fill(HIST("PtMuonUnlike"), p.Pt());
      registry.fill(HIST("MassVsPtUnlike"), p.M(), p.Pt());
      
      // Fill Topology specific 2D histograms
      if (topology == 0) registry.fill(HIST("MassVsPtUnlike_0n0n"), p.M(), p.Pt());
      else if (topology == 1) registry.fill(HIST("MassVsPtUnlike_Xn0n"), p.M(), p.Pt());
      else if (topology == 2) registry.fill(HIST("MassVsPtUnlike_0nXn"), p.M(), p.Pt());
      else if (topology == 3) registry.fill(HIST("MassVsPtUnlike_XnXn"), p.M(), p.Pt());

      if (p.M() >= massMin && p.M() <= massMax) {
        // Fill single muon pT for both tracks if the pair mass is within the region
        registry.fill(HIST("PtEachMuonUnlikeMassRegion"), p1.Pt());
        registry.fill(HIST("PtEachMuonUnlikeMassRegion"), p2.Pt());
      }
    }

    if (isLikeSign) {
      registry.fill(HIST("hSelectionCounter"), 12);
      registry.fill(HIST("MMuonLike"), p.M());
      registry.fill(HIST("PtMuonLike"), p.Pt());
      registry.fill(HIST("MassVsPtLike"), p.M(), p.Pt());
      if (p.M() >= massMin && p.M() <= massMax) {
        registry.fill(HIST("PtEachMuonLikeMassRegion"), p1.Pt());
        registry.fill(HIST("PtEachMuonLikeMassRegion"), p2.Pt());
      }
    }
  }

  //____________________________________________________________________________________________

  // Template that collects all collision IDs and track per collision
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  //____________________________________________________________________________________________

  // process candidates with forward tracks
  void process(UDCollisionsFwd const& eventCandidates,
               fwdtraks const& fwdTracks,
               aod::UDZdcsReduced const& zdcs)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // Build a map for fast lookup of ZDC entries by collision ID
    std::unordered_map<int32_t, int32_t> collisionToZdcMap;
    for (const auto& zdc : zdcs) {
        collisionToZdcMap[zdc.udCollisionId()] = zdc.index();
    }

    // Iterate through candidates. Avoid crashes when not exactly 2 tracks.
    for (const auto& item : tracksPerCand) {
      if (item.second.size() < 2) {
        continue; // Require at least 2 tracks to make a pair
      }

      int32_t candID = item.first;
      const auto& collision = eventCandidates.iteratorAt(candID);
      
      // Look up ZDC energies for this collision
      float E_ZNA = 0.0f;
      float E_ZNC = 0.0f;
      auto itZdc = collisionToZdcMap.find(candID);
      if (itZdc != collisionToZdcMap.end()) {
          const auto& zdc = zdcs.iteratorAt(itZdc->second);
          E_ZNA = zdc.energyCommonZNA();
          E_ZNC = zdc.energyCommonZNC();
      }
      
      // Iterate over pairs if there are >= 2 tracks (or just first 2 like UDTutorial_06)
      for (size_t i = 0; i < item.second.size() - 1; ++i) {
        for (size_t j = i + 1; j < item.second.size(); ++j) {
          int32_t trId1 = item.second[i];
          int32_t trId2 = item.second[j];
          const auto& tr1 = fwdTracks.iteratorAt(trId1);
          const auto& tr2 = fwdTracks.iteratorAt(trId2);
          processCandidate(collision, tr1, tr2, E_ZNA, E_ZNC);
        }
      }
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyUPCTask>(cfgc, TaskName{"my-upc-01"})};
}
