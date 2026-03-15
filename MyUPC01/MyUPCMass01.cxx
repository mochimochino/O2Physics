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

struct MyUPCMass01Task {

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histo"};
  // Configurable<int> nBinsMass{"nBinsMass", 200, "N bins in InvMass histo (50 MeV/c^2 per bin)"};
  // (10 MeV/c^2)
  Configurable<int> nBinsMass{"nBinsMass", 1000, "N bins in InvMass histo (10 MeV/c^2 per bin)"};

  Configurable<float> cutZNAEnergy{"cutZNAEnergy", 1.0f, "ZNA energy threshold [TeV]"};
  Configurable<float> cutZNCEnergy{"cutZNCEnergy", 1.0f, "ZNC energy threshold [TeV]"};
  // targetTopology: 0 = 0n0n, 1 = Xn0n, 2 = 0nXn, 3 = XnXn, -1 = No Cut
  Configurable<int> targetTopology{"targetTopology", -1, "Required neutron topology (3 for XnXn)"};

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
  }

  //____________________________________________________________________________________________

  template <typename TTrack1, typename TTrack2>
  void processCandidate(UDCollisionsFwd::iterator const& collision, TTrack1& tr1, TTrack2& tr2)
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
    }

    if (isLikeSign) {
      registry.fill(HIST("hSelectionCounter"), 12);
      registry.fill(HIST("MMuonLike"), p.M());
      registry.fill(HIST("PtMuonLike"), p.Pt());
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

    // Iterate through candidates. Avoid crashes when not exactly 2 tracks.
    for (const auto& item : tracksPerCand) {
      if (item.second.size() < 2) {
        continue; // Require at least 2 tracks to make a pair
      }

      int32_t candID = item.first;


      if (candID >= zdcs.size()) {
         continue; // データが存在しないインデックスの場合はスキップ
      }
      // ZDC topology selection
      const auto& zdc = zdcs.iteratorAt(candID);
      
      float eZNA = zdc.energyCommonZNA();
      float eZNC = zdc.energyCommonZNC();
      bool hasZNA = eZNA > static_cast<float>(cutZNAEnergy);
      bool hasZNC = eZNC > static_cast<float>(cutZNCEnergy);

      int currentTopology = 0;
      if      (!hasZNA && !hasZNC) currentTopology = 0; // 0n0n
      else if ( hasZNA && !hasZNC) currentTopology = 1; // Xn0n
      else if (!hasZNA &&  hasZNC) currentTopology = 2; // 0nXn
      else                         currentTopology = 3; // XnXn

      if (targetTopology > -1 && currentTopology != targetTopology) {
        continue;
      }
      // --------------------------------------------------

      const auto& collision = eventCandidates.iteratorAt(candID);
      
      // Iterate over pairs if there are >= 2 tracks (or just first 2 like UDTutorial_06)
      for (size_t i = 0; i < item.second.size() - 1; ++i) {
        for (size_t j = i + 1; j < item.second.size(); ++j) {
          int32_t trId1 = item.second[i];
          int32_t trId2 = item.second[j];
          const auto& tr1 = fwdTracks.iteratorAt(trId1);
          const auto& tr2 = fwdTracks.iteratorAt(trId2);
          processCandidate(collision, tr1, tr2);
        }
      }
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyUPCMass01Task>(cfgc, TaskName{"my-upc-mass-01"})};
}
