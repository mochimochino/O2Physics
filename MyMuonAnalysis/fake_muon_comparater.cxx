#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <TMath.h>

using namespace o2;
using namespace o2::framework;

struct FakeMuonComparator {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  AxisSpec axisPt{100, 0.0, 10.0, "p_{T} [GeV/c]"};
  AxisSpec axisEta{50, -4.0, -2.5, "#eta"};
  AxisSpec axisPhi{64, -TMath::Pi(), TMath::Pi(), "#phi [rad]"};

  void init(InitContext const&)
  {
    // matchedQualityCuts histograms
    histos.add("hPt_matchedQualityCuts", "pT (Matched Quality Cuts)", kTH1F, {axisPt});
    histos.add("hEta_matchedQualityCuts", "#eta (Matched Quality Cuts)", kTH1F, {axisEta});
    histos.add("hPhi_matchedQualityCuts", "#phi (Matched Quality Cuts)", kTH1F, {axisPhi});
    histos.add("hEtaPt_matchedQualityCuts", "#eta vs pT (Matched Quality Cuts)", kTH2F, {axisPt, axisEta});

    // All Global Muon (trackType == 0) ---
    histos.add("hPt_GlobalMuon_ALL", "pT (Matched Quality Cuts)", kTH1F, {axisPt});
    histos.add("hEta_GlobalMuon_ALL", "#eta (Matched Quality Cuts)", kTH1F, {axisEta});
    histos.add("hPhi_GlobalMuon_ALL", "#phi (Matched Quality Cuts)", kTH1F, {axisPhi});
    histos.add("hEtaPt_GlobalMuon_ALL", "#eta vs pT (Matched Quality Cuts)", kTH2F, {axisPt, axisEta});

    // True Global Muon (trackType == 0, isFake == false) ---
    histos.add("hPt_GlobalMuon_True", "pT (True Global Muon)", kTH1F, {axisPt});
    histos.add("hEta_GlobalMuon_True", "#eta (True Global Muon)", kTH1F, {axisEta});
    histos.add("hPhi_GlobalMuon_True", "#phi (True Global Muon)", kTH1F, {axisPhi});
    histos.add("hEtaPt_GlobalMuon_True", "#eta vs pT (True Global Muon)", kTH2F, {axisPt, axisEta});

    // Fake Global Muon (trackType == 0, isFake == true) ---
    histos.add("hPt_GlobalMuon_Fake", "pT (Fake Global Muon)", kTH1F, {axisPt});
    histos.add("hEta_GlobalMuon_Fake", "#eta (Fake Global Muon)", kTH1F, {axisEta});
    histos.add("hPhi_GlobalMuon_Fake", "#phi (Fake Global Muon)", kTH1F, {axisPhi});
    histos.add("hEtaPt_GlobalMuon_Fake", "#eta vs pT (Fake Global Muon)", kTH2F, {axisPt, axisEta});
  }

  void process(soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks)
  {
    // document
    const int fakeBit = 0x1 << 7; 

    for (auto& trk : tracks) {
    // matchedQualityCuts
      if (trk.eta() < -4.0 || trk.eta() > -2.5) continue; 
      const float rAbs = trk.rAtAbsorberEnd();
      if (rAbs < 17.6 || rAbs > 89.5) continue; 
      const float pDca = trk.pDca();
      if (pDca < 0.0) continue; 
      if (rAbs < 26.5) {
          if (pDca > 594.0) continue;
      } else {
          if (pDca > 324.0) continue;
      }
      if (trk.chi2() < 0.0 || trk.chi2() > 1e6) continue;
      if (trk.chi2MatchMCHMID() < 0.0 || trk.chi2MatchMCHMID() > 1e6) continue; 
      if (trk.chi2MatchMCHMFT() < 0.0 || trk.chi2MatchMCHMFT() > 1e6) continue;

      // for chacking the cut (Right: matchedQualityCuts = GlobalMuon_ALL)
      histos.fill(HIST("hPt_matchedQualityCuts"), trk.pt());
      histos.fill(HIST("hEta_matchedQualityCuts"), trk.eta());
      histos.fill(HIST("hPhi_matchedQualityCuts"), trk.phi());
      histos.fill(HIST("hEtaPt_matchedQualityCuts"), trk.pt(), trk.eta());  

   
      // Track0: Global Muon
      histos.fill(HIST("hPt_GlobalMuon_ALL"), trk.pt());
      histos.fill(HIST("hEta_GlobalMuon_ALL"), trk.eta());
      histos.fill(HIST("hPhi_GlobalMuon_ALL"), trk.phi());
      histos.fill(HIST("hEtaPt_GlobalMuon_ALL"), trk.pt(), trk.eta());

      bool isFake = (trk.mcMask() & fakeBit) != 0;

      const float pt = trk.pt();
      const float eta = trk.eta();
      const float phi = trk.phi();

      if (isFake) {
        histos.fill(HIST("hPt_GlobalMuon_Fake"), pt);
        histos.fill(HIST("hEta_GlobalMuon_Fake"), eta);
        histos.fill(HIST("hPhi_GlobalMuon_Fake"), phi);
        histos.fill(HIST("hEtaPt_GlobalMuon_Fake"), pt, eta);
      } else {
        histos.fill(HIST("hPt_GlobalMuon_True"), pt);
        histos.fill(HIST("hEta_GlobalMuon_True"), eta);
        histos.fill(HIST("hPhi_GlobalMuon_True"), phi);
        histos.fill(HIST("hEtaPt_GlobalMuon_True"), pt, eta);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<FakeMuonComparator>(cfg)
  };
}