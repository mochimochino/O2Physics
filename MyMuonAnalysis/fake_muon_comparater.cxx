#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <cmath>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;

struct FakeMuonComparator {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> maxChi2Norm{"maxChi2Norm", 3.5, "Max Chi2 Norm"};
  Configurable<float> maxMatchChi2{"maxMatchChi2", 150.0, "Max Match Chi2"};
  Configurable<float> maxPDCA{"maxPDCA", 500.0, "Max p*DCA"};

  AxisSpec axisPt{100, 0.0, 4.0, "p_{T} [GeV/c]"};
  AxisSpec axisP{100, 0.0, 100.0, "p [GeV/c]"};
  AxisSpec axisEta{50, -4.0, -2.5, "#eta"};
  AxisSpec axisPhi{64, -TMath::Pi(), TMath::Pi(), "#phi [rad]"};

  // Quality check
  AxisSpec axisChi2Norm{100, 0.0, 10.0, "#chi^{2}_{MCH} / NDF"};
  AxisSpec axisMatchChi2{100, 0.0, 100.0, "#chi^{2}_{Match} (MCH-MFT)"};
  AxisSpec axisPDCA{100, 0.0, 500.0, "p #times DCA [cm#cdotGeV/c]"};

  void init(InitContext const&)
  {
    // ALL
    histos.add("hPt_GlobalMuon_ALL", "p_{T} (ALL);p_{T} [GeV/c];Counts", kTH1F, {axisPt});
    histos.add("hP_GlobalMuon_ALL", "p (ALL);p [GeV/c];Counts", kTH1F, {axisP});
    histos.add("hEta_GlobalMuon_ALL", "#eta (ALL);#eta;Counts", kTH1F, {axisEta});
    histos.add("hPhi_GlobalMuon_ALL", "#phi (ALL);#phi;Counts", kTH1F, {axisPhi});
    histos.add("hChi2Norm_GlobalMuon_ALL", "Normalized MCH #chi^{2} (ALL);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("hMatchChi2_GlobalMuon_ALL", "MCH-MFT Match #chi^{2} (ALL);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("hPDCA_GlobalMuon_ALL", "p #times DCA (ALL);p #times DCA;Counts", kTH1F, {axisPDCA});

    // True
    histos.add("hPt_GlobalMuon_True", "p_{T} (True);p_{T} [GeV/c];Counts", kTH1F, {axisPt});
    histos.add("hP_GlobalMuon_True", "p (True);p [GeV/c];Counts", kTH1F, {axisP});
    histos.add("hEta_GlobalMuon_True", "#eta (True);#eta;Counts", kTH1F, {axisEta});
    histos.add("hPhi_GlobalMuon_True", "#phi (True);#phi;Counts", kTH1F, {axisPhi});
    histos.add("hChi2Norm_GlobalMuon_True", "Normalized MCH #chi^{2} (True);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("hMatchChi2_GlobalMuon_True", "MCH-MFT Match #chi^{2} (True);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("hPDCA_GlobalMuon_True", "p #times DCA (True);p #times DCA;Counts", kTH1F, {axisPDCA});

    // Fake
    histos.add("hPt_GlobalMuon_Fake", "p_{T} (Fake);p_{T} [GeV/c];Counts", kTH1F, {axisPt});
    histos.add("hP_GlobalMuon_Fake", "p (Fake);p [GeV/c];Counts", kTH1F, {axisP});
    histos.add("hEta_GlobalMuon_Fake", "#eta (Fake);#eta;Counts", kTH1F, {axisEta});
    histos.add("hPhi_GlobalMuon_Fake", "#phi (Fake);#phi;Counts", kTH1F, {axisPhi});
    histos.add("hChi2Norm_GlobalMuon_Fake", "Normalized MCH #chi^{2} (Fake);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("hMatchChi2_GlobalMuon_Fake", "MCH-MFT Match #chi^{2} (Fake);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("hPDCA_GlobalMuon_Fake", "p #times DCA (Fake);p #times DCA;Counts", kTH1F, {axisPDCA});

    // 2D Check
    histos.add("hEtaPt_GlobalMuon_True", "#eta vs pT (True);p_{T};#eta", kTH2F, {axisPt, axisEta});
    histos.add("hEtaPt_GlobalMuon_Fake", "#eta vs pT (Fake);p_{T};#eta", kTH2F, {axisPt, axisEta});

    // ======================================================
    // Chi2 vs Purity for different pT ranges
    // ======================================================
    // pT 0-1 GeV/c
    histos.add("Chi2Purity/hMatchChi2_0to1_ALL", "MCH-MFT Match #chi^{2} (pT 0-1 GeV/c, ALL);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hMatchChi2_0to1_Fake", "MCH-MFT Match #chi^{2} (pT 0-1 GeV/c, Fake);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hChi2Norm_0to1_ALL", "Normalized MCH #chi^{2} (pT 0-1 GeV/c, ALL);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("Chi2Purity/hChi2Norm_0to1_Fake", "Normalized MCH #chi^{2} (pT 0-1 GeV/c, Fake);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    
    // pT 1-2 GeV/c
    histos.add("Chi2Purity/hMatchChi2_1to2_ALL", "MCH-MFT Match #chi^{2} (pT 1-2 GeV/c, ALL);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hMatchChi2_1to2_Fake", "MCH-MFT Match #chi^{2} (pT 1-2 GeV/c, Fake);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hChi2Norm_1to2_ALL", "Normalized MCH #chi^{2} (pT 1-2 GeV/c, ALL);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("Chi2Purity/hChi2Norm_1to2_Fake", "Normalized MCH #chi^{2} (pT 1-2 GeV/c, Fake);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});

    // pT 2-3 GeV/c
    histos.add("Chi2Purity/hMatchChi2_2to3_ALL", "MCH-MFT Match #chi^{2} (pT 2-3 GeV/c, ALL);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hMatchChi2_2to3_Fake", "MCH-MFT Match #chi^{2} (pT 2-3 GeV/c, Fake);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hChi2Norm_2to3_ALL", "Normalized MCH #chi^{2} (pT 2-3 GeV/c, ALL);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("Chi2Purity/hChi2Norm_2to3_Fake", "Normalized MCH #chi^{2} (pT 2-3 GeV/c, Fake);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});

    // pT 4-5 GeV/c
    histos.add("Chi2Purity/hMatchChi2_4to5_ALL", "MCH-MFT Match #chi^{2} (pT 4-5 GeV/c, ALL);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hMatchChi2_4to5_Fake", "MCH-MFT Match #chi^{2} (pT 4-5 GeV/c, Fake);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hChi2Norm_4to5_ALL", "Normalized MCH #chi^{2} (pT 4-5 GeV/c, ALL);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("Chi2Purity/hChi2Norm_4to5_Fake", "Normalized MCH #chi^{2} (pT 4-5 GeV/c, Fake);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});

    // pT 5-20 GeV/c
    histos.add("Chi2Purity/hMatchChi2_5to20_ALL", "MCH-MFT Match #chi^{2} (pT 5-20 GeV/c, ALL);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hMatchChi2_5to20_Fake", "MCH-MFT Match #chi^{2} (pT 5-20 GeV/c, Fake);#chi^{2}_{match};Counts", kTH1F, {axisMatchChi2});
    histos.add("Chi2Purity/hChi2Norm_5to20_ALL", "Normalized MCH #chi^{2} (pT 5-20 GeV/c, ALL);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});
    histos.add("Chi2Purity/hChi2Norm_5to20_Fake", "Normalized MCH #chi^{2} (pT 5-20 GeV/c, Fake);#chi^{2}/NDF;Counts", kTH1F, {axisChi2Norm});

    // ======================================================
    // Selected
    histos.add("Select/hPt_SelectedMuon_ALL", "p_{T} (Selected);p_{T};Counts", kTH1F, {axisPt});
    histos.add("Select/hP_SelectedMuon_ALL", "p (Selected);p;Counts", kTH1F, {axisP});
    histos.add("Select/hEta_SelectedMuon_ALL", "eta (Selected);#eta;Counts", kTH1F, {axisEta});
    histos.add("Select/hPhi_SelectedMuon_ALL", "phi (Selected);#phi;Counts", kTH1F, {axisPhi});
    histos.add("Select/hChi2Norm_SelectedMuon_ALL", "Norm #chi^{2} (Selected);#chi^{2}/NDF", kTH1F, {axisChi2Norm});
    histos.add("Select/hMatchChi2_SelectedMuon_ALL", "MCH-MFT Match #chi^{2} (Selected);#chi^{2}_{match}", kTH1F, {axisMatchChi2});
    histos.add("Select/hPDCA_SelectedMuon_ALL", "p#timesDCA (Selected);p#timesDCA", kTH1F, {axisPDCA});

    // True
    histos.add("Select/hPt_SelectedMuon_True", "p_{T} (Selected/True);p_{T};Counts", kTH1F, {axisPt});
    histos.add("Select/hP_SelectedMuon_True", "p (Selected/True);p;Counts", kTH1F, {axisP});
    histos.add("Select/hEta_SelectedMuon_True", "eta (Selected/True);#eta;Counts", kTH1F, {axisEta});
    histos.add("Select/hPhi_SelectedMuon_True", "phi (Selected/True);#phi;Counts", kTH1F, {axisPhi});
    histos.add("Select/hChi2Norm_SelectedMuon_True", "Norm #chi^{2} (Selected/True);#chi^{2}/NDF", kTH1F, {axisChi2Norm});
    histos.add("Select/hMatchChi2_SelectedMuon_True", "MCH-MFT Match #chi^{2} (Selected/True);#chi^{2}_{match}", kTH1F, {axisMatchChi2});
    histos.add("Select/hPDCA_SelectedMuon_True", "p#timesDCA (Selected/True);p#timesDCA", kTH1F, {axisPDCA});
    // Fake
    histos.add("Select/hPt_SelectedMuon_Fake", "p_{T} (Selected/Fake);p_{T};Counts", kTH1F, {axisPt});
    histos.add("Select/hP_SelectedMuon_Fake", "p (Selected/Fake);p;Counts", kTH1F, {axisP});
    histos.add("Select/hEta_SelectedMuon_Fake", "eta (Selected/Fake);#eta;Counts", kTH1F, {axisEta});
    histos.add("Select/hPhi_SelectedMuon_Fake", "phi (Selected/Fake);#phi;Counts", kTH1F, {axisPhi});
    histos.add("Select/hChi2Norm_SelectedMuon_Fake", "Norm #chi^{2} (Selected/Fake);#chi^{2}/NDF", kTH1F, {axisChi2Norm});
    histos.add("Select/hMatchChi2_SelectedMuon_Fake", "MCH-MFT Match #chi^{2} (Selected/Fake);#chi^{2}_{match}", kTH1F, {axisMatchChi2});
    histos.add("Select/hPDCA_SelectedMuon_Fake", "p#timesDCA (Selected/Fake);p#timesDCA", kTH1F, {axisPDCA}); 
    // ======================================================
  }

  void process(soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks)
  {
    const int fakeBit = 0x1 << 7; 

    for (auto& trk : tracks) {
      // ======================================= 
      // MatchedQUalityCuts
      // =======================================
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
      // ======================================

      float valPt = trk.pt();
      float valEta = trk.eta();
      float valPhi = trk.phi();
      float valP = valPt * std::cosh(valEta);
      float valMatchChi2 = trk.chi2MatchMCHMFT();

      // Chi2 Norm
      float valChi2Norm = -1.0;
      int ndf = 2 * trk.nClusters() - 5;
      if (ndf > 0) {
          valChi2Norm = trk.chi2() / ndf;
      }

      bool isFake = (trk.mcMask() & fakeBit) != 0;

      // ======================================================
      // Fill ALL
      // ======================================================
      histos.fill(HIST("hPt_GlobalMuon_ALL"), valPt);
      histos.fill(HIST("hP_GlobalMuon_ALL"), valP);
      histos.fill(HIST("hEta_GlobalMuon_ALL"), valEta);
      histos.fill(HIST("hPhi_GlobalMuon_ALL"), valPhi);
      if (valChi2Norm > 0) histos.fill(HIST("hChi2Norm_GlobalMuon_ALL"), valChi2Norm);
      histos.fill(HIST("hMatchChi2_GlobalMuon_ALL"), valMatchChi2);
      histos.fill(HIST("hPDCA_GlobalMuon_ALL"), pDca);

      // ======================================================
      // Fill Fake or True 
      // ======================================================
      if (isFake) {
          // --- FAKE ---
          histos.fill(HIST("hPt_GlobalMuon_Fake"), valPt);
          histos.fill(HIST("hP_GlobalMuon_Fake"), valP);
          histos.fill(HIST("hEta_GlobalMuon_Fake"), valEta);
          histos.fill(HIST("hPhi_GlobalMuon_Fake"), valPhi);
          if (valChi2Norm > 0) histos.fill(HIST("hChi2Norm_GlobalMuon_Fake"), valChi2Norm);
          histos.fill(HIST("hMatchChi2_GlobalMuon_Fake"), valMatchChi2);
          histos.fill(HIST("hPDCA_GlobalMuon_Fake"), pDca);

          histos.fill(HIST("hEtaPt_GlobalMuon_Fake"), valPt, valEta);

      } else {
          // --- TRUE ---
          histos.fill(HIST("hPt_GlobalMuon_True"), valPt);
          histos.fill(HIST("hP_GlobalMuon_True"), valP);
          histos.fill(HIST("hEta_GlobalMuon_True"), valEta);
          histos.fill(HIST("hPhi_GlobalMuon_True"), valPhi);
          if (valChi2Norm > 0) histos.fill(HIST("hChi2Norm_GlobalMuon_True"), valChi2Norm);
          histos.fill(HIST("hMatchChi2_GlobalMuon_True"), valMatchChi2);
          histos.fill(HIST("hPDCA_GlobalMuon_True"), pDca);

          histos.fill(HIST("hEtaPt_GlobalMuon_True"), valPt, valEta);
      }
      
      // ======================================================
      // Fill Chi2 vs Purity histograms by pT range
      // ======================================================
      // pT 0-1 GeV/c
      if (valPt >= 0.0 && valPt < 1.0) {
          histos.fill(HIST("Chi2Purity/hMatchChi2_0to1_ALL"), valMatchChi2);
          if (isFake) histos.fill(HIST("Chi2Purity/hMatchChi2_0to1_Fake"), valMatchChi2);
          histos.fill(HIST("Chi2Purity/hChi2Norm_0to1_ALL"), valChi2Norm);
          if (isFake) histos.fill(HIST("Chi2Purity/hChi2Norm_0to1_Fake"), valChi2Norm);
      }
      // pT 1-2 GeV/c
      if (valPt >= 1.0 && valPt < 2.0) {
          histos.fill(HIST("Chi2Purity/hMatchChi2_1to2_ALL"), valMatchChi2);
          if (isFake) histos.fill(HIST("Chi2Purity/hMatchChi2_1to2_Fake"), valMatchChi2);
          histos.fill(HIST("Chi2Purity/hChi2Norm_1to2_ALL"), valChi2Norm);
          if (isFake) histos.fill(HIST("Chi2Purity/hChi2Norm_1to2_Fake"), valChi2Norm);
      }
      // pT 2-3 GeV/c
      if (valPt >= 2.0 && valPt < 3.0) {
          histos.fill(HIST("Chi2Purity/hMatchChi2_2to3_ALL"), valMatchChi2);
          if (isFake) histos.fill(HIST("Chi2Purity/hMatchChi2_2to3_Fake"), valMatchChi2);
          histos.fill(HIST("Chi2Purity/hChi2Norm_2to3_ALL"), valChi2Norm);
          if (isFake) histos.fill(HIST("Chi2Purity/hChi2Norm_2to3_Fake"), valChi2Norm);
      }
      // pT 4-5 GeV/c
      if (valPt >= 4.0 && valPt < 5.0) {
          histos.fill(HIST("Chi2Purity/hMatchChi2_4to5_ALL"), valMatchChi2);
          if (isFake) histos.fill(HIST("Chi2Purity/hMatchChi2_4to5_Fake"), valMatchChi2);
          histos.fill(HIST("Chi2Purity/hChi2Norm_4to5_ALL"), valChi2Norm);
          if (isFake) histos.fill(HIST("Chi2Purity/hChi2Norm_4to5_Fake"), valChi2Norm);
      }
      // pT 5-20 GeV/c
      if (valPt >= 5.0 && valPt < 20.0) {
          histos.fill(HIST("Chi2Purity/hMatchChi2_5to20_ALL"), valMatchChi2);
          if (isFake) histos.fill(HIST("Chi2Purity/hMatchChi2_5to20_Fake"), valMatchChi2);
          histos.fill(HIST("Chi2Purity/hChi2Norm_5to20_ALL"), valChi2Norm);
          if (isFake) histos.fill(HIST("Chi2Purity/hChi2Norm_5to20_Fake"), valChi2Norm);
      }
      // =======================================================
      // Select
      // =======================================================
      // --- Cut parameters ---
      bool isGoodMCH = (valChi2Norm > 0 && valChi2Norm < maxChi2Norm);
      bool isGoodMatch = (valMatchChi2 < maxMatchChi2);
      if (pDca > maxPDCA) continue; 

      if (isGoodMCH && isGoodMatch) {
          histos.fill(HIST("Select/hPt_SelectedMuon_ALL"), valPt);
          histos.fill(HIST("Select/hP_SelectedMuon_ALL"), valP);
          histos.fill(HIST("Select/hEta_SelectedMuon_ALL"), valEta);
          histos.fill(HIST("Select/hPhi_SelectedMuon_ALL"), valPhi);
          histos.fill(HIST("Select/hChi2Norm_SelectedMuon_ALL"), valChi2Norm);
          histos.fill(HIST("Select/hMatchChi2_SelectedMuon_ALL"), valMatchChi2);
          histos.fill(HIST("Select/hPDCA_SelectedMuon_ALL"), pDca);

          if (isFake) {
              histos.fill(HIST("Select/hPt_SelectedMuon_Fake"), valPt);
              histos.fill(HIST("Select/hP_SelectedMuon_Fake"), valP);
              histos.fill(HIST("Select/hEta_SelectedMuon_Fake"), valEta);
              histos.fill(HIST("Select/hPhi_SelectedMuon_Fake"), valPhi);
              histos.fill(HIST("Select/hChi2Norm_SelectedMuon_Fake"), valChi2Norm);
              histos.fill(HIST("Select/hMatchChi2_SelectedMuon_Fake"), valMatchChi2);
              histos.fill(HIST("Select/hPDCA_SelectedMuon_Fake"), pDca);
          } else {
              histos.fill(HIST("Select/hPt_SelectedMuon_True"), valPt);
              histos.fill(HIST("Select/hP_SelectedMuon_True"), valP);
              histos.fill(HIST("Select/hEta_SelectedMuon_True"), valEta);
              histos.fill(HIST("Select/hPhi_SelectedMuon_True"), valPhi);
              histos.fill(HIST("Select/hChi2Norm_SelectedMuon_True"), valChi2Norm);
              histos.fill(HIST("Select/hMatchChi2_SelectedMuon_True"), valMatchChi2);
              histos.fill(HIST("Select/hPDCA_SelectedMuon_True"), pDca);
          }
      }
      // ======================================================
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<FakeMuonComparator>(cfg)
  };
}