#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

struct fwd_dca {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable parameters pT Eta Phi
  Configurable<int> nBinsPt{"nBinsPt", 60, "N bins in pT histo"};
  Configurable<float> minPt{"minPt", 0.0, "min pT"};
  Configurable<float> maxPt{"maxPt", 6.0, "max pT"};
  Configurable<int> nBinsEta{"nBinsEta", 50, "N bins in Eta histo"};
  Configurable<float> minEta{"minEta", -4.0, "min Eta"};
  Configurable<float> maxEta{"maxEta", -2.5, "max Eta"};
  Configurable<int> nBinsPhi{"nBinsPhi", 64, "N bins in Phi histo"};
  Configurable<float> minPhi{"minPhi", - TMath::Pi(), "min Phi"};
  Configurable<float> maxPhi{"maxPhi", TMath::Pi(), "max Phi"};

  Configurable<int> nBinsRes{"nBinsRes", 20, "N bins in resolution histo"};
  Configurable<float> minRes{"minRes", -0.5, "min resolution"};
  Configurable<float> maxRes{"maxRes", 0.5, "max resolution"};

  Configurable<int> nBinsPDca{"nBinsPDca", 200, "N bins in p*Dca histo"};
  Configurable<float> minPDca{"minPDca", 0.0, "min p*Dca"};
  Configurable<float> maxPDca{"maxPDca", 1000.0, "max p*Dca"};

  const std::array<int, 4> trackTypes = {0, 2, 3, 4};

  void init(InitContext const&)
    {
        AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T} [GeV/c]"}; 
        AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
        AxisSpec axisPhi{nBinsPhi, minPhi, maxPhi, "#phi [rad]"};
        AxisSpec axisRes{nBinsRes, minRes, maxRes, "(p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen}"};
        AxisSpec axisEtaRes{100, -0.25, 0.25, "#eta_{reco} - #eta_{gen}"};
        AxisSpec axisPhiRes{100, -0.25, 0.25, "#phi_{reco} - #phi_{gen}"};
        AxisSpec axisPDca{nBinsPDca, minPDca, maxPDca, "p*DCA [GeV/c * cm]"};
        AxisSpec axisChi2Match{100, 0, 100, "Matching #chi^{2}"}; 
        AxisSpec axisRAbs{100, 10, 100, "R at Absorber End [cm]"};



        histos.add("Resolution/pt_resolution", "(pT_reco - pT_gen)/pT_gen", kTH1F, {axisRes});
        histos.add("Resolution/eta_resolution", "eta_reco - eta_gen", kTH1F, {axisEtaRes});
        histos.add("Resolution/phi_resolution", "phi_reco - phi_gen", kTH1F, {axisPhiRes});

        histos.add("Resolution/ptRes_vs_pt", "Resolution vs pT_gen", kTH2F, {axisPt, axisRes});
        histos.add("Resolution/etaRes_vs_pt", "Eta Resolution vs pT_gen", kTH2F, {axisPt, axisEtaRes});
        histos.add("Resolution/phiRes_vs_pt", "Phi Resolution vs pT_gen", kTH2F, {axisPt, axisPhiRes});

        histos.add("pt_gen", "Generated pT (Matched to Reco)", kTH1F, {axisPt});
        histos.add("pt_reco", "Reconstructed pT (Matched to Gen)", kTH1F, {axisPt});
        histos.add("eta_gen", "Generated eta (Matched to Reco)", kTH1F, {axisEta});
        histos.add("eta_reco", "Reconstructed eta (Matched to Gen)", kTH1F, {axisEta});
        histos.add("phi_gen", "Generated phi (Matched to Reco)", kTH1F, {axisPhi});
        histos.add("phi_reco", "Reconstructed phi (Matched to Gen)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all", "Generated pT (All Findable)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all", "Generated eta (All Findable)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all", "Generated phi (All Findable)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all_findable", "Generated pT (All Findable Muons)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all_findable", "Generated eta (All Findable Muons)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all_findable", "Generated phi (All Findable Muons)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all_pion_findable", "Generated pT (All Findable Pions)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all_pion_findable", "Generated eta (All Findable Pions)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all_pion_findable", "Generated phi (All Findable Pions)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all_electron_findable", "Generated pT (All Findable Electrons)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all_electron_findable", "Generated eta (All Findable Electrons)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all_electron_findable", "Generated phi (All Findable Electrons)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all_kaon_findable", "Generated pT (All Findable Kaons)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all_kaon_findable", "Generated eta (All Findable Kaons)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all_kaon_findable", "Generated phi (All Findable Kaons)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all_proton_findable", "Generated pT (All Findable Protons)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all_proton_findable", "Generated eta (All Findable Protons)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all_proton_findable", "Generated phi (All Findable Protons)", kTH1F, {axisPhi});

        histos.add("Findable/pt_gen_all_other_findable", "Generated pT (All Other Findable Particles)", kTH1F, {axisPt});
        histos.add("Findable/eta_gen_all_other_findable", "Generated eta (All Other Findable Particles)", kTH1F, {axisEta});
        histos.add("Findable/phi_gen_all_other_findable", "Generated phi (All Other Findable Particles)", kTH1F, {axisPhi});

        AxisSpec axisParentCat{6, 0.5, 6.5, "Parent Category"};
        histos.add("muon_parent_category", "Category of Muon Parent (from Reco Muon)", kTH1F, {axisParentCat});
        histos.add("pt_reco_primary", "Reco pT (Parent: Primary)", kTH1F, {axisPt});
        histos.add("pt_reco_pion", "Reco pT (Parent: Pion)", kTH1F, {axisPt});
        histos.add("pt_reco_kaon", "Reco pT (Parent: Kaon)", kTH1F, {axisPt});
        histos.add("pt_reco_charm", "Reco pT (Parent: Charm)", kTH1F, {axisPt});
        histos.add("pt_reco_beauty", "Reco pT (Parent: Beauty)", kTH1F, {axisPt});
        histos.add("pt_reco_other", "Reco pT (Parent: Other)", kTH1F, {axisPt});

        histos.add("PDCA/pdca_vs_pt_beauty", "p*DCA vs pT (Parent: Beauty)", kTH2F, {axisPt, axisPDca});
        histos.add("PDCA/pdca_vs_pt_pion", "p*DCA vs pT (Parent: Pion)", kTH2F, {axisPt, axisPDca});
        histos.add("PDCA/pdca_vs_pt_charm", "p*DCA vs pT (Parent: Charm)", kTH2F, {axisPt, axisPDca});

          // 1D pDCA
        histos.add("PDCA/pdca_all", "p*DCA (All Matched Muons)", kTH1F, {axisPDca});
        histos.add("PDCA/pdca_primary", "p*DCA (Parent: Primary)", kTH1F, {axisPDca});
        histos.add("PDCA/pdca_pion", "p*DCA (Parent: Pion)", kTH1F, {axisPDca});
        histos.add("PDCA/pdca_kaon", "p*DCA (Parent: Kaon)", kTH1F, {axisPDca});
        histos.add("PDCA/pdca_charm", "p*DCA (Parent: Charm)", kTH1F, {axisPDca});
        histos.add("PDCA/pdca_beauty", "p*DCA (Parent: Beauty)", kTH1F, {axisPDca});
        histos.add("PDCA/pdca_other", "p*DCA (Parent: Other)", kTH1F, {axisPDca});

        // 2D (vs pT)
        histos.add("PDCA/pdca_vs_pt", "p*DCA vs pT (All Matched Muons)", kTH2F, {axisPt, axisPDca});

        // 1D Inclusive
        histos.add("trk_chi2MatchMCHMFT", "MCH-MFT Match Chi2", kTH1F, {axisChi2Match});
        histos.add("trk_rAbs", "R at Absorber End", kTH1F, {axisRAbs});

        // 1D by Parent (Signal vs BG)
        histos.add("trk_chi2MatchMCHMFT_beauty", "MCH-MFT Match Chi2 (Beauty)", kTH1F, {axisChi2Match});
        histos.add("trk_chi2MatchMCHMFT_pion", "MCH-MFT Match Chi2 (Pion)", kTH1F, {axisChi2Match});
        histos.add("trk_rAbs_beauty", "R at Absorber (Beauty)", kTH1F, {axisRAbs});
        histos.add("trk_rAbs_pion", "R at Absorber (Pion)", kTH1F, {axisRAbs});
        
        // pt_reco_all
        histos.add("pt_reco_all", "Reco pT (All Matched Muons)", kTH1F, {axisPt});

        for (int type : trackTypes) {
        std::string t = std::to_string(type);
        std::string dirName = "TrackType/Type" + t + "/";

        histos.add((dirName + "pt_resolution_Type" + t).c_str(), ("pT resolution (Type " + t + ")").c_str(), kTH1F, {axisRes});
        histos.add((dirName + "eta_resolution_Type" + t).c_str(), ("eta resolution (Type " + t + ")").c_str(), kTH1F, {axisEtaRes});
        histos.add((dirName + "phi_resolution_Type" + t).c_str(), ("phi resolution (Type " + t + ")").c_str(), kTH1F, {axisPhiRes});

        histos.add((dirName + "ptRes_vs_pt_Type" + t).c_str(), ("pT res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisRes});
        histos.add((dirName + "etaRes_vs_pt_Type" + t).c_str(), ("eta res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisEtaRes});
        histos.add((dirName + "phiRes_vs_pt_Type" + t).c_str(), ("phi res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisPhiRes});

        histos.add((dirName + "pt_gen_Type" + t).c_str(), ("Generated pT (Matched, Type " + t + ")").c_str(), kTH1F, {axisPt});
        histos.add((dirName + "pt_reco_Type" + t).c_str(), ("Reconstructed pT (Matched, Type " + t + ")").c_str(), kTH1F, {axisPt});
        histos.add((dirName + "eta_gen_Type" + t).c_str(), ("Generated eta (Matched, Type " + t + ")").c_str(), kTH1F, {axisEta});
        histos.add((dirName + "eta_reco_Type" + t).c_str(), ("Reconstructed eta (Matched, Type " + t + ")").c_str(), kTH1F, {axisEta});
        histos.add((dirName + "phi_gen_Type" + t).c_str(), ("Generated phi (Matched, Type " + t + ")").c_str(), kTH1F, {axisPhi});
        histos.add((dirName + "phi_reco_Type" + t).c_str(), ("Reconstructed phi (Matched, Type " + t + ")").c_str(), kTH1F, {axisPhi});
        }
    } 

  void process(soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks,aod::McParticles const& mcParticles)
   {

    auto normPhi = [](double a) {
      while (a > TMath::Pi()) a -= 2.0 * TMath::Pi();
      while (a <= -TMath::Pi()) a += 2.0 * TMath::Pi();
      return a;
    };

    for (auto& trk : tracks) {
      int mcId = trk.mcParticleId();
      if (mcId < 0 || mcId >= static_cast<int>(mcParticles.size())) continue;

      auto mc = mcParticles.iteratorAt(mcId);
      if (std::abs(mc.pdgCode()) != 13) continue; // muon only

      // matchedQualityCuts
      if (trk.eta() < -4.0 || trk.eta() > -2.5) continue; // Match eta cut
      if (trk.rAtAbsorberEnd() < 17.6 || trk.rAtAbsorberEnd() > 89.5) continue; // Match rAtAbsorberEnd cut
      if (trk.pDca() > 594.0 && trk.rAtAbsorberEnd() < 26.5) continue; // Match pDca cut
      if (trk.pDca() > 324.0 && trk.rAtAbsorberEnd() >= 26.5) continue; // Match pDca cut
      if (trk.chi2() > 1e6) continue; // Match chi2 cut
      if (trk.chi2MatchMCHMID() > 1e6) continue; // Match chi2MatchMCHMID cut
      if (trk.chi2MatchMCHMFT() > 1e6) continue; // Match chi2MatchMCHMFT cut

      float ptReco = trk.pt();
      float ptGen = mc.pt();
      //if (ptGen <= 0) continue;

      float etaReco = trk.eta();
      float etaGen = mc.eta();
      float phiReco = normPhi(trk.phi());
      float phiGen = normPhi(mc.phi());

      //if (ptGen < minPt || ptGen > maxPt) continue;
      //if (etaGen < minEta || etaGen > maxEta) continue;
      //if (phiGen < minPhi || phiGen > maxPhi) continue;

      float resPt = (ptGen > 0) ? (ptReco - ptGen) / ptGen : 0; // Avoid division by zero
      float resEta = etaReco - etaGen;
      float rawDPhi = phiReco - phiGen;
      float resPhi = normPhi(rawDPhi);

      float pdca = trk.pDca();
      float chi2MatchMFT = trk.chi2MatchMCHMFT();
      float rAbs = trk.rAtAbsorberEnd();

      histos.fill(HIST("PDCA/pdca_all"), pdca);
      histos.fill(HIST("PDCA/pdca_vs_pt"), ptReco, pdca);

      histos.fill(HIST("trk_chi2MatchMCHMFT"), chi2MatchMFT);
      histos.fill(HIST("trk_rAbs"), rAbs);
      histos.fill(HIST("pt_reco_all"), ptReco);

      int parentCategory = 6;
      auto motherIds = mc.mothersIds();
      int motherId = -1;
      if (!motherIds.empty()) {
          motherId = motherIds[0];  
      }

      if (motherId < 0) {
          parentCategory = 1; // 1 = Primary
      } else {
          auto mother = mcParticles.iteratorAt(motherId);
          int motherPdg = std::abs(mother.pdgCode());

          bool isPion = (motherPdg == 211); // pi+/-
          // K+/- (321), K_L (130), K_S (310)
          bool isKaon = (motherPdg == 321 || motherPdg == 130 || motherPdg == 310); 
          // Charm (D meson: 4xx, Lambda_c: 4xxx etc)
          bool isCharm = (motherPdg / 100 % 10 == 4) || (motherPdg / 1000 % 10 == 4);
          // Beauty (B meson: 5xx, Lambda_b: 5xxx etc)
          bool isBeauty = (motherPdg / 100 % 10 == 5) || (motherPdg / 1000 % 10 == 5);

          if (isPion) {
              parentCategory = 2; // 2 = Pion
          } else if (isKaon) {
              parentCategory = 3; // 3 = Kaon
          } else if (isCharm) {
              parentCategory = 4; // 4 = Charm
          } else if (isBeauty) {
              parentCategory = 5; // 5 = Beauty
          }
      }

      histos.fill(HIST("muon_parent_category"), (float)parentCategory);

      switch (parentCategory) {
        case 1: // Primary
          histos.fill(HIST("pt_reco_primary"), ptReco);
          histos.fill(HIST("PDCA/pdca_primary"), pdca);
          break;
        case 2: // Pion
          histos.fill(HIST("pt_reco_pion"), ptReco);
          histos.fill(HIST("PDCA/pdca_pion"), pdca);
          histos.fill(HIST("PDCA/pdca_vs_pt_pion"), ptReco, pdca);

          histos.fill(HIST("trk_chi2MatchMCHMFT_pion"), chi2MatchMFT);
          histos.fill(HIST("trk_rAbs_pion"), rAbs);
          break;
        case 3: // Kaon
          histos.fill(HIST("pt_reco_kaon"), ptReco);
          histos.fill(HIST("PDCA/pdca_kaon"), pdca);

          break;
        case 4: // Charm
          histos.fill(HIST("pt_reco_charm"), ptReco);
          histos.fill(HIST("PDCA/pdca_charm"), pdca);
          histos.fill(HIST("PDCA/pdca_vs_pt_charm"), ptReco, pdca);
          break;
        case 5: // Beauty
          histos.fill(HIST("pt_reco_beauty"), ptReco);
          histos.fill(HIST("PDCA/pdca_beauty"), pdca);
          histos.fill(HIST("PDCA/pdca_vs_pt_beauty"), ptReco, pdca);

          histos.fill(HIST("trk_chi2MatchMCHMFT_beauty"), chi2MatchMFT);
          histos.fill(HIST("trk_rAbs_beauty"), rAbs);
          break;
        case 6: // Other
        default:
          histos.fill(HIST("pt_reco_other"), ptReco);
          histos.fill(HIST("PDCA/pdca_other"), pdca);
          break;
      }
      
      int type = trk.trackType();
      
      if (ptGen > 0) { // Avoid filling resolution histos if ptGen is 0
        histos.fill(HIST("Resolution/pt_resolution"), resPt);
        histos.fill(HIST("Resolution/ptRes_vs_pt"), ptGen, resPt);
      }
      histos.fill(HIST("Resolution/eta_resolution"), resEta);
      histos.fill(HIST("Resolution/phi_resolution"), resPhi);
      
      histos.fill(HIST("Resolution/etaRes_vs_pt"), ptGen, resEta);
      histos.fill(HIST("Resolution/phiRes_vs_pt"), ptGen, resPhi);

      // These fill "matched" gen/reco histograms (No path, matches init)
      histos.fill(HIST("pt_gen"), ptGen);
      histos.fill(HIST("pt_reco"), ptReco);
      histos.fill(HIST("eta_gen"), etaGen);
      histos.fill(HIST("eta_reco"), etaReco);
      histos.fill(HIST("phi_gen"), phiGen);
      histos.fill(HIST("phi_reco"), phiReco);


      switch (type) {
        case 0: //global MUON (MFT+MCH+MID)
          if (ptGen > 0) {
            histos.fill(HIST("TrackType/Type0/pt_resolution_Type0"), resPt);
            histos.fill(HIST("TrackType/Type0/ptRes_vs_pt_Type0"), ptGen, resPt);
          }
          histos.fill(HIST("TrackType/Type0/eta_resolution_Type0"), resEta);
          histos.fill(HIST("TrackType/Type0/phi_resolution_Type0"), resPhi);
          histos.fill(HIST("TrackType/Type0/etaRes_vs_pt_Type0"), ptGen, resEta);
          histos.fill(HIST("TrackType/Type0/phiRes_vs_pt_Type0"), ptGen, resPhi);
          histos.fill(HIST("TrackType/Type0/pt_gen_Type0"), ptGen);
          histos.fill(HIST("TrackType/Type0/pt_reco_Type0"), ptReco);
          histos.fill(HIST("TrackType/Type0/eta_gen_Type0"), etaGen);
          histos.fill(HIST("TrackType/Type0/eta_reco_Type0"), etaReco);
          histos.fill(HIST("TrackType/Type0/phi_gen_Type0"), phiGen);
          histos.fill(HIST("TrackType/Type0/phi_reco_Type0"), phiReco);
          break;

        case 2: // MFT+MCH track
         
          break;

        case 3: //standalone MUON (MCH+MID)
          if (ptGen > 0) {
            histos.fill(HIST("TrackType/Type3/pt_resolution_Type3"), resPt);
            histos.fill(HIST("TrackType/Type3/ptRes_vs_pt_Type3"), ptGen, resPt);
          }
          histos.fill(HIST("TrackType/Type3/eta_resolution_Type3"), resEta);
          histos.fill(HIST("TrackType/Type3/phi_resolution_Type3"), resPhi);
          histos.fill(HIST("TrackType/Type3/etaRes_vs_pt_Type3"), ptGen, resEta);
          histos.fill(HIST("TrackType/Type3/phiRes_vs_pt_Type3"), ptGen, resPhi);
          histos.fill(HIST("TrackType/Type3/pt_gen_Type3"), ptGen);
          histos.fill(HIST("TrackType/Type3/pt_reco_Type3"), ptReco);
          histos.fill(HIST("TrackType/Type3/eta_gen_Type3"), etaGen);
          histos.fill(HIST("TrackType/Type3/eta_reco_Type3"), etaReco);
          histos.fill(HIST("TrackType/Type3/phi_gen_Type3"), phiGen);
          histos.fill(HIST("TrackType/Type3/phi_reco_Type3"), phiReco);
          break;

        case 4: //standalone MCH track
       
          break;

        default:
          break;
      }
    }
    
    for (size_t i = 0; i < mcParticles.size(); ++i) {
        auto mc = mcParticles.iteratorAt(i);
        
        float ptGen = mc.pt();
        float etaGen = mc.eta();
        float phiGen = normPhi(mc.phi());

        if (ptGen < minPt || ptGen > maxPt) continue;
        if (etaGen < minEta || etaGen > maxEta) continue;
        // if (phiGen < minPhi || phiGen > maxPhi) continue;

        histos.fill(HIST("Findable/pt_gen_all"), ptGen);
        histos.fill(HIST("Findable/eta_gen_all"), etaGen);
        histos.fill(HIST("Findable/phi_gen_all"), phiGen);

        int pdgCode = std::abs(mc.pdgCode());

        if (pdgCode == 13) { // Muon
            histos.fill(HIST("Findable/pt_gen_all_findable"), ptGen);
            histos.fill(HIST("Findable/eta_gen_all_findable"), etaGen);
            histos.fill(HIST("Findable/phi_gen_all_findable"), phiGen);
        } else if (pdgCode == 211) { // Pion
            histos.fill(HIST("Findable/pt_gen_all_pion_findable"), ptGen);
            histos.fill(HIST("Findable/eta_gen_all_pion_findable"), etaGen);
            histos.fill(HIST("Findable/phi_gen_all_pion_findable"), phiGen);
        } else if (pdgCode == 11) { // Electron
            histos.fill(HIST("Findable/pt_gen_all_electron_findable"), ptGen);
            histos.fill(HIST("Findable/eta_gen_all_electron_findable"), etaGen);
            histos.fill(HIST("Findable/phi_gen_all_electron_findable"), phiGen);
        } else if (pdgCode == 321) { // Kaon
            histos.fill(HIST("Findable/pt_gen_all_kaon_findable"), ptGen);
            histos.fill(HIST("Findable/eta_gen_all_kaon_findable"), etaGen);
            histos.fill(HIST("Findable/phi_gen_all_kaon_findable"), phiGen);
        } else if (pdgCode == 2212) { // Proton
            histos.fill(HIST("Findable/pt_gen_all_proton_findable"), ptGen);
            histos.fill(HIST("Findable/eta_gen_all_proton_findable"), etaGen);
            histos.fill(HIST("Findable/phi_gen_all_proton_findable"), phiGen);
        } else {
            histos.fill(HIST("Findable/pt_gen_all_other_findable"), ptGen);
            histos.fill(HIST("Findable/eta_gen_all_other_findable"), etaGen);
            histos.fill(HIST("Findable/phi_gen_all_other_findable"), phiGen);
        }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<fwd_dca>(cfg)};
}