#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

struct mc_muon_resto_track {
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

  Configurable<int> nBinsRes{"nBinsRes", 100, "N bins in resolution histo"};
  Configurable<float> minRes{"minRes", -0.5, "min resolution"};
  Configurable<float> maxRes{"maxRes", 0.5, "max resolution"};

  const std::array<int, 4> trackTypes = {0, 2, 3, 4};

  void init(InitContext const&)
    {
        AxisSpec axisPt{nBinsPt, minPt, maxPt, "p_{T} [GeV/c]"}; 
        AxisSpec axisEta{nBinsEta, minEta, maxEta, "#eta"};
        AxisSpec axisPhi{nBinsPhi, minPhi, maxPhi, "#phi [rad]"};
        AxisSpec axisRes{nBinsRes, minRes, maxRes, "(p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen}"};
        AxisSpec axisEtaRes{100, -0.25, 0.25, "#eta_{reco} - #eta_{gen}"};
        AxisSpec axisPhiRes{100, -0.25, 0.25, "#phi_{reco} - #phi_{gen}"};

        
        histos.add("pt_resolution", "(pT_reco - pT_gen)/pT_gen", kTH1F, {axisRes});
        histos.add("eta_resolution", "eta_reco - eta_gen", kTH1F, {axisEtaRes});
        histos.add("phi_resolution", "phi_reco - phi_gen", kTH1F, {axisPhiRes});

        histos.add("ptRes_vs_pt", "Resolution vs pT_gen", kTH2F, {axisPt, axisRes});
        histos.add("etaRes_vs_pt", "Eta Resolution vs pT_gen", kTH2F, {axisPt, axisEtaRes});
        histos.add("phiRes_vs_pt", "Phi Resolution vs pT_gen", kTH2F, {axisPt, axisPhiRes});

        histos.add("pt_gen", "Generated pT (Matched to Reco)", kTH1F, {axisPt});
        histos.add("pt_reco", "Reconstructed pT (Matched to Gen)", kTH1F, {axisPt});
        histos.add("eta_gen", "Generated eta (Matched to Reco)", kTH1F, {axisEta});
        histos.add("eta_reco", "Reconstructed eta (Matched to Gen)", kTH1F, {axisEta});
        histos.add("phi_gen", "Generated phi (Matched to Reco)", kTH1F, {axisPhi});
        histos.add("phi_reco", "Reconstructed phi (Matched to Gen)", kTH1F, {axisPhi});

        histos.add("pt_gen_all", "Generated pT (All Findable)", kTH1F, {axisPt});
        histos.add("eta_gen_all", "Generated eta (All Findable)", kTH1F, {axisEta});
        histos.add("phi_gen_all", "Generated phi (All Findable)", kTH1F, {axisPhi});

        histos.add("pt_gen_all_findable", "Generated pT (All Findable Muons)", kTH1F, {axisPt});
        histos.add("eta_gen_all_findable", "Generated eta (All Findable Muons)", kTH1F, {axisEta});
        histos.add("phi_gen_all_findable", "Generated phi (All Findable Muons)", kTH1F, {axisPhi});

        histos.add("pt_gen_all_pion_findable", "Generated pT (All Findable Pions)", kTH1F, {axisPt});
        histos.add("eta_gen_all_pion_findable", "Generated eta (All Findable Pions)", kTH1F, {axisEta});
        histos.add("phi_gen_all_pion_findable", "Generated phi (All Findable Pions)", kTH1F, {axisPhi});

        histos.add("pt_gen_all_electron_findable", "Generated pT (All Findable Electrons)", kTH1F, {axisPt});
        histos.add("eta_gen_all_electron_findable", "Generated eta (All Findable Electrons)", kTH1F, {axisEta});
        histos.add("phi_gen_all_electron_findable", "Generated phi (All Findable Electrons)", kTH1F, {axisPhi});

        histos.add("pt_gen_all_kaon_findable", "Generated pT (All Findable Kaons)", kTH1F, {axisPt});
        histos.add("eta_gen_all_kaon_findable", "Generated eta (All Findable Kaons)", kTH1F, {axisEta});
        histos.add("phi_gen_all_kaon_findable", "Generated phi (All Findable Kaons)", kTH1F, {axisPhi});

        histos.add("pt_gen_all_proton_findable", "Generated pT (All Findable Protons)", kTH1F, {axisPt});
        histos.add("eta_gen_all_proton_findable", "Generated eta (All Findable Protons)", kTH1F, {axisEta});
        histos.add("phi_gen_all_proton_findable", "Generated phi (All Findable Protons)", kTH1F, {axisPhi});

        histos.add("pt_gen_all_other_findable", "Generated pT (All Other Findable Particles)", kTH1F, {axisPt});
        histos.add("eta_gen_all_other_findable", "Generated eta (All Other Findable Particles)", kTH1F, {axisEta});
        histos.add("phi_gen_all_other_findable", "Generated phi (All Other Findable Particles)", kTH1F, {axisPhi});

        for (int type : trackTypes) {
        std::string t = std::to_string(type);

        histos.add(("pt_resolution_Type" + t).c_str(), ("pT resolution (Type " + t + ")").c_str(), kTH1F, {axisRes});
        histos.add(("eta_resolution_Type" + t).c_str(), ("eta resolution (Type " + t + ")").c_str(), kTH1F, {axisEtaRes});
        histos.add(("phi_resolution_Type" + t).c_str(), ("phi resolution (Type " + t + ")").c_str(), kTH1F, {axisPhiRes});

        histos.add(("ptRes_vs_pt_Type" + t).c_str(), ("pT res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisRes});
        histos.add(("etaRes_vs_pt_Type" + t).c_str(), ("eta res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisEtaRes});
        histos.add(("phiRes_vs_pt_Type" + t).c_str(), ("phi res vs pT (Type " + t + ")").c_str(), kTH2F, {axisPt, axisPhiRes});

        histos.add(("pt_gen_Type" + t).c_str(), ("Generated pT (Matched, Type " + t + ")").c_str(), kTH1F, {axisPt});
        histos.add(("pt_reco_Type" + t).c_str(), ("Reconstructed pT (Matched, Type " + t + ")").c_str(), kTH1F, {axisPt});
        histos.add(("eta_gen_Type" + t).c_str(), ("Generated eta (Matched, Type " + t + ")").c_str(), kTH1F, {axisEta});
        histos.add(("eta_reco_Type" + t).c_str(), ("Reconstructed eta (Matched, Type " + t + ")").c_str(), kTH1F, {axisEta});
        histos.add(("phi_gen_Type" + t).c_str(), ("Generated phi (Matched, Type " + t + ")").c_str(), kTH1F, {axisPhi});
        histos.add(("phi_reco_Type" + t).c_str(), ("Reconstructed phi (Matched, Type " + t + ")").c_str(), kTH1F, {axisPhi});
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

      float ptReco = trk.pt();
      float ptGen = mc.pt();
      if (ptGen <= 0) continue;

      float etaReco = trk.eta();
      float etaGen = mc.eta();
      float phiReco = normPhi(trk.phi());
      float phiGen = normPhi(mc.phi());

      if (ptGen < minPt || ptGen > maxPt) continue;
      if (etaGen < minEta || etaGen > maxEta) continue;
      //if (phiGen < minPhi || phiGen > maxPhi) continue;

      float resPt = (ptReco - ptGen) / ptGen;
      float resEta = etaReco - etaGen;
      float rawDPhi = phiReco - phiGen;
      float resPhi = normPhi(rawDPhi);

      int type = trk.trackType();
    
      histos.fill(HIST("pt_resolution"), resPt);
      histos.fill(HIST("eta_resolution"), resEta);
      histos.fill(HIST("phi_resolution"), resPhi);
      
      histos.fill(HIST("ptRes_vs_pt"), ptGen, resPt);
      histos.fill(HIST("etaRes_vs_pt"), ptGen, resEta);
      histos.fill(HIST("phiRes_vs_pt"), ptGen, resPhi);

      // These fill "matched" gen/reco histograms
      histos.fill(HIST("pt_gen"), ptGen);
      histos.fill(HIST("pt_reco"), ptReco);
      histos.fill(HIST("eta_gen"), etaGen);
      histos.fill(HIST("eta_reco"), etaReco);
      histos.fill(HIST("phi_gen"), phiGen);
      histos.fill(HIST("phi_reco"), phiReco);

      switch (type) {
        case 0: //global MUON (MFT+MCH+MID)
          histos.fill(HIST("pt_resolution_Type0"), resPt);
          histos.fill(HIST("eta_resolution_Type0"), resEta);
          histos.fill(HIST("phi_resolution_Type0"), resPhi);
          histos.fill(HIST("ptRes_vs_pt_Type0"), ptGen, resPt);
          histos.fill(HIST("etaRes_vs_pt_Type0"), ptGen, resEta);
          histos.fill(HIST("phiRes_vs_pt_Type0"), ptGen, resPhi);
          histos.fill(HIST("pt_gen_Type0"), ptGen);
          histos.fill(HIST("pt_reco_Type0"), ptReco);
          histos.fill(HIST("eta_gen_Type0"), etaGen);
          histos.fill(HIST("eta_reco_Type0"), etaReco);
          histos.fill(HIST("phi_gen_Type0"), phiGen);
          histos.fill(HIST("phi_reco_Type0"), phiReco);
          break;

        case 2: // MFT+MCH track
          //histos.fill(HIST("pt_resolution_Type2"), resPt);
          //histos.fill(HIST("eta_resolution_Type2"), resEta);
          //histos.fill(HIST("phi_resolution_Type2"), resPhi);
          //histos.fill(HIST("ptRes_vs_pt_Type2"), ptGen, resPt);
          //histos.fill(HIST("etaRes_vs_pt_Type2"), ptGen, resEta);
          //histos.fill(HIST("phiRes_vs_pt_Type2"), ptGen, resPhi);
          //histos.fill(HIST("pt_gen_Type2"), ptGen);
          //histos.fill(HIST("pt_reco_Type2"), ptReco);
          //histos.fill(HIST("eta_gen_Type2"), etaGen);
          //histos.fill(HIST("eta_reco_Type2"), etaReco);
          //histos.fill(HIST("phi_gen_Type2"), phiGen);
          //histos.fill(HIST("phi_reco_Type2"), phiReco);
          break;

        case 3: //standalone MUON (MCH+MID)
          //histos.fill(HIST("pt_resolution_Type3"), resPt);
          //histos.fill(HIST("eta_resolution_Type3"), resEta);
          //histos.fill(HIST("phi_resolution_Type3"), resPhi);
          //histos.fill(HIST("ptRes_vs_pt_Type3"), ptGen, resPt);
          //histos.fill(HIST("etaRes_vs_pt_Type3"), ptGen, resEta);
          //histos.fill(HIST("phiRes_vs_pt_Type3"), ptGen, resPhi);
          //histos.fill(HIST("pt_gen_Type3"), ptGen);
          //histos.fill(HIST("pt_reco_Type3"), ptReco);
          //histos.fill(HIST("eta_gen_Type3"), etaGen);
          //histos.fill(HIST("eta_reco_Type3"), etaReco);
          //histos.fill(HIST("phi_gen_Type3"), phiGen);
          //histos.fill(HIST("phi_reco_Type3"), phiReco);
          break;

        case 4: //standalone MCH track
          //histos.fill(HIST("pt_resolution_Type4"), resPt);
          //histos.fill(HIST("eta_resolution_Type4"), resEta);
          //histos.fill(HIST("phi_resolution_Type4"), resPhi);
          //histos.fill(HIST("ptRes_vs_pt_Type4"), ptGen, resPt);
          //histos.fill(HIST("etaRes_vs_pt_Type4"), ptGen, resEta);
          //histos.fill(HIST("phiRes_vs_pt_Type4"), ptGen, resPhi);
          //histos.fill(HIST("pt_gen_Type4"), ptGen);
          //histos.fill(HIST("pt_reco_Type4"), ptReco);
          //histos.fill(HIST("eta_gen_Type4"), etaGen);
          //histos.fill(HIST("eta_reco_Type4"), etaReco);
          //histos.fill(HIST("phi_gen_Type4"), phiGen);
          //histos.fill(HIST("phi_reco_Type4"), phiReco);
          break;

        default:
          break;
      }
    }
    //all particles
    for (size_t i = 0; i < mcParticles.size(); ++i) {
        auto mc = mcParticles.iteratorAt(i);
        
        float ptGen = mc.pt();
        float etaGen = mc.eta();
        float phiGen = normPhi(mc.phi());

        if (ptGen < minPt || ptGen > maxPt) continue;
        if (etaGen < minEta || etaGen > maxEta) continue;
        // if (phiGen < minPhi || phiGen > maxPhi) continue;


        histos.fill(HIST("pt_gen_all"), ptGen);
        histos.fill(HIST("eta_gen_all"), etaGen);
        histos.fill(HIST("phi_gen_all"), phiGen);
    }

    //muon
    for (size_t i = 0; i < mcParticles.size(); ++i) {
        auto mc = mcParticles.iteratorAt(i);

        if (std::abs(mc.pdgCode()) != 13) continue; 
        
        float ptGen = mc.pt();
        float etaGen = mc.eta();
        float phiGen = normPhi(mc.phi());

        if (ptGen < minPt || ptGen > maxPt) continue;
        if (etaGen < minEta || etaGen > maxEta) continue;
        // if (phiGen < minPhi || phiGen > maxPhi) continue;


        histos.fill(HIST("pt_gen_all_findable"), ptGen);
        histos.fill(HIST("eta_gen_all_findable"), etaGen);
        histos.fill(HIST("phi_gen_all_findable"), phiGen);
    }

    //pion
    for (size_t i = 0; i < mcParticles.size(); ++i) {
        auto mc = mcParticles.iteratorAt(i);

        if (std::abs(mc.pdgCode()) != 211) continue;  //pion +
        
        float ptGen = mc.pt();
        float etaGen = mc.eta();
        float phiGen = normPhi(mc.phi());

        if (ptGen < minPt || ptGen > maxPt) continue;
        if (etaGen < minEta || etaGen > maxEta) continue;
        // if (phiGen < minPhi || phiGen > maxPhi) continue;


        histos.fill(HIST("pt_gen_all_pion_findable"), ptGen);
        histos.fill(HIST("eta_gen_all_pion_findable"), etaGen);
        histos.fill(HIST("phi_gen_all_pion_findable"), phiGen);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<mc_muon_resto_track>(cfg)};
}
