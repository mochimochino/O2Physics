#include <cmath>
#include <TMath.h>
#include <TVector2.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>

#include <TMath.h>
#include <TH1I.h> 
#include <iterator>
#include <map>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"

#include <PWGDQ/Core/VarManager.h>
#include <ReconstructionDataFormats/TrackFwd.h>

using namespace o2;
using namespace o2::framework;
//using namespace o2::aod;

/*using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepStd<double, 5, 5>>;
using SVector5     = ROOT::Math::SVector<double, 5>;*/

struct MchMftResiduals {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /*Configurable<float> maxphires{"maxphires", 1.0, "Max phi residual [rad]"};
  Configurable<float> minphires{"minphires", -1.0, "Min phi residual [rad]"};
  Configurable<float> maxetares{"maxetares", 1.0, "Max eta residual"};
  Configurable<float> minetares{"minetares", -1.0, "Min eta residual"};*/
  AxisSpec axisPtMFT{100, 0.0, 10.0, "MFT p_{T} [GeV/c]"};
  AxisSpec axisPtMCH{100, 0.0, 10.0, "MCH p_{T} [GeV/c]"};
  AxisSpec axisPtDiff{100, 0.0, 5.0, "p_{T} difference (MCH - MFT) [GeV/c]"};
  
  // Position and angle parameters
  AxisSpec axisX{200, -50.0, 50.0, "X position [cm]"};
  AxisSpec axisY{200, -50.0, 50.0, "Y position [cm]"};
  AxisSpec axisPhi{64, -M_PI, M_PI, "#phi [rad]"};
  AxisSpec axisTanl{100, -1.0, 1.0, "tan(#lambda)"};
  
  // Check number of tracks
  AxisSpec axisTracksMFT{300, 0.0, 300.0, "Number of MFT Tracks"};
  AxisSpec axisTracksMCH{20, 0.0, 20.0, "Number of MCH Tracks"};

  void init(InitContext const&)
  {
    //histos.add("Pt/pTMFT_check", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Pt/pTMFT", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Pt/pTMCH", "MCH Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Pt/pT_MFTMCH_2D", "MFT vs MCH Tracks pT; MFT pT [GeV/c]; MCH pT [GeV/c]", kTH2F, {axisPtMFT, axisPtMCH});
    
    // Histograms for closest pT pairs (top 2)
    histos.add("Matching/pTDiff_pair1", "pT difference for closest pair (1st); |p_{T}^{MCH} - p_{T}^{MFT}| [GeV/c]; Counts", kTH1F, {axisPtDiff});
    histos.add("Matching/pTDiff_pair2", "pT difference for 2nd closest pair; |p_{T}^{MCH} - p_{T}^{MFT}| [GeV/c]; Counts", kTH1F, {axisPtDiff});
    histos.add("Matching/pTMCH_pair1", "pT of MCH in closest pair (1st); p_{T}^{MCH} [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Matching/pTMFT_pair1", "pT of MFT in closest pair (1st); p_{T}^{MFT} [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Matching/pTMCH_pair2", "pT of MCH in 2nd closest pair; p_{T}^{MCH} [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Matching/pTMFT_pair2", "pT of MFT in 2nd closest pair; p_{T}^{MFT} [GeV/c]; Counts", kTH1F, {axisPtMFT});
    
    // Track parameters for closest pairs
    histos.add("Matching/X_MCH_pair1", "X position of MCH in closest pair; X [cm]; Counts", kTH1F, {axisX});
    histos.add("Matching/Y_MCH_pair1", "Y position of MCH in closest pair; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("Matching/Phi_MCH_pair1", "Phi of MCH in closest pair; #phi [rad]; Counts", kTH1F, {axisPhi});
    histos.add("Matching/Tanl_MCH_pair1", "Tanl of MCH in closest pair; tan(#lambda); Counts", kTH1F, {axisTanl});
    
    histos.add("Matching/X_MFT_pair1", "X position of MFT in closest pair; X [cm]; Counts", kTH1F, {axisX});
    histos.add("Matching/Y_MFT_pair1", "Y position of MFT in closest pair; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("Matching/Phi_MFT_pair1", "Phi of MFT in closest pair; #phi [rad]; Counts", kTH1F, {axisPhi});
    histos.add("Matching/Tanl_MFT_pair1", "Tanl of MFT in closest pair; tan(#lambda); Counts", kTH1F, {axisTanl});
    
    // Residuals (differences) between MCH and MFT for closest pair
    histos.add("Matching/DeltaX_pair1", "X residual (MCH - MFT) for closest pair; #Delta X [cm]; Counts", kTH1F, {axisX});
    histos.add("Matching/DeltaY_pair1", "Y residual (MCH - MFT) for closest pair; #Delta Y [cm]; Counts", kTH1F, {axisY});
    histos.add("Matching/DeltaPhi_pair1", "Phi residual (MCH - MFT) for closest pair; #Delta #phi [rad]; Counts", kTH1F, {axisPhi});
    histos.add("Matching/DeltaTanl_pair1", "Tanl residual (MCH - MFT) for closest pair; #Delta tan(#lambda); Counts", kTH1F, {axisTanl});
    
    // Check number pf tracks
    histos.add("NumberOfTracks/TracksMFT", "Number of MFT Tracks;N_{tracks}", kTH1F, {axisTracksMFT});
    histos.add("NumberOfTracks/TracksMCH", "Number of MCH Tracks;N_{tracks}", kTH1F, {axisTracksMCH});
    histos.add("NumberOfTracks/TracksMFT-MCH", "2D", kTH2F, {axisTracksMFT, axisTracksMCH});
  }


/*double calculateMatchChi2(const SVector5& mchPars, const SMatrix55Sym& mchCov,
                          const SVector5& mftPars, const SMatrix55Sym& mftCov)
{
  SMatrix55Sym covSum = mftCov + mchCov;
  if (!covSum.Invert()) {
      return 1e10; // return a large chi2 id inversion fails
  }
  SAVector5 resVec = mftPars - mchPars;
  double chi2 = ROOT::Math::Similarity(resVec, covSum);
  return chi2;
}*/


void process(aod::Collisions const& collisions,
               aod::FwdTracks const& mchTracks,
               aod::MFTTracks const& mftTracks
               /*aod::McFwdTrackLabels const& mchLabels,
               aod::McMFTTrackLabels const& mftLabels,
               aod::FwdTrackCovFwd const& mchCovs,
               aod::MFTTrackCovFwd const& mftCovs */)
  {
    // Count MFT and MCH tracks per collision.
    struct Counts { int mft = 0; int mch = 0; };
    std::map<int, Counts> trackCounts;

    // Initialize map keys for all collisions we will consider
    for (auto const& coll : collisions) {
      trackCounts[coll.globalIndex()] = Counts();
    }
    // Count MCH tracks per collision
    for (auto const& mchTr : mchTracks) {
      int mchtrack_type = mchTr.trackType();
      // Keep only standalone MUONs (type 3) as in the original code
      if (mchtrack_type != 3) {
       continue;
      }
      int collId = mchTr.collisionId();
      auto it = trackCounts.find(collId);
      if (it != trackCounts.end()) {
        it->second.mch++;
      }
      // fill overall MCH pT distribution
      //histos.fill(HIST("Pt/pTMCH"), mchTr.pt());
    }

    // Count MFT tracks per collision (apply existing selection used in this file)
    for (auto const& mftTr : mftTracks) {

      int collId = mftTr.collisionId();
      auto it = trackCounts.find(collId);
      if (it != trackCounts.end()) {
        it->second.mft++;
      }
      // fill overall MFT pT distribution
      //histos.fill(HIST("Pt/pTMFT_check"), mftTr.pt());
    }

    // Fill histograms only for collisions that have BOTH MFT and MCH tracks
    // This reduces processing for collisions that don't contain both track types.
    for (auto const& [collId, cnt] : trackCounts) {
      if (cnt.mft > 0 && cnt.mch > 0) {
        // Check number of tracks
        histos.fill(HIST("NumberOfTracks/TracksMFT"), static_cast<float>(cnt.mft));
        histos.fill(HIST("NumberOfTracks/TracksMCH"), static_cast<float>(cnt.mch));
        histos.fill(HIST("NumberOfTracks/TracksMFT-MCH"), static_cast<float>(cnt.mft), static_cast<float>(cnt.mch));
      }
    }

    // Fill selected pT histograms: for tracks that belong to collisions with both MFT and MCH
    for (auto const& mchTr : mchTracks) {
      int collId = mchTr.collisionId();
      auto it = trackCounts.find(collId);
      if (it != trackCounts.end() && it->second.mft > 0 && it->second.mch > 0) {
        int mchtrack_type = mchTr.trackType();
        // Keep only standalone MUONS (type3)
        if (mchtrack_type != 3) {
          continue;
        }
        histos.fill(HIST("Pt/pTMCH"), mchTr.pt());
      }
    }
    for (auto const& mftTr : mftTracks) {
      int collId = mftTr.collisionId();
      auto it = trackCounts.find(collId);
      if (it != trackCounts.end() && it->second.mft > 0 && it->second.mch > 0) {
        histos.fill(HIST("Pt/pTMFT"), mftTr.pt());
      }
    }

    // Find closest pT pairs (top 2) for each collision with both MCH and MFT tracks
    for (auto const& [collId, cnt] : trackCounts) {
      if (cnt.mft > 0 && cnt.mch > 0) {
        // Collect MCH and MFT tracks with their full parameters for this collision
        struct TrackInfo { float pt; float x; float y; float phi; float tanl; };
        std::vector<TrackInfo> mchTracksInfo;
        std::vector<TrackInfo> mftTracksInfo;
        
        for (auto const& mchTr : mchTracks) {
          if (mchTr.collisionId() == collId) {
            int mchtrack_type = mchTr.trackType();
            if (mchtrack_type == 3) {  // Only type 3 (standalone MUON)
              mchTracksInfo.push_back({mchTr.pt(), mchTr.x(), mchTr.y(), mchTr.phi(), mchTr.tgl()});
            }
          }
        }
        
        for (auto const& mftTr : mftTracks) {
          if (mftTr.collisionId() == collId) {
            mftTracksInfo.push_back({mftTr.pt(), mftTr.x(), mftTr.y(), mftTr.phi(), mftTr.tgl()});
          }
        }
        
        // Find closest pairs (by absolute pT difference)
        // Store pairs as (pT_diff, mch_index, mft_index)
        struct Pair { float ptDiff; int mchIdx; int mftIdx; };
        std::vector<Pair> pairs;
        
        for (size_t i = 0; i < mchTracksInfo.size(); i++) {
          for (size_t j = 0; j < mftTracksInfo.size(); j++) {
            float diff = std::abs(mchTracksInfo[i].pt - mftTracksInfo[j].pt);
            pairs.push_back({diff, static_cast<int>(i), static_cast<int>(j)});
          }
        }
        
        //auto propagate_muon = VarManager::PropagateMuon(myMuonTrack, myCollision, targetEndPoint);
        // Sort by pT difference
        std::sort(pairs.begin(), pairs.end(), 
                  [](const Pair& a, const Pair& b) { return a.ptDiff < b.ptDiff; });
        
        // Fill histograms for top 2 pairs
        if (pairs.size() >= 1) {
          const Pair& p1 = pairs[0];
          const TrackInfo& mch1 = mchTracksInfo[p1.mchIdx];
          const TrackInfo& mft1 = mftTracksInfo[p1.mftIdx];
          
          histos.fill(HIST("Matching/pTDiff_pair1"), p1.ptDiff);
          histos.fill(HIST("Matching/pTMCH_pair1"), mch1.pt);
          histos.fill(HIST("Matching/pTMFT_pair1"), mft1.pt);
          
          // Fill track parameters
          histos.fill(HIST("Matching/X_MCH_pair1"), mch1.x);
          histos.fill(HIST("Matching/Y_MCH_pair1"), mch1.y);
          histos.fill(HIST("Matching/Phi_MCH_pair1"), mch1.phi);
          histos.fill(HIST("Matching/Tanl_MCH_pair1"), mch1.tanl);
          
          histos.fill(HIST("Matching/X_MFT_pair1"), mft1.x);
          histos.fill(HIST("Matching/Y_MFT_pair1"), mft1.y);
          histos.fill(HIST("Matching/Phi_MFT_pair1"), mft1.phi);
          histos.fill(HIST("Matching/Tanl_MFT_pair1"), mft1.tanl);
          
          // Fill residuals
          histos.fill(HIST("Matching/DeltaX_pair1"), mch1.x - mft1.x);
          histos.fill(HIST("Matching/DeltaY_pair1"), mch1.y - mft1.y);
          histos.fill(HIST("Matching/DeltaPhi_pair1"), TVector2::Phi_mpi_pi(mch1.phi - mft1.phi));
          histos.fill(HIST("Matching/DeltaTanl_pair1"), mch1.tanl - mft1.tanl);
        }
        
        if (pairs.size() >= 2) {
          const Pair& p2 = pairs[1];
          float diff2 = p2.ptDiff;
          float mchPt2 = mchTracksInfo[p2.mchIdx].pt;
          float mftPt2 = mftTracksInfo[p2.mftIdx].pt;
          histos.fill(HIST("Matching/pTDiff_pair2"), diff2);
          histos.fill(HIST("Matching/pTMCH_pair2"), mchPt2);
          histos.fill(HIST("Matching/pTMFT_pair2"), mftPt2);
        }
      }
    }
  }
/*
      // x,y,phi,tanl,q/pT + Id vector(MCH track)
      SVector5 mchParams(mchTr.x(), mchTr.y(), mchTr.phi(), mchTr.tgl(), mchTr.signed1Pt());
      SMatrix55Sym mchCov;
      mchCov(1,1) = mchCovs.sigmaX();
      mchCov(2,2) = mchCovs.sigmaY();
      mchCov(3,3) = mchCovs.sigmaPhi();
      mchCov(4,4) = mchCovs.sigmaTgl();
      mchCov(5,5) = mchCovs.sigma1Pt();

      SMatrix55Sym mftCov;
      mftCov(1,1) = */





      // From 02 MatchGlobalFwd
      /*
      SMatrix55Sym I = ROOT::Math::SMatrixIdentity(), H_k, V_k;
      SVector5 m_k(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                 mftTrack.getTanl(), mftTrack.getInvQPt()),
        r_k_kminus1;
      SVector5 GlobalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym GlobalMuonTrackCovariances = mchTr.getCovariances();
          V_k(0, 0) = mftTrack.getCovariances()(0, 0);
      V_k(1, 1) = mftTrack.getCovariances()(1, 1);
      V_k(2, 2) = mftTrack.getCovariances()(2, 2);
      V_k(3, 3) = mftTrack.getCovariances()(3, 3);
      V_k(4, 4) = mftTrack.getCovariances()(4, 4);
      H_k(0, 0) = 1.0;
      H_k(1, 1) = 1.0;
      H_k(2, 2) = 1.0;
      H_k(3, 3) = 1.0;
      H_k(4, 4) = 1.0; */
      // x,y,phi,tanl,q/pT + Id vector(MFT track)
      
      // 2 candidates MFT track

      // residuals of chi2 2 candidates

      // histo (MCH-MFT Chi2 residuals vs purity)





      
      /*int iMFT = mchTr.matchMFTTrackId();

      if (iMFT < 0) {
        continue;
      }

      auto const& mftTr = mftTracks.iteratorAt(iMFT);


      int mcIdMCH = -1;
      int mcIdMFT = -1;

      if (mchTr.globalIndex() < mchLabels.size()) {
          auto const& lbl = mchLabels.iteratorAt(mchTr.globalIndex());
          mcIdMCH = lbl.mcParticleId();
      }

      if (iMFT < mftLabels.size()) {
          auto const& lbl = mftLabels.iteratorAt(iMFT);
          mcIdMFT = lbl.mcParticleId();
      }

      bool isTrue = (mcIdMCH >= 0) && (mcIdMFT >= 0) && (mcIdMCH == mcIdMFT);


      float etaMCH = mchTr.eta();
      float phiMCH = mchTr.phi();
      float ptMCH  = mchTr.pt();

      float etaMFT = mftTr.eta();
      float phiMFT = mftTr.phi();

      float dEta = etaMCH - etaMFT;
      float dPhi = TVector2::Phi_mpi_pi(phiMCH - phiMFT);


      if (isTrue) {
          histos.fill(HIST("True/hResEta"), dEta);
          histos.fill(HIST("True/hResPhi"), dPhi);
          histos.fill(HIST("True/hResEtaPhi"), dEta, dPhi);
      } else {
          histos.fill(HIST("Fake/hResEta"), dEta);
          histos.fill(HIST("Fake/hResPhi"), dPhi);
          histos.fill(HIST("Fake/hResEtaPhi"), dEta, dPhi);
          histos.fill(HIST("Fake/hPt_vs_ResPhi"), ptMCH, dPhi);
      }*/
    
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<MchMftResiduals>(cfg)
  };
}