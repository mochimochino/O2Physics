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
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMuonsRealignWithCov = soa::Join<aod::FwdTracksReAlign, aod::FwdTrksCovReAlign, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MFTTrackLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
/*using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepStd<double, 5, 5>>;
using SVector5     = ROOT::Math::SVector<double, 5>;*/

struct MchMftResiduals {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

 struct Counts { int mft = 0; int mch = 0; };
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
    
    // Check number pf tracks
    histos.add("NumberOfTracks/TracksMFT", "Number of MFT Tracks;N_{tracks}", kTH1F, {axisTracksMFT});
    histos.add("NumberOfTracks/TracksMCH", "Number of MCH Tracks;N_{tracks}", kTH1F, {axisTracksMCH});
    histos.add("NumberOfTracks/TracksMFT-MCH", "2D", kTH2F, {axisTracksMFT, axisTracksMCH});

    // VarManager
  }
  

  std::map<int, Counts> countTracksPerCollision(aod::Collisions const& collisions,
                                                aod::FwdTracks const& mchTracks,
                                                aod::MFTTracks const& mftTracks) 
  {
    std::map<int, Counts> countsMap;
    for (auto const& coll : collisions) {
      countsMap[coll.globalIndex()] = Counts();
    }
    for (auto const& mchTr : mchTracks) {
      if (mchTr.trackType() == 3) {
        countsMap[mchTr.collisionId()].mch++;
      }
    }
    for (auto const& mftTr : mftTracks) {
      countsMap[mftTr.collisionId()].mft++;
    }
  
  for (auto const& [collId, cnt] : countsMap) {
    if (cnt.mft > 0 && cnt.mch > 0) {
      histos.fill(HIST("NumberOfTracks/TracksMFT"), static_cast<float>(cnt.mft));
      histos.fill(HIST("NumberOfTracks/TracksMCH"), static_cast<float>(cnt.mch));
      histos.fill(HIST("NumberOfTracks/TracksMFT-MCH"), static_cast<float>(cnt.mft), static_cast<float>(cnt.mch));
    }
  }
  return countsMap;
 }
/*
void fillBasicsHistograms(std::map<int, Counts> const& countsMap,
                          aod::FwdTracks const& mchTracks,
                          aod::MFTTracks const& mftTracks)
{
  for (auto const& mchTr : mchTracks) {
    auto it = countsMap.find(mchTr.collisionId());
    if (it != countsMap.end() && it->second.mft >0 && it->seecond.mch)
  }
}*/

  void processMatchingAnalysis(// std::map<int, Counts> const& countsMap,
                            aod::Collisions const& collisions,
                            aod::FwdTracks const& mchTracks,
                            aod::MFTTracks const& mftTracks,
                            aod::FwdTrackCovFwd const& mchCovs,
                            aod::MFTTrackCovFwd const& mftCovs
                          )
  {
    auto mchJoined = soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>::interator(mchTracks, mchCovs);

    for (auto const& mchTr : mchJoined) {

    auto const& collision = collisions.iteratorAt(mchTr.collisionId());
    auto prop = VarManager::PropagateMuon(mchTr, collision, VarManager::kToMatching);
    
    double x = prop.getX();
    double y = prop.getY();
    double pt = prop.getPt();
    auto cov = prop.getCovariances(); 
    }
    /*auto mchJoined = o2::soa::join(mchTracks, mchCovs);
    for (auto const& mchTr : mchJoined) {
      auto it = countsMap.find(mchTr.collisionId());
      if (it == countsMap.end() || it->second.mft == 0 || it->second.mch ==0) {
        continue;
      }
      if (mchTr.trackType() != 3) {
        continue;
      }

      auto const& collision = collisions[mchTr.collisionId()];

      MCHWrapper wrapper{mchTr.left(), mchTr.right()};
      
      // ★ここで外挿を実行
      o2::dataformats::GlobalFwdTrack propagated_muon = 
          VarManager::PropagateMuon(wrapper, collision, VarManager::kToMatching);
*/
      //auto propX = prop.getX();
      //auto propY = prop.getY();
    //}
  }         



  void process(aod::Collisions const& collisions,
               aod::FwdTracks const& mchTracks,
               aod::MFTTracks const& mftTracks,
               aod::McFwdTrackLabel const& mchLabels,
               aod::McMFTTrackLabel const& mftLabels,
               aod::FwdTrackCovFwd const& mchCovs,
               aod::MFTTrackCovFwd const& mftCovs)
  {
   auto countsMap = countTracksPerCollision(collisions, mchTracks, mftTracks);

   //fillBasicsHistograms(countsMap, mchTracks, mftTracks);

   processMatchingAnalysis(collisions, mchTracks, mftTracks, mchCovs, mftCovs);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<MchMftResiduals>(cfg)
  };
}