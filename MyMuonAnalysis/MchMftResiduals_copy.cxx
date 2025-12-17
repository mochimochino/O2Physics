#include <cmath>
#include <TMath.h>
#include <TVector2.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>

#include <TMath.h>
#include <TH1I.h> 
#include <iterator>
#include <map>

#include "Common/DataModel/EventSelection.h" // for useing aod::EvSels
#include "Common/Core/fwdtrackUtilities.h" //for propagateMuon

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"

#include <DataFormatsParameters/GRPMagField.h>
#include "CCDB/BasicCCDBManager.h" // for accessing to o2::ccdb
#include "Field/MagneticField.h" // for accessing to magnetic field

#include <PWGDQ/Core/VarManager.h>
#include <ReconstructionDataFormats/TrackFwd.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;


using MCHMuons = soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;
//using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
// https://aliceo2group.github.io/analysis-framework/docs/datamodel/joinsAndIterators.html

// For TBC
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>; //BCs: bunch crossing, Timestamps: the timestamp of a BC

//using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
//using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
//using MyMuonsRealignWithCov = soa::Join<aod::FwdTracksReAlign, aod::FwdTrksCovReAlign, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

//using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
//using MFTTrackLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
/*using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepStd<double, 5, 5>>;
using SVector5     = ROOT::Math::SVector<double, 5>;*/

struct MchMftResiduals {

  // =====================================
  // For use CCDB information
  // =====================================
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = 0.0;
  float mBz = 0.0; // from EM/Dilepton/TableProducer/slimmerPrimaryMuon.cxx
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::ccdb::CcdbApi ccdbApi;
  o2::field::MagneticField* fieldB;
  int mRunNumber = 0;
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of hte ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // =====================================
  // Histograms
  // =====================================
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
    // =====================================
    // For use CCDB information
    // =====================================
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
    mBz = 0;

    // =====================================
    // Histograms
    // =====================================
    // histos.add("Pt/pTMFT_check", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Pt/pTMFT", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Pt/pTMCH", "MCH Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Pt/pT_MFTMCH_2D", "MFT vs MCH Tracks pT; MFT pT [GeV/c]; MCH pT [GeV/c]", kTH2F, {axisPtMFT, axisPtMCH});
    
    // Check number of tracks
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

  /*void processMatchingAnalysis(// std::map<int, Counts> const& countsMap,
                            aod::Collisions const& collisions,
                            aod::FwdTracks const& mchTracks,
                            aod::MFTTracks const& mftTracks,
                            aod::FwdTrackCovFwd const& mchCovs,
                            aod::MFTTrackCovFwd const& mftCovs
                          )*/

  // ====================================
  // To take DDDB information
  // ====================================
    template <typename TBC>
    void initCCDB(TBC const& bc)
    {
      if (mRunNumber == bc.runNumber()) {
      return;
      }
      mRunNumber = bc.runNumber();
      std::map<std::string, std::string> metadata;
      auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
      auto ts = soreor.first;
      auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
      }
      o2::mch::TrackExtrap::setField();
      fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      // From EM to add Magnetic information
      //const double centerMFT[3] = {0, 0, -61.4};
      //o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      //mBz = field->getBz(centerMFT); // Get filed at center of MFT
      //LOGF(info, "Bz at center of MFT = %f kZG", mBz);
    }

  // ====================================
  // Matching Analysis
  // ====================================
  // Using VarManager from PWGDQ/Core/VarManager.h
  void processMatchingAnalysis_likeDQ( MCHMuons const& mchJoined,
                                //aod::Collisions const& collisions // will add MyEvent
                                MyEvents const& collisions
    )
  {
    /*for (auto const& collision : collisions) {
      double centerMFT[3] = {0, 0, -61.4};
      auto bc = event.template bc_as<TBC>();
      initCCDB(bc);
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
      LOGF(info, "Bz at center of MFT = %f kZG", Bz);

      for (auto const& mchTr = mchJoined.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex())) {
        int collId = mchTr.collisionId();
        auto const& collRow = collisions.iteratorAt(collId);
        // auto const& collision = collisions.iteratorAt(mchTr.collisionId());
        auto prop = VarManager::PropagateMuon(mchTr, collRow, VarManager::kToMatching);

        //double x = prop.getX();
        //double y = prop.getY();
        double pt = prop.getPt();
        //auto cov = prop.getCovariances(); 
        histos.fill(HIST("Pt/pTMCH"), pt);
      }
    }*/
    for (auto const& mchTr : mchJoined) {
    int collId = mchTr.collisionId();
    auto const& collRow = collisions.iteratorAt(collId);
    // auto const& collision = collisions.iteratorAt(mchTr.collisionId());
    auto prop = VarManager::PropagateMuon(mchTr, collRow, VarManager::kToMatching);

    //double x = prop.getX();
    //double y = prop.getY();
    double pt = prop.getPt();
    //auto cov = prop.getCovariances(); 
    histos.fill(HIST("Pt/pTMCH"), pt);
    }
  }  

  // Using fwdtrackUtilities from Common/Core/fwdtrackUtilities.h
  /*template <bool withMFTCov, typename TFwdTracks, typename TMFTTracks, typename TCollision, typename TFwdTrack, typename TMFTTracksCov>
  void processMatchingAnalysis_likeEM(TCollision const& collision, TFwdTrack fwdtrack, TMFTTracksCov const& mftCovs, const bool isAmbiguous)
  {
    const auto& mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    const auto& mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
    
    auto muonAtMP = propagateMuon(mchtrack, mchtrack, collision, propagattionPoint::kToMatchingPlane, matchingZ, mBz);
    double x = muonAtMP.getEta();

  }  */     



  void process(MCHMuons const& mchJoined,
               aod::Collisions const& collisions,
               MyEvents const& collision,
               aod::FwdTracks const& mchTracks,
               aod::MFTTracks const& mftTracks
               //aod::McFwdTrackLabel const& mchLabels,
               //aod::McMFTTrackLabel const& mftLabels,
               //aod::FwdTrackCovFwd const& mchCovs,
               //aod::MFTTrackCovFwd const& mftCovs
               )
  {
   auto countsMap = countTracksPerCollision(collisions, mchTracks, mftTracks);

   //fillBasicsHistograms(countsMap, mchTracks, mftTracks);

   //processMatchingAnalysis(collisions, mchTracks, mftTracks, mchCovs, mftCovs);
   processMatchingAnalysis_likeDQ(mchJoined, collision);
   // processMatchingAnalysis_likeEM();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<MchMftResiduals>(cfg)
  };
}