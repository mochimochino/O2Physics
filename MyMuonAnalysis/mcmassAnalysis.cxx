#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <unordered_map>

#include <TMath.h>
#include <TVector2.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TH1I.h> 
#include <iterator>

// O2 Headers
#include "Common/DataModel/EventSelection.h" // for useing aod::EvSels
#include "Common/Core/fwdtrackUtilities.h" //for propagateMuon

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h" // for accessing to o2::ccdb
#include "Field/MagneticField.h" // for accessing to magnetic field
#include "DetectorsBase/GeometryManager.h" // GeometryManager
#include "DetectorsBase/Propagator.h"

// ML
#include "PWGDQ/Core/MuonMatchingMlResponse.h"

#include <PWGDQ/Core/VarManager.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <TLorentzVector.h> // Mass

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;


// True or Fake counter
using MCHMuons = soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov, o2::aod::McFwdTrackLabels>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

// https://aliceo2group.github.io/analysis-framework/docs/datamodel/joinsAndIterators.html
using MFTTracks = o2::aod::MFTTracks;
using MFTCovs = o2::aod::MFTTracksCov;
// For TBC
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>; //BCs: bunch crossing, Timestamps: the timestamp of a BC


// Particle Information Analysis
using ParticleInfo = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
using ParticleInfo_mft = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;



struct mcmassAnalysis {

  // =====================================
  // For use CCDB information
  // =====================================
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::ccdb::CcdbApi ccdbApi;
  o2::field::MagneticField* fieldB;
  int mRunNumber = 0;


  int mCurrentRun = -1;
  float mMagField = 0.0;
  float mBz = 0.0; // from EM/Dilepton/TableProducer/slimmerPrimaryMuon.cxx
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of hte ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // For finding matching MFTtrack and mftCov
  std::unordered_map<int64_t, int32_t> map_mfttrackcovs;

  struct : ConfigurableGroup {
    // Track related options
    Configurable<bool> fPropTrack{"cfgPropTrack", true, "Propagate tracks to primary vertex"};
    // Muon related options
    Configurable<bool> fPropMuon{"cfgPropMuon", true, "Propagate muon tracks through absorber (do not use if applying pairing)"};
    Configurable<bool> fRefitGlobalMuon{"cfgRefitGlobalMuon", true, "Correct global muon parameters"};
    Configurable<bool> fKeepBestMatch{"cfgKeepBestMatch", false, "Keep only the best match global muons in the skimming"};
    Configurable<bool> fUseML{"cfgUseML", false, "Import ONNX model from ccdb to decide which matching candidates to keep"};
    Configurable<float> fMuonMatchEtaMin{"cfgMuonMatchEtaMin", -4.0f, "Definition of the acceptance of muon tracks to be matched with MFT"};
    Configurable<float> fMuonMatchEtaMax{"cfgMuonMatchEtaMax", -2.5f, "Definition of the acceptance of muon tracks to be matched with MFT"};
    Configurable<float> fzMatching{"cfgzMatching", -77.5f, "Plane for MFT-MCH matching"};
    Configurable<std::vector<std::string>> fModelPathsCCDB{"fModelPathsCCDB", std::vector<std::string>{"Users/m/mcoquet/MLTest"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> fInputFeatures{"cfgInputFeatures", std::vector<std::string>{"chi2MCHMFT"}, "Names of ML model input features"};
    Configurable<std::vector<std::string>> fModelNames{"cfgModelNames", std::vector<std::string>{"model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  } fConfigVariousOptions;

  std::map<uint32_t, bool> fBestMatch;

  o2::analysis::MlResponseMFTMuonMatch<float> matchingMlResponse;
  // =====================================
  // Histograms
  // =====================================
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Mass (configurable bins)
  Configurable<int> nBinsMass{"nBinsMass", 100, "N bins in mass histo"};
  Configurable<float> minMass{"minMass", 0.0, "min mass in mass histo"};
  Configurable<float> maxMass{"maxMass", 6.0, "max mass in mass histo"};
  AxisSpec axisMass{nBinsMass, minMass, maxMass, "Invariant Mass [GeV/c^{2}]"};
  AxisSpec axisPairType{3, -0.5, 2.5, "Pair Type (0:FF, 1:TF, 2:TT)"};


  Configurable<int> nBinspT{"nBinspT", 100, "N bins in pT histo"};
  Configurable<float> minpT{"minpT", 0, "min pT in pT histo"};
  Configurable<float> maxpT{"maxpT", 10, "max pT in pT histo"};
  AxisSpec axispT{nBinspT, minpT, maxpT, "p_{T} [GeV/c]"};

  // 親粒子見る？
  // 生成点を見るのもあり？
  // chi2もみるべき？



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

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    VarManager::SetDefaultVarNames(); // Important that this is called before DefineCuts() !!!

    
    mRunNumber = 0;
    mBz = 0;

    // =====================================
    // Histograms
    // =====================================
    // pT
    histos.add("SingleMuon/Pt_All", "Single Muon p_{T} - All; p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("SingleMuon/Pt_Pos", "Single Muon p_{T} - Positive; p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("SingleMuon/Pt_Neg", "Single Muon p_{T} - Negative; p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    // Mass
    histos.add("Mass/Global_Mass_Unlike", "Global Dimuon Mass (Unlike Sign) - All; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});
    histos.add("Mass/Global_Mass_Like", "Global Dimuon Mass (Like Sign) - All; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});
    histos.add("Mass/ND/Global_Mass_Unlike", "Global Dimuon Mass (Unlike Sign) nD; M_{#mu#mu} [GeV/c^{2}; Counts",  kTHnSparseD, {axisMass, axispT});


    histos.add("Mass/Global_Mass_Unlike_TT", "Global Mass (Unlike) - True-True; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});
    
    histos.add("Mass/Global_Mass_Unlike_TF", "Global Mass (Unlike) - True-Fake; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});
    
    histos.add("Mass/Global_Mass_Unlike_FF", "Global Mass (Unlike) - Fake-Fake; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});

    histos.add("Mass/Global_Mass_LikePP", "Global Dimuon Mass (Like Sign ++) ; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});
    histos.add("Mass/Global_Mass_LikeMM", "Global Dimuon Mass (Like Sign --) ; M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});


    histos.add("Dimuon/Pt_InJpsiMassRegion", "Dimuon p_{T} in J/#psi Mass Region (2.9 < M_{#mu#mu} < 3.3 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_InJpsiMassRegion_TT", "Dimuon p_{T} in J/#psi Mass Region - True-True (2.9 < M_{#mu#mu} < 3.3 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_InJpsiMassRegion_TF", "Dimuon p_{T} in J/#psi Mass Region - True-Fake (2.9 < M_{#mu#mu} < 3.3 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_InJpsiMassRegion_FF", "Dimuon p_{T} in J/#psi Mass Region - Fake-Fake (2.9 < M_{#mu#mu} < 3.3 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_LowMassRegion", "Dimuon p_{T} in Low Mass Region (0.2 < M_{#mu#mu} < 0.5 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_LowMassRegion_TF", "Dimuon p_{T} in Low Mass Region - True-Fake (0.2 < M_{#mu#mu} < 0.5 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_LowMassRegion_FF", "Dimuon p_{T} in Low Mass Region - Fake-Fake (0.2 < M_{#mu#mu} < 0.5 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
    histos.add("Dimuon/Pt_LowMassRegion_TT", "Dimuon p_{T} in Low Mass Region - True-True (0.2 < M_{#mu#mu} < 0.5 GeV/c^{2}); p_{T} [GeV/c]; Counts", kTH1F, {axispT});
  }


  // ====================================
  // To take DDDB information (not using now)
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

  // Mass Analysis for SImple Dimuon
  // ====================================
  // To Do
  // Add MC True Dimuon Analysis
  // Finding Matching relation from Mass information
  // ====================================
  void processSimpleDimuons(MCHMuons const& mchJoined)
  {
    std::map<int, std::vector<int>> collToTrackMap;
    
    for (int i = 0; i < mchJoined.size(); ++i) {
      auto const& track = mchJoined.iteratorAt(i);
      
      if (track.trackType() != 0) continue; // MFT-MCH-MID
      
      // =======================================
      // Acceptance Cuts
      // =======================================
      if (track.eta() < -3.6 || track.eta() > -2.5) continue;  // < -4.0 is right cut
      const float rAbs = track.rAtAbsorberEnd();
      if (rAbs < 17.6 || rAbs > 89.5) continue; 
      const float pDca = track.pDca();
      if (pDca < 0.0) continue; 
      if (rAbs < 26.5) {
          if (pDca > 594.0) continue;
      } else {
          if (pDca > 324.0) continue;
      }
      if (track.chi2() < 0.0 || track.chi2() > 1e6) continue;
      if (track.chi2MatchMCHMID() < 0.0 || track.chi2MatchMCHMID() > 1e6) continue; 
      if (track.chi2MatchMCHMFT() < 0.0 || track.chi2MatchMCHMFT() > 1e6) continue;

      // pTCut
      if (track.pt() < minpT) continue;
      

      histos.fill(HIST("SingleMuon/Pt_All"), track.pt());

      if (track.sign() > 0) {
           histos.fill(HIST("SingleMuon/Pt_Pos"), track.pt());
       } else {
           histos.fill(HIST("SingleMuon/Pt_Neg"), track.pt());
       }

      collToTrackMap[track.collisionId()].push_back(i);
    }

    const double mMu = 0.105658; // [GeV/c2]

    for (auto const& [collId, trackIndices] : collToTrackMap) {
      if (trackIndices.size() < 2) continue;

      for (size_t i = 0; i < trackIndices.size(); ++i) {
        for (size_t j = i + 1; j < trackIndices.size(); ++j) {
          
          auto const& tr1 = mchJoined.iteratorAt(trackIndices[i]);
          auto const& tr2 = mchJoined.iteratorAt(trackIndices[j]);

          // =============================================
          // True / Fake 
          // =============================================
          int MatchingLabel1 = tr1.mcMask();
          int MatchingLabel2 = tr2.mcMask();

          bool isTrue1 = (MatchingLabel1 == 0);
          bool isTrue2 = (MatchingLabel2 == 0);
          int pairType = 0; // 0: FF, 1: TF, 2: TT
          if (isTrue1 && isTrue2) {
            pairType = 2;
          } else if (isTrue1 || isTrue2) {
            pairType = 1;
          } else {
            pairType = 0;
          }

          float pt1 = 0, pt2 = 0;
          if (std::abs(tr1.signed1Pt()) > 1e-8) pt1 = 1.0f / std::abs(tr1.signed1Pt());
          if (std::abs(tr2.signed1Pt()) > 1e-8) pt2 = 1.0f / std::abs(tr2.signed1Pt());

          TLorentzVector v1, v2, vPair;
          v1.SetPtEtaPhiM(pt1, tr1.eta(), tr1.phi(), mMu);
          v2.SetPtEtaPhiM(pt2, tr2.eta(), tr2.phi(), mMu);
          
          vPair = v1 + v2;
          float mass = vPair.M();
          float pairPt = vPair.Pt();

          // ここで全部入れたほうが楽
          histos.fill(HIST("Mass/ND/Global_Mass_Unlike"), mass, vPair.Pt()); //?????????


          if (tr1.sign() * tr2.sign() < 0) {
 
              histos.fill(HIST("Mass/Global_Mass_Unlike"), mass);

              double massMin = 2.5;
              double massMax = 3.5;

              if (mass >= massMin && mass <= massMax) {
                  histos.fill(HIST("Dimuon/Pt_InJpsiMassRegion"), pairPt);
                        
              // さらにTT/TF/FFで分けたい場合
              if (pairType == 2) histos.fill(HIST("Dimuon/Pt_InJpsiMassRegion_TT"), pairPt);
              } else if (pairType == 0) {
                  // Fake-Fake (FF)
                  histos.fill(HIST("Dimuon/Pt_InJpsiMassRegion_FF"), pairPt);
              } else {
                  // True-Fake (TF)
                  histos.fill(HIST("Dimuon/Pt_InJpsiMassRegion_TF"), pairPt);
              }
                    
              // 別の領域例: Low Mass (< 2.0 GeV/c2)
              if (mass < 2.0) {
                  histos.fill(HIST("Dimuon/Pt_LowMassRegion"), pairPt);
                  if (pairType == 2) {
                      histos.fill(HIST("Dimuon/Pt_LowMassRegion_TT"), pairPt);
                  } else if (pairType == 0) {
                      // Fake-Fake (FF)
                      histos.fill(HIST("Dimuon/Pt_LowMassRegion_FF"), pairPt);
                  } else {
                      // True-Fake (TF)
                      histos.fill(HIST("Dimuon/Pt_LowMassRegion_TF"), pairPt);
                  }
              }
              
              if (pairType == 2) {
                  // True-True (TT)
                  histos.fill(HIST("Mass/Global_Mass_Unlike_TT"), mass);
              } else if (pairType == 0) {
                  // Fake-Fake (FF)
                  histos.fill(HIST("Mass/Global_Mass_Unlike_FF"), mass);
              } else {
                  // True-Fake (TF)
                  histos.fill(HIST("Mass/Global_Mass_Unlike_TF"), mass);
              }

          } else {
              // Like Sign
              histos.fill(HIST("Mass/Global_Mass_Like"), mass);

              if (tr1.sign() > 0) {
                histos.fill(HIST("Mass/Global_Mass_LikePP"), mass);
              } else {
                histos.fill(HIST("Mass/Global_Mass_LikeMM"), mass);
              }
          }
      }
    }
    }
  }

  void process(MCHMuons const& mchJoined,
               aod::Collisions const& rawcollisions,
               MyEvents const& collisions,
               aod::FwdTracks const& mchTracks,
               MFTTracks const& mftTracks,
               ExtBCs const& bcs,
               ParticleInfo const& mctracks,
               ParticleInfo_mft const& mcmfttracks,
               aod::McParticles const& mcParticles
               )
  {
  if (bcs.size() > 0) {
      int runNumber = bcs.begin().runNumber();
      
      if (runNumber != mCurrentRun) {
        long timestamp = bcs.begin().timestamp();
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, timestamp);
        
        if (grpmag != nullptr) {
           o2::base::Propagator::initFieldFromGRP(grpmag);
           
           VarManager::SetMagneticField(grpmag->getNominalL3Field());

           VarManager::SetMatchingPlane(-77.5); // for matching
        }
        

        VarManager::SetupMuonMagField();
        
        mCurrentRun = runNumber;
        LOG(info) << "Run " << mCurrentRun << ": Magnetic field updated for VarManager.";
      }
    }
   processSimpleDimuons(mchJoined); // mass analysis
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcmassAnalysis>(cfg)
  };
}