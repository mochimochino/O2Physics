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

#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h" // for accessing to o2::ccdb
#include "Field/MagneticField.h" // for accessing to magnetic field
#include "DetectorsBase/GeometryManager.h" // GeometryManager
#include "DetectorsBase/Propagator.h"

// ML
#include "PWGDQ/Core/MuonMatchingMlResponse.h"

#include <unordered_map> 

#include <PWGDQ/Core/VarManager.h>
#include <ReconstructionDataFormats/TrackFwd.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;


// True or Fake counter
using MCHMuons = soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov, o2::aod::McFwdTrackLabels>;


// using MCHMuons = soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;
//using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
// For skimBestMuonMatchesML
// using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;
// using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
// using MFTTrackLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;


// https://aliceo2group.github.io/analysis-framework/docs/datamodel/joinsAndIterators.html
using MFTTracks = o2::aod::MFTTracks;
using MFTCovs = o2::aod::MFTTracksCov;
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

 struct Counts { int mft = 0; int mch = 0; };
  AxisSpec axisPtMFT{100, 0.0, 10.0, "MFT p_{T} [GeV/c]"};
  AxisSpec axisPtMCH{100, 0.0, 10.0, "MCH p_{T} [GeV/c]"};
  AxisSpec axisPtDiff{100, 0.0, 5.0, "p_{T} difference (MCH - MFT) [GeV/c]"};
  AxisSpec axisP{100, 0.0, 100.0, "Total Momentum p [GeV/c]"};
  
  // Position and angle parameters
  AxisSpec axisX{200, -50.0, 50.0, "X position [cm]"};
  AxisSpec axisY{200, -50.0, 50.0, "Y position [cm]"};
  AxisSpec axisXRabs{200, -100.0, 100.0, "X position (Absorber/Matching) [cm]"};
  AxisSpec axisYRabs{200, -100.0, 100.0, "Y position (Absorber/Matching) [cm]"};

  AxisSpec axisZ{200, -100.0, 100.0, "Z position [cm]"};
  AxisSpec axisPhi{64, -M_PI, M_PI, "#phi [rad]"};
  AxisSpec axisEta{100, -4.5, -2.0, "#eta"};
  AxisSpec axisTanl{100, -1.0, 1.0, "tan(#lambda)"};
  AxisSpec axisCharge{3, -1.5, 1.5, "Charge q/p_{T}"};

  AxisSpec axisZ_Absorber{200, -600.0, -400.0, "Z position (Absorber/Matching) [cm]"}; 
  AxisSpec axisZ_Vertex{200, -30.0, 30.0, "Z position (Vertex/DCA) [cm]"};

  // Error parameters
  AxisSpec axisErrPtRel{100, 0.0, 0.5, "#sigma_{p_{T}} / p_{T}"}; // Relative pT error
  AxisSpec axisErrXY{100, 0.0, 2.0, "#sigma_{xy} [cm]"}; // Position error
  
  // Check number of tracks
  AxisSpec axisTracksMFT{300, 0.0, 300.0, "Number of MFT Tracks"};
  AxisSpec axisTracksMCH{20, 0.0, 20.0, "Number of MCH Tracks"};

  // Residuals
  AxisSpec axisDeltaX{500, -50, 50, "Delta X [cm]"};
  AxisSpec axisDeltaY{500, -50, 50, "Delta Y [cm]"};

  // Fake or true counts for check
  AxisSpec axisCountsTrueorFake{3, -1.5, 1.5, "Counts (True = 1, Fake = -1)"};

  // Number of Track types
  AxisSpec axisTrackTypes{5, -0.5, 4.5, "Track Types (0: MFT-MCH-MID, 1: none, 2: MFT-MCH, 3: MCH-MID, 4: MCH)"};
  AxisSpec axisNTracksMFT{3, -1.5, 1.5, "Number of MFT Tracks"};

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
    // histos.add("Pt/pTMFT_check", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Pt/pTMCH_kToMatching", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    

    histos.add("Pt/pTMFT", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Pt/pTMCH", "MCH Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Pt/pT_MFTMCH_2D", "MFT vs MCH Tracks pT; MFT pT [GeV/c]; MCH pT [GeV/c]", kTH2F, {axisPtMFT, axisPtMCH});


    histos.add("MFT/X", "MFT TrackX", kTH1F, {axisX});
    histos.add("MFT/Y", "MFT TrackY", kTH1F, {axisY});
    histos.add("MFT/Z", "MFT TrackZ", kTH1F, {axisZ});

    histos.add("MFT/XY", "MFT Track XY", kTH2F, {axisX, axisY});
  
    
    // kToMatching
    histos.add("kToMatching/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToMatching/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("kToMatching/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("kToMatching/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ});
    histos.add("kToMatching/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    histos.add("kToMatching/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ});
    histos.add("kToMatching/XY-pT", "MCH Track pT vs XY at Matching Plane; X [cm]; Y[cm]; pT [GeV/c]", kTH3F, {axisX, axisY, axisPtMCH});
    histos.add("kToMatching/XY-p", "MCH Track p vs XY at Matching Plane; X [cm]; Y[cm]; p [GeV/c]", kTH3F, {axisX, axisY, axisP});

    histos.add("kToMatching/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("kToMatching/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("kToMatching/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToMatching/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // Error
    histos.add("kToMatching/ErrPtRel", "Relative pT Error; #sigma_{p_{T}}/p_{T}", kTH1F, {axisErrPtRel});
    histos.add("kToMatching/ErrX", "X Position Error; #sigma_{X} [cm]", kTH1F, {axisErrXY});
    histos.add("kToMatching/ErrY", "Y Position Error; #sigma_{Y} [cm]", kTH1F, {axisErrXY});

    // k kToVertex
    histos.add("kToVertex/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToVertex/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("kToVertex/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("kToVertex/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ_Vertex});
    histos.add("kToVertex/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    histos.add("kToVertex/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ_Vertex});
    histos.add("kToVertex/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("kToVertex/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("kToVertex/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToVertex/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // Error
    histos.add("kToVertex/ErrPtRel", "Relative pT Error; #sigma_{p_{T}}/p_{T}", kTH1F, {axisErrPtRel});
    histos.add("kToVertex/ErrX", "X Position Error; #sigma_{X} [cm]", kTH1F, {axisErrXY});
    histos.add("kToVertex/ErrY", "Y Position Error; #sigma_{Y} [cm]", kTH1F, {axisErrXY});


    // k kToDCA
    histos.add("kToDCA/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToDCA/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("kToDCA/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("kToDCA/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ_Vertex});
    histos.add("kToDCA/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    histos.add("kToDCA/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ_Vertex});
    histos.add("kToDCA/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("kToDCA/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("kToDCA/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToDCA/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // Error
    histos.add("kToDCA/ErrPtRel", "Relative pT Error; #sigma_{p_{T}}/p_{T}", kTH1F, {axisErrPtRel});
    histos.add("kToDCA/ErrX", "X Position Error; #sigma_{X} [cm]", kTH1F, {axisErrXY});
    histos.add("kToDCA/ErrY", "Y Position Error; #sigma_{Y} [cm]", kTH1F, {axisErrXY});

    // kToRabs
    histos.add("kToRabs/pTMCH", "MCH Track pT at After Absorber (z = -505 cm); pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToRabs/X_MCH", "MCH Track X at After Absorber (z = -505 cm); X [cm]; Counts", kTH1F, {axisXRabs});
    histos.add("kToRabs/Y_MCH", "MCH Track Y at After Absorber (z = -505 cm); Y [cm]; Counts", kTH1F, {axisYRabs});
    histos.add("kToRabs/Z_MCH", "MCH Track Z at After Absorber (z = -505 cm); Z [cm]; Counts", kTH1F, {axisZ_Absorber});
    histos.add("kToRabs/X-Y_MCH", "MCH Track X vs Y at After Absorber (z = -505 cm); X [cm]; Y[cm]", kTH2F, {axisXRabs, axisYRabs});
    histos.add("kToRabs/XYZ", "3D Position at After Absorber (z = -505 cm);X [cm];Y [cm];Z [cm]", kTH3F, {axisXRabs, axisYRabs, axisZ_Absorber});
    histos.add("kToRabs/pMCH", "MCH Total Momentum at After Absorber (z = -505 cm); p [GeV/c]", kTH1F, {axisP});
    histos.add("kToRabs/EtaMCH", "MCH #eta at After Absorber (z = -505 cm); #eta", kTH1F, {axisEta});
    histos.add("kToRabs/PhiMCH", "MCH #phi at After Absorber (z = -505 cm); #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToRabs/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // Error
    histos.add("kToRabs/ErrPtRel", "Relative pT Error; #sigma_{p_{T}}/p_{T}", kTH1F, {axisErrPtRel});
    histos.add("kToRabs/ErrX", "X Position Error; #sigma_{X} [cm]", kTH1F, {axisErrXY});
    histos.add("kToRabs/ErrY", "Y Position Error; #sigma_{Y} [cm]", kTH1F, {axisErrXY});
    // =====================================
    // For Check
    // k ToMatching
    histos.add("ToMatching/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("ToMatching/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("ToMatching/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("ToMatching/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ});
    histos.add("ToMatching/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    histos.add("ToMatching/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ});
    histos.add("ToMatching/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("ToMatching/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("ToMatching/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("ToMatching/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // Error
    histos.add("ToMatching/ErrPtRel", "Relative pT Error; #sigma_{p_{T}}/p_{T}", kTH1F, {axisErrPtRel});
    histos.add("ToMatching/ErrX", "X Position Error; #sigma_{X} [cm]", kTH1F, {axisErrXY});
    histos.add("ToMatching/ErrY", "Y Position Error; #sigma_{Y} [cm]", kTH1F, {axisErrXY});
    // =====================================

    // Check number of tracks
    histos.add("NumberOfTracks/TracksMFT", "Number of MFT Tracks;N_{tracks}", kTH1F, {axisTracksMFT});
    histos.add("NumberOfTracks/TracksMCH", "Number of MCH Tracks;N_{tracks}", kTH1F, {axisTracksMCH});
    histos.add("NumberOfTracks/TracksMFT-MCH", "2D", kTH2F, {axisTracksMFT, axisTracksMCH});

    // ressidual analysis
    histos.add("Residuals/DeltaX_MFT-MCH_kToMatching", "Delta X at Matching Plane; #Delta X [cm]; Counts", kTH1F, {axisDeltaX});
    histos.add("Residuals/DeltaY_MFT-MCH_kToMatching", "Delta Y at Matching Plane; #Delta Y [cm]; Counts", kTH1F, {axisDeltaY});

    // True of Fake counter
    histos.add("Check/CountsTrueFake", "True vs Fake Tracks; Status (True=1, Fake=-1); Counts", kTH1F, {axisCountsTrueorFake});
    // Number of Track types
    histos.add("Check/CountTrackTypeswithoutCut", "Track Types (0: MFT-MCH-MID, 1: none, 2: MFT-MCH, 3: MCH-MID, 4: MCH)", kTH1F, {axisTrackTypes});
    histos.add("Check/CountTrackTypeswithCut", "Track Types (0: MFT-MCH-MID, 1: none, 2: MFT-MCH, 3: MCH-MID, 4: MCH)", kTH1F, {axisTrackTypes});
    histos.add("Check/CountMFTTrackswithCut", "Number of MFT Tracks per Event", kTH1F, {axisTracksMFT});
    histos.add("Check/CountMFTTrackswithoutCut", "Number of MFT Tracks per Event", kTH1F, {axisNTracksMFT});

    // Matching Analysis
    // kToMatching
    histos.add("Matching/kToMatching/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Matching/kToMatching/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("Matching/kToMatching/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("Matching/kToMatching/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ});
    histos.add("Matching/kToMatching/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    histos.add("Matching/kToMatching/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ});
    histos.add("Matching/kToMatching/XY-pT", "MCH Track pT vs XY at Matching Plane; X [cm]; Y[cm]; pT [GeV/c]", kTH3F, {axisX, axisY, axisPtMCH});
    histos.add("Matching/kToMatching/XY-p", "MCH Track p vs XY at Matching Plane; X [cm]; Y[cm]; p [GeV/c]", kTH3F, {axisX, axisY, axisP});
    histos.add("Matching/kToMatching/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("Matching/kToMatching/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("Matching/kToMatching/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("Matching/kToMatching/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // kToRabs
    histos.add("Matching/kToRabs/pTMCH", "MCH Track pT at After Absorber (z = -505 cm); pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Matching/kToRabs/X_MCH", "MCH Track X at After Absorber (z = -505 cm); X [cm]; Counts", kTH1F, {axisXRabs});
    histos.add("Matching/kToRabs/Y_MCH", "MCH Track Y at After Absorber (z = -505 cm); Y [cm]; Counts", kTH1F, {axisYRabs});
    histos.add("Matching/kToRabs/Z_MCH", "MCH Track Z at After Absorber (z = -505 cm); Z [cm]; Counts", kTH1F, {axisZ_Absorber});
    histos.add("Matching/kToRabs/X-Y_MCH", "MCH Track X vs Y at After Absorber (z = -505 cm); X [cm]; Y[cm]", kTH2F, {axisXRabs, axisYRabs});
    histos.add("Matching/kToRabs/XYZ", "3D Position at After Absorber (z = -505 cm);X [cm];Y [cm];Z [cm]", kTH3F, {axisXRabs, axisYRabs, axisZ_Absorber});
    histos.add("Matching/kToRabs/pMCH", "MCH Total Momentum at After Absorber (z = -505 cm); p [GeV/c]", kTH1F, {axisP});
    histos.add("Matching/kToRabs/EtaMCH", "MCH #eta at After Absorber (z = -505 cm); #eta", kTH1F, {axisEta});
    histos.add("Matching/kToRabs/PhiMCH", "MCH #phi at After Absorber (z = -505 cm); #phi [rad]", kTH1F, {axisPhi});
    histos.add("Matching/kToRabs/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});
    // VarManager

    // Purity
    histos.add("Purity/Total_vs_MFTMult", "Total Global Muons vs MFT Multiplicity;MFT Multiplicity;Counts", kTH1F, {axisTracksMFT});
    histos.add("Purity/True_vs_MFTMult", "True Global Muons vs MFT Multiplicity;MFT Multiplicity;Counts", kTH1F, {axisTracksMFT});
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




  /*void processMatchingAnalysis(// std::map<int, Counts> const& countsMap,
                            aod::Collisions const& collisions,
                            aod::FwdTracks const& mchTracks,
                            aod::MFTTracks const& mftTracks,
                            aod::FwdTrackCovFwd const& mchCovs,
                            aod::MFTTrackCovFwd const& mftCovs
                          )*/

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

  // ====================================
  // Matching Analysis
  // ====================================
  // Using VarManager from PWGDQ/Core/VarManager.h
  void processMatchingAnalysis_likeDQ(MCHMuons const& mchJoined,
                                MyEvents const& collisions,
                                MFTTracks const& mfttracks,
                                MFTCovs const& mftCovs,
                                aod::FwdTracks const& mchTracksRaw
    )
  {
    map_mfttrackcovs.clear();

    for (auto& mfttrackConv : mftCovs) {
        map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }


    for (auto const& globalTr : mchJoined) {
      // track
      if (globalTr.trackType() != 0) continue;
      // ======================================= 
      // Acceptance cuts 
      // =======================================
      if (globalTr.eta() < -3.6 || globalTr.eta() > -2.5) continue;  // < -4.0 is right cut
      const float rAbs = globalTr.rAtAbsorberEnd();
      if (rAbs < 17.6 || rAbs > 89.5) continue; 
      const float pDca = globalTr.pDca();
      if (pDca < 0.0) continue; 
      if (rAbs < 26.5) {
          if (pDca > 594.0) continue;
      } else {
          if (pDca > 324.0) continue;
      }
      // if (globalTr.chi2() < 0.0 || globalTr.chi2() > 1e6) continue;
      // if (globalTr.chi2MatchMCHMID() < 0.0 || globalTr.chi2MatchMCHMID() > 1e6) continue; 
      // if (globalTr.chi2MatchMCHMFT() < 0.0 || globalTr.chi2MatchMCHMFT() > 1e6) continue;

      int collId = globalTr.collisionId();

      if (collId < 0 || collId >= collisions.size()) {
        continue;
      }

      auto const& collRow = collisions.iteratorAt(collId);
      // This information has Track 0 or 1 type only (including MFT track)
      int mchID = globalTr.matchMCHTrackId(); // MCH (Type 3)ID
      int mftID = globalTr.matchMFTTrackId(); // MFT ID


      // LOG(info) << "MFT Cov size: " << mftCovs.size();
      // LOG(info) << "matchMFTTrackId" << globalTr.matchMFTTrackId();
      // if (globalTr.matchMFTTrackId() >= 0) {  // To check Matched MFT track
      // int muonId = globalTr.matchMFTTrackId();
      LOG(info) << "MFT track ID: " << mftID;
      LOG(info) << "MCH track ID: " << mchID;
      // auto mchTr = mchTracksRaw.rawIteratorAt(mchID);

      auto mftTr = mfttracks.rawIteratorAt(mftID);

      // Propagate MCH tracks to each plane
      // "kToMatching", "kToVertex", "kToDCA", "kToRabs"
      /*switch (globalTr.trackType()) {
        case 3:
          LOG(info) << "MCH Track Type 3";
          auto globalTr3 = globalTr;
          break;
        default:
          LOG(info) << "Other Track Type";
          break;
      }*/
      auto propAtMP = VarManager::PropagateMuon(globalTr, collRow, VarManager::kToMatching);
      auto propAtVX = VarManager::PropagateMuon(globalTr, collRow, VarManager::kToVertex);
      auto propAtDCA = VarManager::PropagateMuon(globalTr, collRow, VarManager::kToDCA);
      auto propAtRabs = VarManager::PropagateMuon(globalTr, collRow, VarManager::kToRabs);

      double xAtMP = propAtMP.getX();
      double yAtMP = propAtMP.getY();
      double zAtMP = propAtMP.getZ();
      // LOG(info) << "MP Z: " << zAtMP << ", MP same Z: " << propAtMP.getZ();

      double phiAtMP = propAtMP.getPhi();
      double tanlAtMP = propAtMP.getTanl();
      double ptAtMP = propAtMP.getPt();

      double etaAtMP = std::asinh(tanlAtMP);
      double pAtMP = ptAtMP * std::sqrt(1.0 + tanlAtMP * tanlAtMP);

      double q2ptAtMP = propAtMP.getInvQPt(); // signed 1/pT
      double chargeAtMP = (q2ptAtMP > 0) ? 1.0 : -1.0;

      double xAtVX = propAtVX.getX();
      double yAtVX = propAtVX.getY();
      double zAtVX = propAtVX.getZ();

      double phiAtVX = propAtVX.getPhi();
      double tanlAtVX = propAtVX.getTanl();
      double ptAtVX = propAtVX.getPt();

      double etaAtVX = std::asinh(tanlAtVX);
      double pAtVX = ptAtVX * std::sqrt(1.0 + tanlAtVX * tanlAtVX);

      double q2ptAtVX = propAtVX.getInvQPt(); // signed 1/pT
      double chargeAtVX = (q2ptAtVX > 0) ? 1.0 : -1.0;

      double xAtDCA = propAtDCA.getX();
      double yAtDCA = propAtDCA.getY();
      double zAtDCA = propAtDCA.getZ();

      double phiAtDCA = propAtDCA.getPhi();
      double tanlAtDCA = propAtDCA.getTanl();
      double ptAtDCA = propAtDCA.getPt();

      double etaAtDCA = std::asinh(tanlAtDCA);
      double pAtDCA = ptAtDCA * std::sqrt(1.0 + tanlAtDCA * tanlAtDCA);

      double q2ptAtDCA = propAtDCA.getInvQPt(); // signed 1/pT
      double chargeAtDCA = (q2ptAtDCA > 0) ? 1.0 : -1.0;

      double xAtRabs = propAtRabs.getX();
      double yAtRabs = propAtRabs.getY();
      double zAtRabs = propAtRabs.getZ();

      double phiAtRabs = propAtRabs.getPhi();
      double tanlAtRabs = propAtRabs.getTanl();
      double ptAtRabs = propAtRabs.getPt();

      double etaAtRabs = std::asinh(tanlAtRabs);
      double pAtRabs = ptAtRabs * std::sqrt(1.0 + tanlAtRabs * tanlAtRabs);

      double q2ptAtRabs = propAtRabs.getInvQPt(); // signed 1/pT
      double chargeAtRabs = (q2ptAtRabs > 0) ? 1.0 : -1.0;

      //auto cov = prop.getCovariances(); 
      histos.fill(HIST("Pt/pTMCH_kToMatching"), ptAtMP);
      histos.fill(HIST("ToMatching/pMCH"), pAtMP);
      histos.fill(HIST("ToMatching/EtaMCH"), etaAtMP);
      histos.fill(HIST("ToMatching/PhiMCH"), phiAtMP);
      histos.fill(HIST("ToMatching/Charge"), chargeAtMP);
      histos.fill(HIST("ToMatching/pTMCH"), ptAtMP);
      histos.fill(HIST("ToMatching/X_MCH"), xAtMP);
      histos.fill(HIST("ToMatching/Y_MCH"), yAtMP);
      histos.fill(HIST("ToMatching/Z_MCH"), zAtMP);
      histos.fill(HIST("ToMatching/X-Y_MCH"), xAtMP, yAtMP);
      histos.fill(HIST("ToMatching/XYZ"), xAtMP, yAtMP, zAtMP);
      
      // kToMatchingPlane
      if (zAtMP == -77.5) {

        histos.fill(HIST("kToMatching/pMCH"), pAtMP);
        histos.fill(HIST("kToMatching/EtaMCH"), etaAtMP);
        histos.fill(HIST("kToMatching/PhiMCH"), phiAtMP);
        histos.fill(HIST("kToMatching/Charge"), chargeAtMP);
        histos.fill(HIST("kToMatching/pTMCH"), ptAtMP);
        histos.fill(HIST("kToMatching/X_MCH"), xAtMP);
        histos.fill(HIST("kToMatching/Y_MCH"), yAtMP);
        histos.fill(HIST("kToMatching/Z_MCH"), zAtMP);
        histos.fill(HIST("kToMatching/X-Y_MCH"), xAtMP, yAtMP);
        histos.fill(HIST("kToMatching/XYZ"), xAtMP, yAtMP, zAtMP);

        // XY vs pT and p 3D histograms
        histos.fill(HIST("kToMatching/XY-pT"), xAtMP, yAtMP, ptAtMP);
        histos.fill(HIST("kToMatching/XY-p"), xAtMP, yAtMP, pAtMP);

      }

      // kToVertex

      histos.fill(HIST("kToVertex/pMCH"), pAtVX);
      histos.fill(HIST("kToVertex/EtaMCH"), etaAtVX);
      histos.fill(HIST("kToVertex/PhiMCH"), phiAtVX);
      histos.fill(HIST("kToVertex/Charge"), chargeAtVX);
      histos.fill(HIST("kToVertex/pTMCH"), ptAtVX);
      histos.fill(HIST("kToVertex/X_MCH"), xAtVX);
      histos.fill(HIST("kToVertex/Y_MCH"), yAtVX);
      histos.fill(HIST("kToVertex/Z_MCH"), zAtVX);
      histos.fill(HIST("kToVertex/X-Y_MCH"), xAtVX, yAtVX);
      histos.fill(HIST("kToVertex/XYZ"), xAtVX, yAtVX, zAtVX);

      // kToDCA
      histos.fill(HIST("kToDCA/pMCH"), pAtDCA);
      histos.fill(HIST("kToDCA/EtaMCH"), etaAtDCA);
      histos.fill(HIST("kToDCA/PhiMCH"), phiAtDCA);
      histos.fill(HIST("kToDCA/Charge"), chargeAtDCA);
      histos.fill(HIST("kToDCA/pTMCH"), ptAtDCA);
      histos.fill(HIST("kToDCA/X_MCH"), xAtDCA);
      histos.fill(HIST("kToDCA/Y_MCH"), yAtDCA);
      histos.fill(HIST("kToDCA/Z_MCH"), zAtDCA);
      histos.fill(HIST("kToDCA/X-Y_MCH"), xAtDCA, yAtDCA);
      histos.fill(HIST("kToDCA/XYZ"), xAtDCA, yAtDCA, zAtDCA);

      // kToRabs
      if (zAtRabs == -505) {
        histos.fill(HIST("kToRabs/pMCH"), pAtRabs);
        histos.fill(HIST("kToRabs/EtaMCH"), etaAtRabs);
        histos.fill(HIST("kToRabs/PhiMCH"), phiAtRabs);
        histos.fill(HIST("kToRabs/Charge"), chargeAtRabs);
        histos.fill(HIST("kToRabs/pTMCH"), ptAtRabs);
        histos.fill(HIST("kToRabs/X_MCH"), xAtRabs);
        histos.fill(HIST("kToRabs/Y_MCH"), yAtRabs);
        histos.fill(HIST("kToRabs/Z_MCH"), zAtRabs);
        histos.fill(HIST("kToRabs/X-Y_MCH"), xAtRabs, yAtRabs);
        histos.fill(HIST("kToRabs/XYZ"), xAtRabs, yAtRabs, zAtRabs);
       }

        
        // Propagate MFT tracks to each plane
        // https://github.com/AliceO2Group/O2Physics/blob/master/PWGDQ/Core/VarManager.h#L1517-L1542
        
        bool hasCov = (map_mfttrackcovs.find(mftID) != map_mfttrackcovs.end());

        if (hasCov) {
            int32_t covIndex = map_mfttrackcovs[mftID];

            auto const& mftTrcov = mftCovs.rawIteratorAt(covIndex);

            auto mftpropAtMP = VarManager::PropagateFwd(mftTr, mftTrcov, -77.5);
            
            double xAtMP_MFT = mftpropAtMP.getX();
            double yAtMP_MFT = mftpropAtMP.getY();
            double zAtMP_MFT = mftpropAtMP.getZ();
            // LOG(info) << "MFT propagated Z;" << zAtMP_MFT;
            // LOG(info) << "MCH propagated Z;" << propAtMP.getZ();

            if (propAtMP.getZ() != zAtMP_MFT) {
              LOG(warn) << "Z at Matching Plane mismatch between MCH and MFT: MCH Z = " << propAtMP.getZ() << ", MFT Z = " << zAtMP_MFT;
            }
            double deltaX = propAtMP.getX() - mftpropAtMP.getX();
            double deltaY = propAtMP.getY() - mftpropAtMP.getY();
            // LOG(info) << "Z position for residuals" << propAtMP.getZ();
            histos.fill(HIST("Residuals/DeltaX_MFT-MCH_kToMatching"), deltaX);
            histos.fill(HIST("Residuals/DeltaY_MFT-MCH_kToMatching"), deltaY);
            // MTF
            histos.fill(HIST("MFT/X"), xAtMP_MFT);
            histos.fill(HIST("MFT/Y"), yAtMP_MFT);
            histos.fill(HIST("MFT/Z"), zAtMP_MFT);
            histos.fill(HIST("MFT/XY"), xAtMP_MFT, yAtMP_MFT);
        }
      
       
    }
  }
  // Using fwdtrackUtilities from Common/Core/fwdtrackUtilities.h
  /*template <bool withMFTCov, typename TFwdTracks, typename TMFTTracks, typename TCollision, typename TFwdTrack, typename TMFTTracksCov>
  void processMatchingAnalysis_likeEM(TCollision const& collision, TFwdTrack fwdtrack, TMFTTracksCov const& mftCovs, const bool isAmbiguous)
  {
    const auto& globalTrack = fwdtrack.template matchglobalTrack_as<TFwdTracks>(); // MCH-MID
    const auto& mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
    
    auto muonAtMP = propagateMuon(globalTrack, globalTrack, collision, propagattionPoint::kToMatchingPlane, matchingZ, mBz);
    double x = muonAtMP.getEta();

  }  */  
  template <typename MCHMuons,
            typename MyEvents,
            typename MFTTracks,
            typename MFTCovs,
            typename FwdTracks>
  void MatchingIDCompare( MCHMuons const& mchJoined,
                                MyEvents const& collisions,
                                MFTTracks const& mfttracks,
                                MFTCovs const& mftCovs,
                                FwdTracks const& mchTracksRaw)
  {
    for (auto const& globalTr : mchJoined) {
      if (globalTr.trackType() != 0) continue; // Only MFT-MCH-MID tracks (Global Muon tracks)
      int globalID = globalTr.globalIndex(); // Global Track ID
      int collID = globalTr.collisionId(); // Collision ID
      int mchTrID = globalTr.matchMCHTrackID(); // MCH ID
      int mftTrID = globalTr.matchMFTTrackId(); // MFT ID
      LOG(info) << " Global Track ID: " << globalID << ", Collision ID: " << collID << ", MCH Track ID: " << mchTrID << ", MFT Track ID: " << mftTrID;

    }
  }


  template <typename TMuons, typename TMFTTracks, typename TMFTCovs, typename TEvent>
  void skimBestMuonMatchesML(TMuons const& muons, TMFTTracks const& /*mfttracks*/, TMFTCovs const& mfCovs, TEvent const& collision)
  {
    std::unordered_map<int, std::pair<float, int>> mCandidates;
    for (const auto& muon : muons) {
      if (static_cast<int>(muon.trackType()) < 2) {
        auto muonID = muon.matchMCHTrackId();
        auto muontrack = muon.template matchMCHTrack_as<TMuons>();
        auto mfttrack = muon.template matchMFTTrack_as<TMFTTracks>();
        auto const& mfttrackcov = mfCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
        o2::track::TrackParCovFwd mftprop = VarManager::FwdToTrackPar(mfttrack, mfttrackcov);
        o2::dataformats::GlobalFwdTrack muonprop = VarManager::FwdToTrackPar(muontrack, muontrack);
        if (fConfigVariousOptions.fzMatching.value < 0.) {
          mftprop = VarManager::PropagateFwd(mfttrack, mfttrackcov, fConfigVariousOptions.fzMatching.value);
          muonprop = VarManager::PropagateMuon(muontrack, collision, VarManager::kToMatching);
        }
        std::vector<float> output;
        std::vector<float> inputML = matchingMlResponse.getInputFeaturesGlob(muon, muonprop, mftprop, collision);
        matchingMlResponse.isSelectedMl(inputML, 0, output);
        float score = output[0];
        if (mCandidates.find(muonID) == mCandidates.end()) {
          mCandidates[muonID] = {score, muon.globalIndex()};
        } else {
          if (score < mCandidates[muonID].first) {
            mCandidates[muonID] = {score, muon.globalIndex()};
          }
        }
      }
    }
    for (auto& pairCand : mCandidates) {
      fBestMatch[pairCand.second.second] = true;
    }
  }

  // ====================================
  // Purity (Number of particles dependence)
  // 横軸:MFTトラックのイベントごとの粒子数
  // 縦軸:MFT-MCH-MIDトラックのPurity
  // ====================================
  template <typename MCHMuons>
  void PurityCounter(MCHMuons const& tracks, std::map<int, Counts> const& countsMap)
  {
    for (auto const& track : tracks) {
      if (track.trackType() != 0) continue; //MFT-MCH-MID track only
      // ======================================= 
      // Acceptance cuts 
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

      // Count number of tarcks per event
      int collId = track.collisionId();
      if (countsMap.find(collId) == countsMap.end()) continue;
      float mftMult = static_cast<float>(countsMap.at(collId).mft);

      bool isTrue = (track.mcMask() == 0); // 0 is Ture
      histos.fill(HIST("Purity/Total_vs_MFTMult"), mftMult);
      if (isTrue) {
        histos.fill(HIST("Purity/True_vs_MFTMult"), mftMult);
      }

    }
  }


  void TrackTypeCounter(MCHMuons const& tracks)
  {
    for (auto const& track : tracks) {
      int trackType = track.trackType();
      histos.fill(HIST("Check/CountTrackTypeswithoutCut"), static_cast<float>(trackType));
      // ======================================= 
      // Acceptance cuts 
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
      histos.fill(HIST("Check/CountTrackTypeswithCut"), static_cast<float>(trackType));
    }
  }

  void MFTTrackCounter(MFTTracks const& tracks)
  {
    for (auto const& track : tracks) {
    
      histos.fill(HIST("Check/CountMFTTrackswithoutCut"), 1.0);
      // ======================================= 
      // Acceptance cuts 
      // =======================================
      if (track.eta() < -3.6 || track.eta() > -2.5) continue;  // < -4.0 is right cut
      histos.fill(HIST("Check/CountMFTTrackswithCut"), -1.0);
    }
  }


  void TrueorFakeCounter(MCHMuons const& tracks)
  {
    // const int fakeBit = 0.0; // Bit 0: Wrongly 
   
    for (auto const& track : tracks) {
      // ======================================= 
      // Acceptance cuts 
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


      if (track.trackType() != 0) continue; // MFT-MCH-MID tracks only

      int mcLabel = track.mcMask();
      if (mcLabel == 0) {
        histos.fill(HIST("Check/CountsTrueFake"), 1.0);
      } else {
        histos.fill(HIST("Check/CountsTrueFake"), -1.0);
      }
    }
  }   



  void process(MCHMuons const& mchJoined,
               aod::Collisions const& rawcollisions,
               MyEvents const& collisions,
               aod::FwdTracks const& mchTracks,
               MFTTracks const& mftTracks,
               MFTCovs const& mftCovs,
               ExtBCs const& bcs
               // MyEventsWithMults const& collisionsML,
               // MyMuonsWithCov const& mchtrackML,
               // MFTTrackLabeled const& mfttrackML
               //aod::McFwdTrackLabel const& mchLabels,
               //aod::McMFTTrackLabel const& mftLabels,
               //aod::FwdTrackCovFwd const& mchCovs,
               //aod::MFTTrackCovFwd const& mftCovs
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

  
   auto countsMap = countTracksPerCollision(rawcollisions, mchTracks, mftTracks);

   //fillBasicsHistograms(countsMap, mchTracks, mftTracks);

   // MatchingIDCompare(mchJoined, collisions, mfttracks, mftCovs, mchTracksRaw);
   // skimBestMuonMatchesML(mchtrackML, mfttrackML, aod::MFTTracksCov const& mftCovs, collisionsML);
   PurityCounter(mchJoined, countsMap);
   TrackTypeCounter(mchJoined);
   TrueorFakeCounter(mchJoined);
   MFTTrackCounter(mftTracks);
   //processMatchingAnalysis(collisions, mchTracks, mftTracks, mchCovs, mftCovs);
   processMatchingAnalysis_likeDQ(mchJoined, collisions, mftTracks, mftCovs, mchTracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<MchMftResiduals>(cfg)
  };
}