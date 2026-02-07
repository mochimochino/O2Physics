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



struct mcParticleAnalysis {

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

  struct Counts { int mft = 0; int mch = 0; };
  AxisSpec axisPtMFT{100, 0.0, 10.0, "MFT p_{T} [GeV/c]"};
  AxisSpec axisPtMCH{100, 0.0, 10.0, "MCH p_{T} [GeV/c]"};
  AxisSpec axisPtDiff{100, 0.0, 5.0, "p_{T} difference (MCH - MFT) [GeV/c]"};
  AxisSpec axisP{100, 0.0, 100.0, "Total Momentum p [GeV/c]"};
  
  // Position and angle parameters
  // [OPTIMIZATION] Binning reduced slightly to save memory
  AxisSpec axisX{100, -50.0, 50.0, "X position [cm]"};
  AxisSpec axisY{100, -50.0, 50.0, "Y position [cm]"};
  AxisSpec axisXRabs{100, -100.0, 100.0, "X position (Absorber/Matching) [cm]"};
  AxisSpec axisYRabs{100, -100.0, 100.0, "Y position (Absorber/Matching) [cm]"};

  AxisSpec axisZ{200, -100.0, 100.0, "Z position [cm]"};
  AxisSpec axisPhi{64, -M_PI, M_PI, "#phi [rad]"};
  AxisSpec axisEta{100, -4.5, -2.0, "#eta"};
  AxisSpec axisTanl{100, -1.0, 1.0, "tan(#lambda)"};
  AxisSpec axisCharge{3, -1.5, 1.5, "Charge q/p_{T}"};

  AxisSpec axisZ_Absorber{200, -600.0, -400.0, "Z position (Absorber/Matching) [cm]"}; 
  AxisSpec axisZ_Vertex{200, -30.0, 30.0, "Z position (Vertex/DCA) [cm]"};

  AxisSpec axisMatchChi2{100, 0.0, 100.0, "#chi^{2}_{Match} (MCH-MFT)"};

  // Error parameters
  AxisSpec axisErrPtRel{100, 0.0, 0.5, "#sigma_{p_{T}} / p_{T}"}; // Relative pT error
  AxisSpec axisErrXY{100, 0.0, 2.0, "#sigma_{xy} [cm]"}; // Position error
  
  // Check number of tracks
  AxisSpec axisTracksMFT{300, 0.0, 300.0, "Number of MFT Tracks"};
  AxisSpec axisTracksMCH{20, 0.0, 20.0, "Number of MCH Tracks"};

  AxisSpec axisIndex{2000, 0.0, 10000.0, "Global Index"}; // [OPTIMIZATION] Binning reduced

  // Residuals
  AxisSpec axisDeltaX{500, -50, 50, "Delta X [cm]"};
  AxisSpec axisDeltaY{500, -50, 50, "Delta Y [cm]"};

  // Fake or true counts for check
  AxisSpec axisCountsTrueorFake{3, -1.5, 1.5, "Counts (True = 1, Fake = -1)"};

  // Number of Track types
  AxisSpec axisTrackTypes{5, -0.5, 4.5, "Track Types (0: MFT-MCH-MID, 1: none, 2: MFT-MCH, 3: MCH-MID, 4: MCH)"};
  AxisSpec axisNTracksMFT{3, -1.5, 1.5, "Number of MFT Tracks"};

  AxisSpec axisIsTrue{2, -0.5, 1.5, "Is True Candidate (0:No, 1:Yes)"};


  // Particle information
  // AxisSpec axisSource{6, -0.5, 5.5, "Source (0:Pi, 1:K, 2:HF, 3:Res, 4:Other, 5:Fake)"};
  AxisSpec axisPDGCode{1000, 0.0, 1000.0, "PDG Code"};
  AxisSpec axisVx{200, -100.0, 100.0, "Prodcution vertex vx [cm]"};
  AxisSpec axisVy{200, -100.0, 100.0, "Prodcution vertex vy [cm]"};
  AxisSpec axisVz{500, -800.0, 100.0, "Prodcution vertex vz [cm]"};
  AxisSpec axisR{200, -100.0, 100.0, "Production vertex R [cm]"};
  AxisSpec axisTrackTime{500, -1000.0, 1000.0, "Track Time [ns]"};
  AxisSpec axisNClsMFT{15, -0.5, 14.5, "Number of MFT Clusters"};

  // Mass
  AxisSpec axisMass{200, 0.0, 4.0, "Invariant Mass [GeV/c^{2}]"};

  AxisSpec axisIsCandidateTrue{2, -0.5, 1.5, "Is True? (0:Fake, 1:True)"};

  // 0 Best-Fake Second-Fake
  // 1 Best-True Second-Fake
  // 2 Best-Fake Second-True
  // 3 Best-True Second-True (empty bin)
  AxisSpec axisMatchStatus{4, -0.5, 3.5, "Match Status (0:FF, 1:TF, 2:FT, 3:TT)"};



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
    //histos.add("Pt/pTMCH_kToMatching", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    

    //histos.add("Pt/pTMFT", "MFT Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    //histos.add("Pt/pTMCH", "MCH Track pT; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    //histos.add("Pt/pT_MFTMCH_2D", "MFT vs MCH Tracks pT; MFT pT [GeV/c]; MCH pT [GeV/c]", kTH2F, {axisPtMFT, axisPtMCH});
    /*
    // kToMatching
    histos.add("kToMatching/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToMatching/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("kToMatching/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("kToMatching/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ});
    histos.add("kToMatching/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    
    // [MEMORY FIX] TH3F removed/replaced with 2D+1D
    // histos.add("kToMatching/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ});
    // histos.add("kToMatching/XY-pT", "MCH Track pT vs XY at Matching Plane; X [cm]; Y[cm]; pT [GeV/c]", kTH3F, {axisX, axisY, axisPtMCH});
    // histos.add("kToMatching/XY-p", "MCH Track p vs XY at Matching Plane; X [cm]; Y[cm]; p [GeV/c]", kTH3F, {axisX, axisY, axisP});
    // [MEMORY FIX] Alternative 2D plots for pT/p
    histos.add("kToMatching/X_pT", "X vs pT; X [cm]; pT [GeV/c]", kTH2F, {axisX, axisPtMCH});
    histos.add("kToMatching/Y_pT", "Y vs pT; Y [cm]; pT [GeV/c]", kTH2F, {axisY, axisPtMCH});


    histos.add("kToMatching/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("kToMatching/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("kToMatching/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToMatching/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // k kToVertex
    histos.add("kToVertex/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToVertex/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("kToVertex/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("kToVertex/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ_Vertex});
    histos.add("kToVertex/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    // [MEMORY FIX] Removed 3D
    // histos.add("kToVertex/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ_Vertex});
    
    histos.add("kToVertex/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("kToVertex/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("kToVertex/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToVertex/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // k kToDCA
    histos.add("kToDCA/pTMCH", "MCH Track pT at Matching Plane; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToDCA/X_MCH", "MCH Track X at Matching Plane; X [cm]; Counts", kTH1F, {axisX});
    histos.add("kToDCA/Y_MCH", "MCH Track Y at Matching Plane; Y [cm]; Counts", kTH1F, {axisY});
    histos.add("kToDCA/Z_MCH", "MCH Track Z at Matching Plane; Z [cm]; Counts", kTH1F, {axisZ_Vertex});
    histos.add("kToDCA/X-Y_MCH", "MCH Track X vs Y at Matching Plane; X [cm]; Y[cm]", kTH2F, {axisX, axisY});
    // [MEMORY FIX] Removed 3D
    // histos.add("kToDCA/XYZ", "3D Position at Matching Plane;X [cm];Y [cm];Z [cm]", kTH3F, {axisX, axisY, axisZ_Vertex});
    
    histos.add("kToDCA/pMCH", "MCH Total Momentum at Matching Plane; p [GeV/c]", kTH1F, {axisP});
    histos.add("kToDCA/EtaMCH", "MCH #eta at Matching Plane; #eta", kTH1F, {axisEta});
    histos.add("kToDCA/PhiMCH", "MCH #phi at Matching Plane; #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToDCA/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});

    // kToRabs
    histos.add("kToRabs/pTMCH", "MCH Track pT at After Absorber (z = -505 cm); pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("kToRabs/X_MCH", "MCH Track X at After Absorber (z = -505 cm); X [cm]; Counts", kTH1F, {axisXRabs});
    histos.add("kToRabs/Y_MCH", "MCH Track Y at After Absorber (z = -505 cm); Y [cm]; Counts", kTH1F, {axisYRabs});
    histos.add("kToRabs/Z_MCH", "MCH Track Z at After Absorber (z = -505 cm); Z [cm]; Counts", kTH1F, {axisZ_Absorber});
    histos.add("kToRabs/X-Y_MCH", "MCH Track X vs Y at After Absorber (z = -505 cm); X [cm]; Y[cm]", kTH2F, {axisXRabs, axisYRabs});
    // [MEMORY FIX] Removed 3D
    // histos.add("kToRabs/XYZ", "3D Position at After Absorber (z = -505 cm);X [cm];Y [cm];Z [cm]", kTH3F, {axisXRabs, axisYRabs, axisZ_Absorber});
    
    histos.add("kToRabs/pMCH", "MCH Total Momentum at After Absorber (z = -505 cm); p [GeV/c]", kTH1F, {axisP});
    histos.add("kToRabs/EtaMCH", "MCH #eta at After Absorber (z = -505 cm); #eta", kTH1F, {axisEta});
    histos.add("kToRabs/PhiMCH", "MCH #phi at After Absorber (z = -505 cm); #phi [rad]", kTH1F, {axisPhi});
    histos.add("kToRabs/Charge", "MCH Charge; Charge", kTH1F, {axisCharge});
    */

    // Check number of tracks
    histos.add("NumberOfTracks/TracksMFT", "Number of MFT Tracks;N_{tracks}", kTH1F, {axisTracksMFT});
    histos.add("NumberOfTracks/TracksMCH", "Number of MCH Tracks;N_{tracks}", kTH1F, {axisTracksMCH});
    histos.add("NumberOfTracks/TracksMFT-MCH", "2D", kTH2F, {axisTracksMFT, axisTracksMCH});

    // Number of Track types
    histos.add("Check/CountTrackTypeswithoutCut", "Track Types (0: MFT-MCH-MID, 1: none, 2: MFT-MCH, 3: MCH-MID, 4: MCH)", kTH1F, {axisTrackTypes});
    histos.add("Check/CountTrackTypeswithCut", "Track Types (0: MFT-MCH-MID, 1: none, 2: MFT-MCH, 3: MCH-MID, 4: MCH)", kTH1F, {axisTrackTypes});
    histos.add("Check/CountMFTTrackswithCut", "Number of MFT Tracks per Event", kTH1F, {axisTracksMFT});
    histos.add("Check/CountMFTTrackswithoutCut", "Number of MFT Tracks per Event", kTH1F, {axisNTracksMFT});

    // Particle information
    // histos.add("Particle/Source/SourceParticle", "Source Particle; Source (0:Pi, 1:K, 2:HF, 3:Res, 4:Other, 5:Fake); Counts", kTH1F, {axisSource});
    histos.add("Particle/Detect/Particle_PDGCode", "Detected Muon PDG Code; PDG Code; Counts", kTH1F, {axisPDGCode});

    histos.add("Particle/Detect/Particle_vx", "Detected Particle X production vertex; vx [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/Particle_vy", "Detected Particle Y production vertex; vy [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/Particle_vz", "Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Particle_vxy", "Detected Particle XY production vertex; vx [cm]; vy [cm]", kTH2F, {axisVx, axisVy});
    histos.add("Particle/Detect/TruePrimaryParticle_pT", "True Primary Detected Particle pT; pT [GeV/c]; Counts", kTH1F, {axisPtMCH});
    // [MEMORY FIX] Removed 3D
    // histos.add("Particle/Detect/Particle_vxyz", "Detected Particle XYZ production vertex; vx [cm]; vy [cm]; vz [cm]", kTH3F, {axisVx, axisVy, axisVz});
    
    histos.add("Particle/Detect/Particle_vz_vs_chi2", "Detected Particle vz vs Matching Chi2; vz [cm]; #chi^{2}", kTH2F, {axisVz, axisMatchChi2});
    histos.add("Particle/Detect/Particle_vz_vs_tracktime", "Detected Particle vz vs Track Time; vz [cm]; Track Time [ns]", kTH2F, {axisVz, axisTrackTime});

    histos.add("Particle/Detect/TrueParticles_vz", "True Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    // [MEMORY FIX] Removed 3D
    // histos.add("Particle/Detect/TrueParticle_vxyz", "True Detected Particle XYZ production vertex; vx [cm]; vy [cm]; vz [cm]", kTH3F, {axisVx, axisVy, axisVz});
    histos.add("Particle/Detect/TrueParticle_vz_vs_chi2", "True Detected Particle vz vs Matching Chi2; vz [cm]; #chi^{2}", kTH2F, {axisVz, axisMatchChi2});
    // histos.add("Particle/Detect/TrueParticle_vz_vs_tracktime", "True Detected Particle vz vs Track Time; vz [cm]; Track Time [ns]", kTH2F, {axisVz, axisTrackTime});
    histos.add("Particle/Detect/TruePrimaryParticle_vz", "True Primary Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/TrueNotPrimaryParticle_vz", "True Not Primary Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/TrueNotPrimaryParticleFromMaterial_vz", "True Not Primary Detected Particle from Material Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/TrueNotPrimaryParticleFromDecay_vz", "True Not Primary Detected Particle not from Material Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/TrueNotPrimaryParticleFromMaterial_pT", "True Not Primary Detected Particle from Material pT; pt [GeV/c]; Counts", kTH1F, {axisPtMCH});

    histos.add("Particle/Detect/FakeParticles_vz", "Fake Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    // [MEMORY FIX] Removed 3D
    // histos.add("Particle/Detect/FakeParticle_vxyz", "Fake Detected Particle XYZ production vertex; vx [cm]; vy [cm]; vz [cm]", kTH3F, {axisVx, axisVy, axisVz});
    histos.add("Particle/Detect/FakeParticle_vz_vs_chi2", "Fake Detected Particle vz vs Matching Chi2; vz [cm]; #chi^{2}", kTH2F, {axisVz, axisMatchChi2});
    // histos.add("Particle/Detect/FakeParticle_vz_vs_tracktime", "Fake Detected Particle vz vs Track Time; vz [cm]; Track Time [ns]", kTH2F, {axisVz, axisTrackTime});
    histos.add("Particle/Detect/FakePrimaryParticle_vz", "Fake Primary Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/FakeNotPrimaryParticle_vz", "Fake Not Primary Detected Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/FakeNotPrimaryParticleFromMaterial_vz", "Fake Not Primary Detected Particle from Material Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/FakeNotPrimaryParticleFromDecay_vz", "Fake Not Primary Detected Particle not from Material Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Material/Self_PDG", "Fake Not Primary Detected particle from Material PDG Code; PDG Code; Counts", kTH1F, {axisPDGCode});
    histos.add("Particle/Detect/Material/Mother_PDG", "Fake Particle from Material Mother PDG Code; PDG Code; Counts", kTH1F, {axisPDGCode});
    histos.add("Particle/Detect/Material/Mother_pT", "Fake Particle from Material Mother pT; pt [GeV/c]; Counts", kTH1F, {axisPtMCH});
    histos.add("Particle/Detect/Material/Mother_vs_Self_PDG", "Fake Particle from Material Mother vs Self PDG Code; PDG Code; Counts", kTH2F, {axisPDGCode, axisPDGCode});
    histos.add("Particle/Detect/Material/Mother_vz", "Fake Particle from Material Mother vz; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Material/Mother_vxy", "Fake Particle from Material Mother vxy; vx [cm]; vy [cm]", kTH2F, {axisVx, axisVy});
    histos.add("Particle/Detect/Material/Mother_R", "Fake Particle from Material Mother R; R [cm]; Counts", kTH1F, {axisR});
    histos.add("Particle/Detect/Material/FakeNotPrimaryParticleFromMaterial_pT", "Fake Not Primary Detected Particle from Material pT; pt [GeV/c]; Counts", kTH1F, {axisPtMCH});

    // propagate
    /*histos.add("Particle/Detect/Propagate/Particle_PropagatedToVertex_vx", "Propagated Detected Particle X production vertex to Vertex; vx [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/Propagate/Particle_PropagatedToVertex_vy", "Propagated Detected Particle Y production vertex to Vertex; vy [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/Propagate/Particle_PropagatedToVertex_vz", "Propagated Detected Particle Z production vertex to Vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Propagate/Particle_PropagatedToDCA_x", "Propagated Detected Particle X production vertex to DCA; vx [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/Propagate/Particle_PropagatedToDCA_y", "Propagated Detected Particle Y production vertex to DCA; vy [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/Propagate/Particle_PropagatedToDCA_z", "Propagated Detected Particle Z production vertex to DCA; vz [cm]; Counts", kTH1F, {axisVz});
*/
    // MFT track
    histos.add("Particle/Detect/MFT/Particle_PDGCode", "Detected MFT Particle PDG Code; PDG Code; Counts", kTH1F, {axisPDGCode});
    histos.add("Particle/Detect/MFT/Particle_vz", "Detected MFT Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});

    histos.add("Particle/Detect/MFT/PrimaryParticle_vx", "Primary Detected MFT Particle X production vertex; vx [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/MFT/NotPrimaryParticle_vx", "Not Primary Detected MFT Particle X production vertex; vx [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/MFT/PrimaryParticle_vy", "Primary Detected MFT Particle Y production vertex; vy [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/MFT/NotPrimaryParticle_vy", "Not Primary Detected MFT Particle Y production vertex; vy [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/MFT/PrimaryParticle_vz", "Primary Detected MFT Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/MFT/NotPrimaryParticle_vz", "Not Primary Detected MFT Particle Z production vertex; vz [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/MFT/PrimaryParticle_PDGCode", "Primary Detected MFT Particle PDG Code; PDG Code; Counts", kTH1F, {axisPDGCode});
    histos.add("Particle/Detect/MFT/NotPrimaryParticle_PDGCode", "Not Primary Detected MFT Particle PDG Code; PDG Code; Counts", kTH1F, {axisPDGCode});

    histos.add("Particle/Detect/MFT/NumberOfMFTClusters", "Number of MFT Clusters Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    histos.add("Particle/Detect/MFT/PrimaryNumberOfMFTClusters", "Number of MFT Clusters Primary Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    histos.add("Particle/Detect/MFT/NotPrimaryNumberOfMFTClusters", "Number of MFT Clusters Not Primary Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    
    // nClusters (I can't find any dependency)
    /*histos.add("Particle/Detect/TrueNumberOfMFTClusters", "Number of MFT Clusters Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    histos.add("Particle/Detect/TruePrimaryNumberOfMFTClusters", "Number of MFT Clusters True Primary Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    histos.add("Particle/Detect/TrueNotPrimaryNumberOfMFTClusters", "Number of MFT Clusters True Not Primary Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    histos.add("Particle/Detect/FakeNumberOfMFTClusters", "Number of MFT Clusters Fake Detected Particles; N_{clusters}; Counts", kTH1F, {axisNClsMFT});
    */

    /*
    // z
    histos.add("Particle/Detect/True/PrimaryMFTProductionPoint_z", "Primary Detected Particle MFT Production Point z; z [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/True/SecondaryMFTProductionPoint_z", "Secondary Detected Particle MFT Production Point z; z [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Fake/MFTProductionPoint_z", "Fake Detected Particle MFT Production Point z; z [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Fake/PrimaryMFTProductionPoint_z", "Primary Fake Detected Particle MFT Production Point z; z [cm]; Counts", kTH1F, {axisVz});
    histos.add("Particle/Detect/Fake/SecondaryMFTProductionPoint_z", "Secondary Fake Detected Particle MFT Production Point z; z [cm]; Counts", kTH1F, {axisVz});

    // x
    histos.add("Particle/Detect/True/PrimaryMFTProductionPoint_x", "Primary Detected Particle MFT Production Point x; x [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/True/SecondaryMFTProductionPoint_x", "Secondary Detected Particle MFT Production Point x; x [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/Fake/MFTProductionPoint_x", "Fake Detected Particle MFT Production Point x; x [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/Fake/PrimaryMFTProductionPoint_x", "Primary Fake Detected Particle MFT Production Point x; x [cm]; Counts", kTH1F, {axisVx});
    histos.add("Particle/Detect/Fake/SecondaryMFTProductionPoint_x", "Secondary Fake Detected Particle MFT Production Point x; x [cm]; Counts", kTH1F, {axisVx});

    // y
    histos.add("Particle/Detect/True/PrimaryMFTProductionPoint_y", "Primary Detected Particle MFT Production Point y; y [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/True/SecondaryMFTProductionPoint_y", "Secondary Detected Particle MFT Production Point y; y [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/Fake/MFTProductionPoint_y", "Fake Detected Particle MFT Production Point y; y [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/Fake/PrimaryMFTProductionPoint_y", "Primary Fake Detected Particle MFT Production Point y; y [cm]; Counts", kTH1F, {axisVy});
    histos.add("Particle/Detect/Fake/SecondaryMFTProductionPoint_y", "Secondary Fake Detected Particle MFT Production Point y; y [cm]; Counts", kTH1F, {axisVy});

    // p
    histos.add("Particle/Detect/True/PrimaryMFTProductionPoint_p", "Primary Detected Particle MFT Production Point p; p [GeV/c]; Counts", kTH1F, {axisP});
    histos.add("Particle/Detect/True/SecondaryMFTProductionPoint_p", "Secondary Detected Particle MFT Production Point p; p [GeV/c]; Counts", kTH1F, {axisP});
    histos.add("Particle/Detect/Fake/MFTProductionPoint_p", "Fake Detected Particle MFT Production Point p; p [GeV/c]; Counts", kTH1F, {axisP});
    histos.add("Particle/Detect/Fake/PrimaryMFTProductionPoint_p", "Primary Fake Detected Particle MFT Production Point p; p [GeV/c]; Counts", kTH1F, {axisP});
    histos.add("Particle/Detect/Fake/SecondaryMFTProductionPoint_p", "Secondary Fake Detected Particle MFT Production Point p; p [GeV/c]; Counts", kTH1F, {axisP});

    // pT
    histos.add("Particle/Detect/True/PrimaryMFTProductionPoint_pt", "Primary Detected Particle MFT Production Point pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Particle/Detect/True/SecondaryMFTProductionPoint_pt", "Secondary Detected Particle MFT Production Point pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Particle/Detect/Fake/MFTProductionPoint_pt", "Fake Detected Particle MFT Production Point pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Particle/Detect/Fake/PrimaryMFTProductionPoint_pt", "Primary Fake Detected Particle MFT Production Point pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});
    histos.add("Particle/Detect/Fake/SecondaryMFTProductionPoint_pt", "Secondary Fake Detected Particle MFT Production Point pT; pT [GeV/c]; Counts", kTH1F, {axisPtMFT});

    // eta
    histos.add("Particle/Detect/True/PrimaryMFTProductionPoint_eta", "Primary Detected Particle MFT Production Point eta; eta; Counts", kTH1F, {axisEta});
    histos.add("Particle/Detect/True/SecondaryMFTProductionPoint_eta", "Secondary Detected Particle MFT Production Point eta; eta; Counts", kTH1F, {axisEta});
    histos.add("Particle/Detect/Fake/MFTProductionPoint_eta", "Fake Detected Particle MFT Production Point eta; eta; Counts", kTH1F, {axisEta});
    histos.add("Particle/Detect/Fake/PrimaryMFTProductionPoint_eta", "Primary Fake Detected Particle MFT Production Point eta; eta; Counts", kTH1F, {axisEta});
    histos.add("Particle/Detect/Fake/SecondaryMFTProductionPoint_eta", "Secondary Fake Detected Particle MFT Production Point eta; eta; Counts", kTH1F, {axisEta});
    */

    // Ambiguity
    histos.add("Ambiguity/N_Candidates", "Number of Global Muon Candidates per MCH Track; N_{candidates}; Counts", kTH1F, {axisNTracksMFT});
    histos.add("Ambiguity/Chi2_Best", "Best Match #chi^{2}; #chi^{2}_{best}; Counts", kTH1F, {axisMatchChi2});
    histos.add("Ambiguity/Chi2_Second", "Second Best Match #chi^{2}; #chi^{2}_{2nd}; Counts", kTH1F, {axisMatchChi2});
    histos.add("Ambiguity/Chi2_Correlation", "Best vs Second Best #chi^{2}; #chi^{2}_{best}; #chi^{2}_{2nd}", kTH2F, {axisMatchChi2, axisMatchChi2});
    histos.add("Ambiguity/DeltaChi2", "Difference in #chi^{2} (2nd - Best); #Delta#chi^{2}; Counts", kTH1F, {axisMatchChi2});
    
    histos.add("Ambiguity/MC_Best_vs_Second", "MC Truth: Best vs 2nd Candidate; Best is True?; 2nd is True?", kTH2F, {axisIsCandidateTrue, axisIsCandidateTrue});
    histos.add("Ambiguity/DeltaChi2_vs_BestTrue", "#Delta#chi^{2} vs Best Candidate Truth; #Delta#chi^{2} (2nd - Best); Best is True?", kTH2F, {axisMatchChi2, axisIsCandidateTrue});
    
    // pT vs Delta Chi2
    histos.add("Ambiguity/Pt_vs_DeltaChi2", "p_{T} vs #Delta#chi^{2}; p_{T} [GeV/c]; #Delta#chi^{2} (2nd - Best)", kTH2F, {axisPtMCH, axisMatchChi2});
    histos.add("Ambiguity/Pt_vs_MatchStatus", "p_{T} vs Match Status; p_{T} [GeV/c]; Status (0:FF, 1:TF, 2:FT, 3:TT)", kTH2F, {axisPtMCH, axisMatchStatus});
    histos.add("Ambiguity/Pt_vs_BestIsTrue", "p_{T} vs Best is True; p_{T} [GeV/c]; Best is True (0/1)", kTH2F, {axisPtMCH, axisIsCandidateTrue});

    // Mass
    histos.add("Mass/Global_Mass_Unlike", "Global Dimuon Mass (Unlike Sign); M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});
    histos.add("Mass/Global_Mass_Like", "Global Dimuon Mass (Like Sign); M_{#mu#mu} [GeV/c^{2}]; Counts", kTH1F, {axisMass});

  }

  // [MEMORY FIX] Use std::vector instead of std::map for better performance and memory
  std::vector<Counts> countTracksPerCollision(aod::Collisions const& collisions,
                                              aod::FwdTracks const& mchTracks,
                                              aod::MFTTracks const& mftTracks) 
  {
    std::vector<Counts> countsVec(collisions.size()); 

    for (auto const& mchTr : mchTracks) {
      if (mchTr.trackType() == 3) {
         int collId = mchTr.collisionId();
         if(collId >= 0 && collId < countsVec.size()) countsVec[collId].mch++;
      }
    }
    for (auto const& mftTr : mftTracks) {
      int collId = mftTr.collisionId();
      if(collId >= 0 && collId < countsVec.size()) countsVec[collId].mft++;
    }
  
    for (auto const& cnt : countsVec) {
      if (cnt.mft > 0 && cnt.mch > 0) {
        histos.fill(HIST("NumberOfTracks/TracksMFT"), static_cast<float>(cnt.mft));
        histos.fill(HIST("NumberOfTracks/TracksMCH"), static_cast<float>(cnt.mch));
        histos.fill(HIST("NumberOfTracks/TracksMFT-MCH"), static_cast<float>(cnt.mft), static_cast<float>(cnt.mch));
      }
    }
    return countsVec;
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
  /*void processMatchingAnalysis_likeDQ(MCHMuons const& mchJoined,
                                MyEvents const& collisions)
  {
    for (auto const& globalTr : mchJoined) {
      // track
      if (globalTr.trackType() != 3) continue;
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

      int collId = globalTr.collisionId();

      if (collId < 0 || collId >= collisions.size()) {
        continue;
      }

      auto const& collRow = collisions.iteratorAt(collId);
      // This information has Track 0 or 1 type only (including MFT track)

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
        
        // [MEMORY FIX] Removed 3D fill
        // histos.fill(HIST("kToMatching/XYZ"), xAtMP, yAtMP, zAtMP);
        // histos.fill(HIST("kToMatching/XY-pT"), xAtMP, yAtMP, ptAtMP);
        // histos.fill(HIST("kToMatching/XY-p"), xAtMP, yAtMP, pAtMP);
        
        // [MEMORY FIX] Fill 2D projections instead
        histos.fill(HIST("kToMatching/X_pT"), xAtMP, ptAtMP);
        histos.fill(HIST("kToMatching/Y_pT"), yAtMP, ptAtMP);
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
      // [MEMORY FIX] Removed 3D fill
      // histos.fill(HIST("kToVertex/XYZ"), xAtVX, yAtVX, zAtVX);

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
      // [MEMORY FIX] Removed 3D fill
      // histos.fill(HIST("kToDCA/XYZ"), xAtDCA, yAtDCA, zAtDCA);

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
        // [MEMORY FIX] Removed 3D fill
        // histos.fill(HIST("kToRabs/XYZ"), xAtRabs, yAtRabs, zAtRabs);
       }
    }
  }*/


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
      if (track.chi2() < 0.0 || track.chi2() > 1e6) continue;
      if (track.chi2MatchMCHMID() < 0.0 || track.chi2MatchMCHMID() > 1e6) continue; 
      if (track.chi2MatchMCHMFT() < 0.0 || track.chi2MatchMCHMFT() > 1e6) continue;
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

  // Particle infformation analysis
  // process(soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks,aod::McParticles const& mcParticles)
  template <typename MCFwdTracks, typename MCMFTTracks, typename mcParticles>
  void ParticleInfoAnalysis(MCFwdTracks const& fwdTracks, MCMFTTracks const& mftTracks, mcParticles const& particles)
  {
    for (auto const& track : fwdTracks) {
      // =======================================
      // Forward Tracks Analysis
      // =======================================
      // Select track type
      if (track.trackType() != 0) continue; // MFT-MCH-MID tracks only
      // if (track.trackType() != 0) continue; // MCH-MID tracks only

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
      // 0.3 GeV/c でCut,dqの人たち

      
      int mcId = track.mcParticleId();
      auto mcParticle = particles.iteratorAt(mcId);
      int pdgCode = mcParticle.pdgCode();

      // Count number fo particle and particle names
      histos.fill(HIST("Particle/Detect/Particle_PDGCode"), pdgCode);

      // Production Point
      float vx = mcParticle.vx();
      float vy = mcParticle.vy();
      float vz = mcParticle.vz();

      float pT = track.pt();

      float chi2_matching = track.chi2MatchMCHMFT();
      float tracktime = track.trackTime();

      histos.fill(HIST("Particle/Detect/Particle_vx"), vx);
      histos.fill(HIST("Particle/Detect/Particle_vy"), vy);
      histos.fill(HIST("Particle/Detect/Particle_vz"), vz);
      histos.fill(HIST("Particle/Detect/Particle_vxy"), vx, vy);
      // [MEMORY FIX] Removed 3D fill
      // histos.fill(HIST("Particle/Detect/Particle_vxyz"), vx, vy, vz);
      histos.fill(HIST("Particle/Detect/Particle_vz_vs_chi2"), vz, chi2_matching);
      histos.fill(HIST("Particle/Detect/Particle_vz_vs_tracktime"), vz, tracktime);

      // Propagate to vertex and DCA (FwdTracks)
      /*int collId = track.collisionId();

      if (collId < 0 || collId >= collisions.size()) {
        continue;
      }

      auto const& collRow = collisions.iteratorAt(collId);
      auto propAtVertex = VarManager::PropagateMuon(track, collRow, VarManager::kToVertex);
      // auto propAtDCA = VarManager::PropagateMuon(track, collRow, VarManager::kToDCA);
      float vx_propagate = propAtVertex.getX();
      float vy_propagate = propAtVertex.getY();
      float vz_propagate = propAtVertex.getZ();
      // float xDCA = propAtDCA.getX();
      // float yDCA = propAtDCA.getY();
      // float zDCA = propAtDCA.getZ();

      histos.fill(HIST("Particle/Detect/Propagate/Particle_PropagatedToVertex_vx"), vx_propagate);
      histos.fill(HIST("Particle/Detect/Propagate/Particle_PropagatedToVertex_vy"), vy_propagate);
      histos.fill(HIST("Particle/Detect/Propagate/Particle_PropagatedToVertex_vz"), vz_propagate);
      // histos.fill(HIST("Particle/Detect/Propagate/Particle_PropagatedToDCA_x"), xDCA);
      // histos.fill(HIST("Particle/Detect/Propagate/Particle_PropagatedToDCA_y"), yDCA);
      // histos.fill(HIST("Particle/Detect/Propagate/Particle_PropagatedToDCA_z"), zDCA);
*/
      int mftTrackID = track.matchMFTTrackId();
      auto mfttrack = mftTracks.rawIteratorAt(mftTrackID);
      // int nClustersMFT = mfttrack.nClusters();
      /*float z = mfttrack.z();
      float x = mfttrack.x();
      float y = mfttrack.y();
      float p = mfttrack.p();
      float pT = mfttrack.pt();
      float eta = mfttrack.eta();*/


      int mcLabel = track.mcMask();
      if (mcLabel == 0) {
        histos.fill(HIST("Particle/Detect/TrueParticles_vz"), vz);

        // histos.fill(HIST("Particle/Detect/TrueNumberOfMFTClusters"), nClustersMFT);

        // [MEMORY FIX] Removed 3D fill
        // histos.fill(HIST("Particle/Detect/TrueParticle_vxyz"), vx, vy, vz);
        histos.fill(HIST("Particle/Detect/TrueParticle_vz_vs_chi2"), vz, chi2_matching);
        // histos.fill(HIST("Particle/Detect/TrueParticle_vz_vs_tracktime"), vz, tracktime);
        if (mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("Particle/Detect/TruePrimaryParticle_vz"), vz);
          histos.fill(HIST("Particle/Detect/TruePrimaryParticle_pT"), pT);

          // histos.fill(HIST("Particle/Detect/TruePrimaryNumberOfMFTClusters"), nClustersMFT);
          /*histos.fill(HIST("Particle/Detect/True/PrimaryMFTProductionPoint_z"), z);
          histos.fill(HIST("Particle/Detect/True/PrimaryMFTProductionPoint_x"), x);
          histos.fill(HIST("Particle/Detect/True/PrimaryMFTProductionPoint_y"), y);
          histos.fill(HIST("Particle/Detect/True/PrimaryMFTProductionPoint_p"), p);
          histos.fill(HIST("Particle/Detect/True/PrimaryMFTProductionPoint_pt"), pT);
          histos.fill(HIST("Particle/Detect/True/PrimaryMFTProductionPoint_eta"), eta);*/
        } else {
          histos.fill(HIST("Particle/Detect/TrueNotPrimaryParticle_vz"), vz);
          // histos.fill(HIST("Particle/Detect/TrueNotPrimaryNumberOfMFTClusters"), nClustersMFT);
          /*histos.fill(HIST("Particle/Detect/True/SecondaryMFTProductionPoint_z"), z);
          histos.fill(HIST("Particle/Detect/True/SecondaryMFTProductionPoint_x"), x);
          histos.fill(HIST("Particle/Detect/True/SecondaryMFTProductionPoint_y"), y);
          histos.fill(HIST("Particle/Detect/True/SecondaryMFTProductionPoint_p"), p);
          histos.fill(HIST("Particle/Detect/True/SecondaryMFTProductionPoint_pt"), pT);
          histos.fill(HIST("Particle/Detect/True/SecondaryMFTProductionPoint_eta"), eta);*/
          if (mcParticle.producedByGenerator()) {
            // For Check
            histos.fill(HIST("Particle/Detect/TrueNotPrimaryParticleFromDecay_vz"), vz); // This HIST should be empty!!!!
          } else {
            histos.fill(HIST("Particle/Detect/TrueNotPrimaryParticleFromMaterial_vz"), vz);
            histos.fill(HIST("Particle/Detect/TrueNotPrimaryParticleFromMaterial_pT"), pT);
          }
        }
      } else {
        histos.fill(HIST("Particle/Detect/FakeParticles_vz"), vz);
        // [MEMORY FIX] Removed 3D fill
        // histos.fill(HIST("Particle/Detect/FakeParticle_vxyz"), vx, vy, vz);
        histos.fill(HIST("Particle/Detect/FakeParticle_vz_vs_chi2"), vz, chi2_matching);
        // histos.fill(HIST("Particle/Detect/FakeParticle_vz_vs_tracktime"), vz, tracktime);
        // histos.fill(HIST("Particle/Detect/FakeNumberOfMFTClusters"), nClustersMFT);
        /*histos.fill(HIST("Particle/Detect/Fake/MFTProductionPoint_z"), z);
        histos.fill(HIST("Particle/Detect/Fake/MFTProductionPoint_x"), x);
        histos.fill(HIST("Particle/Detect/Fake/MFTProductionPoint_y"), y);
        histos.fill(HIST("Particle/Detect/Fake/MFTProductionPoint_p"), p);
        histos.fill(HIST("Particle/Detect/Fake/MFTProductionPoint_pt"), pT);
        histos.fill(HIST("Particle/Detect/Fake/MFTProductionPoint_eta"), eta);*/
        if (mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("Particle/Detect/FakePrimaryParticle_vz"), vz);
          /*histos.fill(HIST("Particle/Detect/Fake/PrimaryMFTProductionPoint_z"), z);
          histos.fill(HIST("Particle/Detect/Fake/PrimaryMFTProductionPoint_x"), x);
          histos.fill(HIST("Particle/Detect/Fake/PrimaryMFTProductionPoint_y"), y);
          histos.fill(HIST("Particle/Detect/Fake/PrimaryMFTProductionPoint_p"), p);
          histos.fill(HIST("Particle/Detect/Fake/PrimaryMFTProductionPoint_pt"), pT);
          histos.fill(HIST("Particle/Detect/Fake/PrimaryMFTProductionPoint_eta"), eta);*/
        } else {
          histos.fill(HIST("Particle/Detect/FakeNotPrimaryParticle_vz"), vz);
          /*histos.fill(HIST("Particle/Detect/Fake/SecondaryMFTProductionPoint_z"), z);
          histos.fill(HIST("Particle/Detect/Fake/SecondaryMFTProductionPoint_x"), x);
          histos.fill(HIST("Particle/Detect/Fake/SecondaryMFTProductionPoint_y"), y);
          histos.fill(HIST("Particle/Detect/Fake/SecondaryMFTProductionPoint_p"), p);
          histos.fill(HIST("Particle/Detect/Fake/SecondaryMFTProductionPoint_pt"), pT);
          histos.fill(HIST("Particle/Detect/Fake/SecondaryMFTProductionPoint_eta"), eta);*/
          if (mcParticle.producedByGenerator()) {
            histos.fill(HIST("Particle/Detect/FakeNotPrimaryParticleFromDecay_vz"), vz);
          } else {
            histos.fill(HIST("Particle/Detect/FakeNotPrimaryParticleFromMaterial_vz"), vz);
            histos.fill(HIST("Particle/Detect/Material/FakeNotPrimaryParticleFromMaterial_pT"), pT);
            int selfPdg = mcParticle.pdgCode();
            histos.fill(HIST("Particle/Detect/Material/Self_PDG"), selfPdg);

            auto mothers = mcParticle.template mothers_as<o2::aod::McParticles>();

            for (auto const& mother : mothers) {
                int motherPdg = mother.pdgCode();

                float motherVx = mother.vx();
                float motherVy = mother.vy();
                float motherVz = mother.vz();
                float motherpT = mother.pt();
                float motherR = std::sqrt(motherVx*motherVx + motherVy*motherVy);
                
                histos.fill(HIST("Particle/Detect/Material/Mother_PDG"), motherPdg);
                histos.fill(HIST("Particle/Detect/Material/Mother_vs_Self_PDG"), motherPdg, selfPdg);
                histos.fill(HIST("Particle/Detect/Material/Mother_pT"), motherpT);
                histos.fill(HIST("Particle/Detect/Material/Mother_vz"), motherVz);
                histos.fill(HIST("Particle/Detect/Material/Mother_vxy"), motherVx, motherVy);
                histos.fill(HIST("Particle/Detect/Material/Mother_R"), motherR);
            }
          }
          
        }
      }
    }
    // How about doing same thing to MFT tracks?
    // =======================================
    // MFT Tracks Analysis
    // =======================================
    for (auto const& track : mftTracks) {
      if (track.eta() < -3.6 || track.eta() > -2.5) continue;

      int mcId = track.mcParticleId();

      if (mcId < 0) continue;

      auto mcParticle = particles.iteratorAt(mcId);
      int pdgCode = mcParticle.pdgCode();

      histos.fill(HIST("Particle/Detect/MFT/Particle_PDGCode"), pdgCode);
      
      int nClustersMFT = track.nClusters();
      histos.fill(HIST("Particle/Detect/MFT/NumberOfMFTClusters"), nClustersMFT);

      float vx = mcParticle.vx();
      float vy = mcParticle.vy();
      float vz = mcParticle.vz();


      histos.fill(HIST("Particle/Detect/MFT/Particle_vz"), vz);
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("Particle/Detect/MFT/PrimaryParticle_vx"), vx);
        histos.fill(HIST("Particle/Detect/MFT/PrimaryParticle_vy"), vy);
        histos.fill(HIST("Particle/Detect/MFT/PrimaryParticle_vz"), vz);
        histos.fill(HIST("Particle/Detect/MFT/PrimaryParticle_PDGCode"), pdgCode);
        histos.fill(HIST("Particle/Detect/MFT/PrimaryNumberOfMFTClusters"), nClustersMFT);
      } else {
        histos.fill(HIST("Particle/Detect/MFT/NotPrimaryParticle_vx"), vx);
        histos.fill(HIST("Particle/Detect/MFT/NotPrimaryParticle_vy"), vy);
        histos.fill(HIST("Particle/Detect/MFT/NotPrimaryParticle_vz"), vz);
        histos.fill(HIST("Particle/Detect/MFT/NotPrimaryParticle_PDGCode"), pdgCode);
        histos.fill(HIST("Particle/Detect/MFT/NotPrimaryNumberOfMFTClusters"), nClustersMFT);
      }
      // To cut secondary particles from material - ncluster? 

    } 
   }



  // ====================================
  // Ambiguity Analysis
  // ====================================
  /*void processAmbiguityCheck(MCHMuons const& mchJoined)
  {
    struct Candidate {
      int index;
      float chi2;
      bool isTrue;
    };

    std::map<std::pair<int, int>, std::vector<Candidate>> candidatesMap;

    for (int i = 0; i < mchJoined.size(); ++i) {
      auto const& track = mchJoined.iteratorAt(i);

      if (track.trackType() != 0) continue; 

      if (track.eta() < -3.6 || track.eta() > -2.5) continue;
      
      int mchID = track.matchMCHTrackId(); 
      int collID = track.collisionId();
      float chi2Match = track.chi2MatchMCHMFT();

      int mcLabel = track.mcMask();
      bool isTrue = (mcLabel == 0);

      candidatesMap[{collID, mchID}].push_back({i, chi2Match, isTrue});
    }

    for (auto& [key, candidates] : candidatesMap) {
      int nCandidates = candidates.size();
      histos.fill(HIST("Ambiguity/N_Candidates"), nCandidates);

      if (nCandidates >= 2) {
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
          return a.chi2 < b.chi2;
        });

        float chi2_best = candidates[0].chi2;
        bool isTrue_best = candidates[0].isTrue;

        float chi2_second = candidates[1].chi2;
        bool isTrue_second = candidates[1].isTrue;

        float delta_chi2 = chi2_second - chi2_best;
        

        histos.fill(HIST("Ambiguity/Chi2_Best"), chi2_best);
        histos.fill(HIST("Ambiguity/Chi2_Second"), chi2_second);
        histos.fill(HIST("Ambiguity/Chi2_Correlation"), chi2_best, chi2_second);
        histos.fill(HIST("Ambiguity/DeltaChi2"), delta_chi2);

        histos.fill(HIST("Ambiguity/MC_Best_vs_Second"), 
                    isTrue_best ? 1.0 : 0.0, 
                    isTrue_second ? 1.0 : 0.0);

        histos.fill(HIST("Ambiguity/DeltaChi2_vs_BestTrue"), 
                    delta_chi2, 
                    isTrue_best ? 1.0 : 0.0);
      }
    }
  }*/

  // ====================================
  // Ambiguity Analysis
  // ====================================
  void processAmbiguityCheck(MCHMuons const& mchJoined)
  {
    struct Candidate {
      int index;
      float chi2;
      bool isTrue;
      float pt;
    };

    std::map<std::pair<int, int>, std::vector<Candidate>> candidatesMap;

    for (int i = 0; i < mchJoined.size(); ++i) {
      auto const& track = mchJoined.iteratorAt(i);

      if (track.trackType() != 0) continue; 
      // Acceptance cuts
      if (track.eta() < -3.6 || track.eta() > -2.5) continue;
      
      int mchID = track.matchMCHTrackId(); 
      int collID = track.collisionId();
      float chi2Match = track.chi2MatchMCHMFT();
      
      float pt = track.pt(); 

      int mcLabel = track.mcMask();
      bool isTrue = (mcLabel == 0);

      candidatesMap[{collID, mchID}].push_back({i, chi2Match, isTrue, pt});
    }

    for (auto& [key, candidates] : candidatesMap) {
      int nCandidates = candidates.size();
      histos.fill(HIST("Ambiguity/N_Candidates"), nCandidates);

      if (nCandidates >= 2) {
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
          return a.chi2 < b.chi2;
        });

        float chi2_best = candidates[0].chi2;
        bool isTrue_best = candidates[0].isTrue;
        float pt_best = candidates[0].pt;

        float chi2_second = candidates[1].chi2;
        bool isTrue_second = candidates[1].isTrue;
        // float pt_second = candidates[1].pt;

        float delta_chi2 = chi2_second - chi2_best;

        histos.fill(HIST("Ambiguity/Chi2_Best"), chi2_best);
        histos.fill(HIST("Ambiguity/DeltaChi2"), delta_chi2);

        histos.fill(HIST("Ambiguity/Pt_vs_DeltaChi2"), pt_best, delta_chi2);

        // 0: F-F, 1: T-F(Best=True), 2: F-T(Best=Fake), 3: T-T
        int status = 0;
        if (isTrue_best && !isTrue_second) status = 1;      
        else if (!isTrue_best && isTrue_second) status = 2;
        else if (isTrue_best && isTrue_second) status = 3;
        else status = 0;

        histos.fill(HIST("Ambiguity/Pt_vs_MatchStatus"), pt_best, static_cast<float>(status));

        histos.fill(HIST("Ambiguity/Pt_vs_BestIsTrue"), pt_best, isTrue_best ? 1.0 : 0.0);
      }
    }
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
      

      collToTrackMap[track.collisionId()].push_back(i);
    }

    const double mMu = 0.105658; // [GeV/c2]

    for (auto const& [collId, trackIndices] : collToTrackMap) {
      if (trackIndices.size() < 2) continue;

      for (size_t i = 0; i < trackIndices.size(); ++i) {
        for (size_t j = i + 1; j < trackIndices.size(); ++j) {
          
          auto const& tr1 = mchJoined.iteratorAt(trackIndices[i]);
          auto const& tr2 = mchJoined.iteratorAt(trackIndices[j]);

          float pt1 = 0, pt2 = 0;
          if (std::abs(tr1.signed1Pt()) > 1e-8) pt1 = 1.0f / std::abs(tr1.signed1Pt());
          if (std::abs(tr2.signed1Pt()) > 1e-8) pt2 = 1.0f / std::abs(tr2.signed1Pt());

          TLorentzVector v1, v2, vPair;
          v1.SetPtEtaPhiM(pt1, tr1.eta(), tr1.phi(), mMu);
          v2.SetPtEtaPhiM(pt2, tr2.eta(), tr2.phi(), mMu);
          
          vPair = v1 + v2;
          float mass = vPair.M();

          // Unlike Sign +/-, Like Sign ++/--
          if (tr1.sign() * tr2.sign() < 0) {
              histos.fill(HIST("Mass/Global_Mass_Unlike"), mass);
          } else {
              histos.fill(HIST("Mass/Global_Mass_Like"), mass);
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

  
   auto countsVec = countTracksPerCollision(rawcollisions, mchTracks, mftTracks);

   TrackTypeCounter(mchJoined);
   MFTTrackCounter(mftTracks);
   ParticleInfoAnalysis(mchJoined, mcmfttracks, mcParticles); //mctracks


   // processAmbiguityCheck(mchJoined); // ambiguity analysis: (Type0 = x5 file)

   processSimpleDimuons(mchJoined); // mass analysis


   // processMatchingAnalysis_likeDQ(mchJoined, collisions);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcParticleAnalysis>(cfg)
  };
}