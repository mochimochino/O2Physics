// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

/// \file   globalMuonQA.cxx
/// \brief  QA task for Global Muon (MCH-MFT) matching evaluation in MC
/// \author Takuma

#include "PWGUD/DataModel/UDTables.h"

#include <Common/Core/fwdtrackUtilities.h>

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>
#include <SimulationDataFormat/MCCompLabel.h>
#include <SimulationDataFormat/MCTruthContainer.h>

#include <TGeoGlobalMagField.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h> // TMath::Pi() のために追加

#include <array>
#include <cmath>
#include <map>

using namespace o2::framework;

struct GlobalMuonQAPars {
  bool enableMFT{true};
  float maxChi2MatchMCHMFT{100.f};
  float maxDCAxy{999.f};
  float zShift{0.0f};
};

struct GlobalMuonQA {
  Configurable<bool> cfgEnableMFT{"QA_enableMFT", true, "Enable MFT processing"};
  Configurable<float> cfgMaxChi2{"QA_maxChi2", 100.f, "Max chi2 for MCH-MFT matching"};
  Configurable<float> cfgMaxDCAxy{"QA_maxDCAxy", 999.f, "Max DCAxy"};
  Configurable<float> cfgZShift{"QA_zShift", 0.0f, "z-shift"};

  GlobalMuonQAPars fPars;
  HistogramRegistry histRegistry{"QAHistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  int fRun{0};
  float fBz{0};
  static constexpr double fCcenterMFT[3] = {0, 0, -61.4};
  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  void init(InitContext&)
  {
    fPars.enableMFT = (bool)cfgEnableMFT;
    fPars.maxChi2MatchMCHMFT = (float)cfgMaxChi2;
    fPars.maxDCAxy = (float)cfgMaxDCAxy;
    fPars.zShift = (float)cfgZShift;

    fCCDB->setURL("http://alice-ccdb.cern.ch");
    fCCDB->setCaching(true);
    fCCDBApi.init("http://alice-ccdb.cern.ch");

    // ==========================================
    // 既存の QAヒストグラム
    // ==========================================
    const AxisSpec axisChi2{200, 0., 100., "chi2 MCH-MFT"};
    histRegistry.add("hChi2_All", "chi2 MCH-MFT (All)", kTH1F, {axisChi2});
    histRegistry.add("hChi2_True", "chi2 MCH-MFT (True Match)", kTH1F, {axisChi2});
    histRegistry.add("hChi2_Fake", "chi2 MCH-MFT (Fake Match)", kTH1F, {axisChi2});

    const AxisSpec axisPt{100, 0., 20., "p_{T} (GeV/c)"};
    histRegistry.add("hChi2_vs_Pt_All", "chi2 vs pT (All)", kTH2F, {axisPt, axisChi2});
    histRegistry.add("hChi2_vs_Pt_True", "chi2 vs pT (True Match)", kTH2F, {axisPt, axisChi2});
    histRegistry.add("hChi2_vs_Pt_Fake", "chi2 vs pT (Fake Match)", kTH2F, {axisPt, axisChi2});

    const AxisSpec axisEta{100, -4.0, -2.0, "eta"};
    histRegistry.add("hChi2_vs_Eta_All", "chi2 vs eta (All)", kTH2F, {axisEta, axisChi2});

    const AxisSpec axisDCAxy{200, 0., 10., "DCA_{xy} (cm)"};
    const AxisSpec axisDCAz{200, -10., 10., "DCA_{z} (cm)"};
    histRegistry.add("hDCAxy_True", "DCAxy (True Match)", kTH1F, {axisDCAxy});
    histRegistry.add("hDCAxy_Fake", "DCAxy (Fake Match)", kTH1F, {axisDCAxy});
    histRegistry.add("hDCAz_True", "DCAz (True Match)", kTH1F, {axisDCAz});
    histRegistry.add("hDCAz_Fake", "DCAz (Fake Match)", kTH1F, {axisDCAz});

    // ==========================================
    // 新規追加: レゾリューション（分解能）ヒストグラム
    // ==========================================
    const AxisSpec axisPtRes{200, -0.5, 0.5, "(p_{T}^{reco} - p_{T}^{MC}) / p_{T}^{MC}"};
    histRegistry.add("hPtRes_True", "p_{T} Resolution (True Match)", kTH2F, {axisPt, axisPtRes});
    histRegistry.add("hPtRes_Fake", "p_{T} Resolution (Fake Match)", kTH2F, {axisPt, axisPtRes});

    const AxisSpec axisEtaRes{200, -0.1, 0.1, "#eta^{reco} - #eta^{MC}"};
    histRegistry.add("hEtaRes_True", "#eta Resolution (True Match)", kTH2F, {axisEta, axisEtaRes});
    histRegistry.add("hEtaRes_Fake", "#eta Resolution (Fake Match)", kTH2F, {axisEta, axisEtaRes});

    const AxisSpec axisPhi{100, 0., 2 * TMath::Pi(), "#phi"};
    const AxisSpec axisPhiRes{200, -0.2, 0.2, "#phi^{reco} - #phi^{MC} (rad)"};
    histRegistry.add("hPhiRes_True", "#phi Resolution (True Match)", kTH2F, {axisPhi, axisPhiRes});
    histRegistry.add("hPhiRes_Fake", "#phi Resolution (Fake Match)", kTH2F, {axisPhi, axisPhiRes});
  }

  bool isTrueGlobalMatch(int64_t trackId, const o2::aod::McFwdTrackLabels& mcLabels) const
  {
    const auto& label = mcLabels.iteratorAt(trackId);
    return label.mcMask() == 0; // 0 が True, フラグ立ち(128等) が Fake
  }

  void fillQAHistograms(const ForwardTracks::iterator& track, bool isTrueMatch, double dcaXY, double dcaZ)
  {
    float pt = track.pt();
    float eta = track.eta();
    float chi2MCHMFT = track.chi2MatchMCHMFT();

    histRegistry.fill(HIST("hChi2_All"), chi2MCHMFT);
    histRegistry.fill(HIST("hChi2_vs_Pt_All"), pt, chi2MCHMFT);
    histRegistry.fill(HIST("hChi2_vs_Eta_All"), eta, chi2MCHMFT);

    if (isTrueMatch) {
      histRegistry.fill(HIST("hChi2_True"), chi2MCHMFT);
      histRegistry.fill(HIST("hChi2_vs_Pt_True"), pt, chi2MCHMFT);
      histRegistry.fill(HIST("hDCAxy_True"), dcaXY);
      histRegistry.fill(HIST("hDCAz_True"), dcaZ);
    } else {
      histRegistry.fill(HIST("hChi2_Fake"), chi2MCHMFT);
      histRegistry.fill(HIST("hChi2_vs_Pt_Fake"), pt, chi2MCHMFT);
      histRegistry.fill(HIST("hDCAxy_Fake"), dcaXY);
      histRegistry.fill(HIST("hDCAz_Fake"), dcaZ);
    }
  }

  // MC情報を引数として追加（o2::aod::McParticles）
  void process(ForwardTracks const& fwdTracks,
               o2::aod::BCs const& bcs,
               o2::aod::Collisions const& collisions,
               o2::aod::McFwdTrackLabels const& mcFwdTrackLabels,
               o2::aod::McParticles const& mcParticles)
  {
    if (bcs.size() == 0)
      return;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != fRun && runNumber > 0) {
      fRun = runNumber;
      std::map<std::string, std::string> metadata;
      auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, fRun);
      auto ts = soreor.second;
      auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", metadata, ts);
      o2::base::Propagator::initFieldFromGRP(grpmag);
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      if (field)
        fBz = field->getBz(fCcenterMFT);
    }

    std::map<int64_t, std::array<double, 3>> colMap;
    for (const auto& col : collisions) {
      colMap[col.globalIndex()] = {col.posX(), col.posY(), col.posZ()};
    }

    for (const auto& track : fwdTracks) {
      auto trackType = track.trackType();

      if (trackType != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack &&
          trackType != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack) {
        continue;
      }

      if (track.chi2MatchMCHMFT() < 0 || track.chi2MatchMCHMFT() > fPars.maxChi2MatchMCHMFT) {
        continue;
      }

      auto trackId = track.globalIndex();
      bool isTrueMatch = isTrueGlobalMatch(trackId, mcFwdTrackLabels);

      // ==========================================
      // 分解能の計算とFill
      // ==========================================
      const auto& label = mcFwdTrackLabels.iteratorAt(trackId);
      int mcPartId = label.mcParticleId();

      // 正しいMC Particle IDを持っている場合のみ分解能を計算
      if (mcPartId >= 0 && mcPartId < mcParticles.size()) {
        auto mcPart = mcParticles.iteratorAt(mcPartId);
        float ptMC = mcPart.pt();
        float etaMC = mcPart.eta();
        float phiMC = mcPart.phi();

        float ptReco = track.pt();
        float etaReco = track.eta();
        float phiReco = track.phi();

        if (ptMC > 0.f) {
          float ptRes = (ptReco - ptMC) / ptMC; // 相対分解能
          float etaRes = etaReco - etaMC;       // 絶対分解能
          float phiRes = phiReco - phiMC;       // 絶対分解能

          // phiRes を -pi から pi の範囲に補正
          while (phiRes > TMath::Pi())
            phiRes -= 2 * TMath::Pi();
          while (phiRes < -TMath::Pi())
            phiRes += 2 * TMath::Pi();

          if (isTrueMatch) {
            histRegistry.fill(HIST("hPtRes_True"), ptMC, ptRes);
            histRegistry.fill(HIST("hEtaRes_True"), etaMC, etaRes);
            histRegistry.fill(HIST("hPhiRes_True"), phiMC, phiRes);
          } else {
            histRegistry.fill(HIST("hPtRes_Fake"), ptMC, ptRes);
            histRegistry.fill(HIST("hEtaRes_Fake"), etaMC, etaRes);
            histRegistry.fill(HIST("hPhiRes_Fake"), phiMC, phiRes);
          }
        }
      }

      // DCA計算
      double colX = 0., colY = 0., colZ = 0.;
      int64_t colId = track.collisionId();
      if (colMap.find(colId) != colMap.end()) {
        colX = colMap[colId][0];
        colY = colMap[colId][1];
        colZ = colMap[colId][2];
      }

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(track, fPars.zShift);
      std::array<double, 3> dcaOrig;
      trackPar.propagateToDCAhelix(fBz, {colX, colY, colZ}, dcaOrig);

      double dcaXY = std::sqrt(dcaOrig[0] * dcaOrig[0] + dcaOrig[1] * dcaOrig[1]);
      double dcaZ = dcaOrig[2];

      fillQAHistograms(track, isTrueMatch, dcaXY, dcaZ);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GlobalMuonQA>(cfgc)};
}
