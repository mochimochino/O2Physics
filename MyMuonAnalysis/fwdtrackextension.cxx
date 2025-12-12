
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/DataTypes.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct FwdTrackExtension {

    Produces<aod::FwdTracksDCA> extendedTrackQuantities;

    HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
    AxisSpec axisDCAX{200, 0, 0.1, "X position [cm]"};
    AxisSpec axisDCAY{200, 0, 0.1, "Y position [cm]"};


   void init(InitContext const&) {
    histos.add("DCAx", "DCA x position", kTH1F, {axisDCAX});
    histos.add("DCAy", "DCA y position", kTH1F, {axisDCAY});
    histos.add("DCAxDCAy2D", "DCA X-Y", kTH2F, {axisDCAX, axisDCAY});
   }

    void process (aod::FwdTracks const& tracks, aod::Collisions const&) {
        for (auto& track : tracks) {
            
      float dcaX = -999;
      float dcaY = -999;
      if (track.has_collision()) {
        if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack || track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack || track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {

          auto const& collision = track.collision();
          double chi2 = track.chi2();
          SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
          std::vector<double> v1;
          SMatrix55 tcovs(v1.begin(), v1.end());
          o2::track::TrackParCovFwd pars1{track.z(), tpars, tcovs, chi2};
          pars1.propagateToZlinear(collision.posZ());

          dcaX = (pars1.getX() - collision.posX());
          dcaY = (pars1.getY() - collision.posY());

          histos.fill(HIST("DCAx"), dcaX);
          histos.fill(HIST("DCAy"), dcaY);
          histos.fill(HIST("DCAxDCAy2D"), dcaX, dcaY);
        }
      }
      extendedTrackQuantities(dcaX, dcaY);
        
    }
}
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FwdTrackExtension>(cfgc)};
  return workflow;
}