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

using namespace o2;
using namespace o2::framework;

struct collision_associate {
    HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    Configurable<int> nBinsNTrk{"nBinsNTrk", 201, "N bins in N FwdTracks histo"};
    Configurable<float> minNTrk{"minNTrk", -0.5, "min N FwdTracks"};
    Configurable<float> maxNTrk{"maxNTrk", 200.5, "max N FwdTracks"};

    void init(InitContext const&)
    {
        AxisSpec axisNTrk{nBinsNTrk, minNTrk, maxNTrk, "Number of FwdTracks per Collision"};
        
        histos.add("hNfwdTracksPerColl", 
                   "Forward Track Multiplicity per Collision; N_{FwdTracks / Collision}; Counts", 
                   kTH1I,
                   {axisNTrk});
    } 

    void process(aod::Collisions const& collisions, aod::FwdTracks const& fwdTracks)
    {
        std::map<int, int> trackCounts;
。
        for (auto& coll : collisions) {
     
            trackCounts[coll.globalIndex()] = 0;
        }

        for (auto& track : fwdTracks) {
            
            int collId = track.collisionId(); 

            auto it = trackCounts.find(collId);
            if (it != trackCounts.end()) {
                it->second++;
            }
        }

        for (auto const& [collId, count] : trackCounts) {
            histos.fill(HIST("hNfwdTracksPerColl"), (float)count);
        }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<collision_associate>(cfg)};
}