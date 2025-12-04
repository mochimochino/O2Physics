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

struct collision_associate_central {
    HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    Configurable<int> nBinsNCTrk{"nBinsNCTrk", 5, "N bins in N FwdTracks histo"};
    Configurable<float> minNCTrk{"minNCTrk", -0.5, "min N FwdTracks"};
    Configurable<float> maxNCTrk{"maxNCTrk", 5.0, "max N FwdTracks"};

    void init(InitContext const&)
    {
        //AxisSpec axisNTrk{nBinsNTrk, minNTrk, maxNTrk, "Number of FwdTracks per Collision"};
        AxisSpec axisNCTrk{nBinsNCTrk, minNCTrk, maxNCTrk, "Number of Tracks per Collision"};
        
        //histos.add("hNfwdTracksPerColl", 
        //           "Forward Track Multiplicity per Collision; N_{FwdTracks / Collision}; Counts", 
         //          kTH1I,
          //         {axisNTrk});
        histos.add("hNTracksPerColl", 
                   "Track Multiplicity per Collision; N_{Tracks / Collision}; Counts", 
                   kTH1I,
                   {axisNCTrk});
    } 

    void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
    {
        std::map<int, int> trackCounts;
        for (auto& coll : collisions) {
     
            trackCounts[coll.globalIndex()] = 0;
        }

        for (auto& track : tracks) {
            
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
    adaptAnalysisTask<collision_associate_central>(cfg)};
}