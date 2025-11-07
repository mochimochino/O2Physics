#include "MymuonCut.h"
#include <TMath.h>

bool isGoodMuonTrack(const o2::aod::FwdTrack& track)
{
    // Example cut criteria
    if (track.pt() < 1.0) {
        return false; // pT cut
    }
    if (track.eta() < -4.0 || track.eta() > -2.5) {
        return false; // Eta cut
    }
    if (track.chi2() > 5.0) {
        return false; // Chi2 cut
    }
    if (track.trackType() != 0 && track.trackType() != 3) {
        return false; // Type cut
    }
    // Add more criteria as needed

    return true; // Track passes all cuts
}