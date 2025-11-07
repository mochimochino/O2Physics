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

bool matchedQualityCuts(const o2::aod::FwdTrack& track)
{
    // Example matched quality criteria
    if (track.chi2() > 3.0) {
        return false; // Match Chi2 cut
    }

    // Add more criteria as needed

    return true; // Track passes all matched quality cuts
}
   // cut->AddCut(VarManager::kEta, -4.0, -2.5);
    //cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    //cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    //cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    //cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    //cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    //cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    //return cut;