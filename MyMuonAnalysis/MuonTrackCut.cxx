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

    return true; // Track passes all cuts
}

bool matchedQualityCuts(const o2::aod::FwdTrack& track) //from PWGDQ/Core/CutsLibrary.cxx
{
    if (track.eta() < -4.0 || track.eta() > -2.5) {
        return false; // Match eta cut
    }
    if (track.rAtAbsorberEnd() < 17.6 || track.rAtAbsorberEnd() > 89.5) {
        return false; // Match rAtAbsorberEnd cut
    }
    if (track.pDca() > 594.0 && track.rAtAbsorberEnd() < 26.5) {
        return false; // Match pDca cut
    }
    if (track.pDca() > 324.0 && track.rAtAbsorberEnd() >= 26.5) {
        return false; // Match pDca cut
    }
    if (track.chi2() > 1e6) {
        return false; // Match chi2 cut
    }
    if (track.chi2MatchMCHMID() > 1e6) {
        return false; // Match chi2MatchMCHMID cut
    }
    if (track.chi2MatchMCHMFT() > 1e6) {
        return false; // Match chi2MatchMCHMFT cut
    }
    return true;
}
   // cut->AddCut(VarManager::kEta, -4.0, -2.5);
    //cut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    //cut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    //cut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    //cut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    //cut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    //cut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH
    //return cut;