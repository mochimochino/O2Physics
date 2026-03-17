// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Custom task to compute number of events per run for Luminosity calculation.
/// \author Takuma
/// \date 2026

// O2 headers
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

// ROOT headers
#include "TH1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct MyUPCLuminosityTask {

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // We loop over UDCollisions to count events per RunNumber based on Triggers
  using UDCollisionsFwd = o2::aod::UDCollisions;

  void init(InitContext const&)
  {
    // Define an axis for RunNumbers. Assuming a range that covers typical Run3 PP or PbPb runs.
    const AxisSpec axisRun{50000, 500000.5, 550000.5, "Run Number"};

    // Create histograms for different trigger selections
    registry.add("hEvents_All", "All Skimmed Events; Run Number; Counts", kTH1D, {axisRun});
  }

  //____________________________________________________________________________________________
  void process(UDCollisionsFwd const& eventCandidates)
  {
    for (const auto& col : eventCandidates) {
      int runNumber = col.runNumber();

      // Count ALL Collisions in the skimmed AO2D
      registry.fill(HIST("hEvents_All"), runNumber);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyUPCLuminosityTask>(cfgc, TaskName{"my-upc-luminosity"})};
}
