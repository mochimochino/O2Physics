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
/// \brief Simple task to read ZDC (Zero Degree Calorimeter) information
///        from UDZdcsReduced and fill histograms for EnergyCommonZNA and EnergyCommonZNC.
///        These energies indicate electromagnetic dissociation of the nuclei in UPC events.
/// \author Takuma
/// \date 2026

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// O2Physics headers
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

struct MyZDCTask {

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // ZNA energy axis: up to 200 TeV (typical scale for Pb-Pb ZDC); adjust if needed
    const AxisSpec axisZNEnergy{500, 0.0, 250.0, "E_{ZN} [TeV]"};
    const AxisSpec axisZNTime{200, -10.0, 10.0, "t_{ZN} [ns]"};

    // 1D histograms
    registry.add("hEnergyZNA", "ZNA Common Energy;E_{ZNA} [TeV];Counts", kTH1F, {axisZNEnergy});
    registry.add("hEnergyZNC", "ZNC Common Energy;E_{ZNC} [TeV];Counts", kTH1F, {axisZNEnergy});
    registry.add("hTimeZNA",   "ZNA Time;t_{ZNA} [ns];Counts",           kTH1F, {axisZNTime});
    registry.add("hTimeZNC",   "ZNC Time;t_{ZNC} [ns];Counts",           kTH1F, {axisZNTime});

    // 2D correlation between ZNA and ZNC energy (useful for topology classification)
    registry.add("hEnergyZNAvsZNC", "ZNA vs ZNC Energy;E_{ZNC} [TeV];E_{ZNA} [TeV]",
                 kTH2F, {axisZNEnergy, axisZNEnergy});

    // Topology counter: events classified by which side had neutron(s)
    // 0 = 0n0n, 1 = Xn0n (ZNA active), 2 = 0nXn (ZNC active), 3 = XnXn (both active)
    auto hTopology = registry.add<TH1>("hTopology", "Neutron topology;;Counts", kTH1I, {{4, -0.5, 3.5}});
    hTopology->GetXaxis()->SetBinLabel(1, "0n0n");
    hTopology->GetXaxis()->SetBinLabel(2, "Xn0n");
    hTopology->GetXaxis()->SetBinLabel(3, "0nXn");
    hTopology->GetXaxis()->SetBinLabel(4, "XnXn");
  }

  // Energy threshold to classify a ZDC side as "active" (neutron detected)
  // The default is 1.0 TeV - adjust based on your specific analysis
  Configurable<float> cutZNAEnergy{"cutZNAEnergy", 1.0f, "ZNA energy threshold [TeV] for neutron detection"};
  Configurable<float> cutZNCEnergy{"cutZNCEnergy", 1.0f, "ZNC energy threshold [TeV] for neutron detection"};

  void process(aod::UDZdcsReduced const& zdcs)
  {
    for (const auto& zdc : zdcs) {
      float eZNA = zdc.energyCommonZNA();
      float eZNC = zdc.energyCommonZNC();
      float tZNA = zdc.timeZNA();
      float tZNC = zdc.timeZNC();

      registry.fill(HIST("hEnergyZNA"), eZNA);
      registry.fill(HIST("hEnergyZNC"), eZNC);
      registry.fill(HIST("hTimeZNA"),   tZNA);
      registry.fill(HIST("hTimeZNC"),   tZNC);
      registry.fill(HIST("hEnergyZNAvsZNC"), eZNC, eZNA);

      // Neutron topology classification
      bool hasZNA = eZNA > static_cast<float>(cutZNAEnergy);
      bool hasZNC = eZNC > static_cast<float>(cutZNCEnergy);

      if      (!hasZNA && !hasZNC) registry.fill(HIST("hTopology"), 0); // 0n0n
      else if ( hasZNA && !hasZNC) registry.fill(HIST("hTopology"), 1); // Xn0n
      else if (!hasZNA &&  hasZNC) registry.fill(HIST("hTopology"), 2); // 0nXn
      else                         registry.fill(HIST("hTopology"), 3); // XnXn
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MyZDCTask>(cfgc, TaskName{"my-zdc-01"})};
}
