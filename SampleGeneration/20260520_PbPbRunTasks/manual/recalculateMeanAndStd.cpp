#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

#include "utility.h"

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  std::string InputFile = CL.Get("Input", "output/HFDistributions_First_100.root");
  std::string OutputFile = CL.Get("Output", "output/HFDistributions_First_100_recalculated.root");

  TFile Input(InputFile.c_str(), "READ");
  if (Input.IsOpen() == false) {
    std::cerr << "Failed to open input file: " << InputFile << std::endl;
    return -1;
  }

  TFile Output(OutputFile.c_str(), "RECREATE");
  if (Output.IsOpen() == false) {
    std::cerr << "Failed to open output file: " << OutputFile << std::endl;
    return -1;
  }

  std::vector<TH2D *> negEnergyHists = loadEnergyDistributionHistograms(Input, "neg", false);
  std::vector<TH2D *> posEnergyHists = loadEnergyDistributionHistograms(Input, "pos", true);

  std::vector<std::vector<TH2D *>> H36Negative_meta;
  std::vector<std::vector<TH2D *>> H36Positive_meta;
  std::vector<std::vector<TH2D *>> H18NegativeMain_meta;
  std::vector<std::vector<TH2D *>> H18NegativeWrapLow_meta;
  std::vector<std::vector<TH2D *>> H18NegativeWrapHigh_meta;
  std::vector<std::vector<TH2D *>> H18PositiveMain_meta;
  std::vector<std::vector<TH2D *>> H18PositiveWrapLow_meta;
  std::vector<std::vector<TH2D *>> H18PositiveWrapHigh_meta;

  for (int iMode = 0; iMode < modeLabels.size(); iMode++) {
    H36Negative_meta.push_back(std::vector<TH2D *>());
    H36Positive_meta.push_back(std::vector<TH2D *>());
    H18NegativeMain_meta.push_back(std::vector<TH2D *>());
    H18NegativeWrapLow_meta.push_back(std::vector<TH2D *>());
    H18NegativeWrapHigh_meta.push_back(std::vector<TH2D *>());
    H18PositiveMain_meta.push_back(std::vector<TH2D *>());
    H18PositiveWrapLow_meta.push_back(std::vector<TH2D *>());
    H18PositiveWrapHigh_meta.push_back(std::vector<TH2D *>());
    for (int iID = 0; iID < idLabels.size(); iID++) {
      H36Negative_meta[iMode].push_back(
          new TH2D(("h36Negative_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "36-phi negative-eta occupancy;pf#eta;pf#phi", HFPFBinEdges::EtaEdges36Negative().size() - 1,
                   HFPFBinEdges::EtaEdges36Negative().data(), HFPFBinEdges::PhiEdges36().size() - 1,
                   HFPFBinEdges::PhiEdges36().data()));
      H36Positive_meta[iMode].push_back(
          new TH2D(("h36Positive_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "36-phi positive-eta occupancy;pf#eta;pf#phi", HFPFBinEdges::EtaEdges36Positive().size() - 1,
                   HFPFBinEdges::EtaEdges36Positive().data(), HFPFBinEdges::PhiEdges36().size() - 1,
                   HFPFBinEdges::PhiEdges36().data()));
      H18NegativeMain_meta[iMode].push_back(
          new TH2D(("h18NegativeMain_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "18-phi negative-eta occupancy main;pf#eta;pf#phi", HFPFBinEdges::EtaEdges18Negative().size() - 1,
                   HFPFBinEdges::EtaEdges18Negative().data(), HFPFBinEdges::PhiEdges18Main().size() - 1,
                   HFPFBinEdges::PhiEdges18Main().data()));
      H18NegativeWrapLow_meta[iMode].push_back(
          new TH2D(("h18NegativeWrapLow_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "18-phi negative-eta occupancy wrap low;pf#eta;pf#phi",
                   HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
                   HFPFBinEdges::PhiEdges18WrapLow().size() - 1, HFPFBinEdges::PhiEdges18WrapLow().data()));
      H18NegativeWrapHigh_meta[iMode].push_back(
          new TH2D(("h18NegativeWrapHigh_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "18-phi negative-eta occupancy wrap high;pf#eta;pf#phi",
                   HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
                   HFPFBinEdges::PhiEdges18WrapHigh().size() - 1, HFPFBinEdges::PhiEdges18WrapHigh().data()));
      H18PositiveMain_meta[iMode].push_back(
          new TH2D(("h18PositiveMain_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "18-phi positive-eta occupancy main;pf#eta;pf#phi", HFPFBinEdges::EtaEdges18Positive().size() - 1,
                   HFPFBinEdges::EtaEdges18Positive().data(), HFPFBinEdges::PhiEdges18Main().size() - 1,
                   HFPFBinEdges::PhiEdges18Main().data()));
      H18PositiveWrapLow_meta[iMode].push_back(
          new TH2D(("h18PositiveWrapLow_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "18-phi positive-eta occupancy wrap low;pf#eta;pf#phi",
                   HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
                   HFPFBinEdges::PhiEdges18WrapLow().size() - 1, HFPFBinEdges::PhiEdges18WrapLow().data()));
      H18PositiveWrapHigh_meta[iMode].push_back(
          new TH2D(("h18PositiveWrapHigh_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str(),
                   "18-phi positive-eta occupancy wrap high;pf#eta;pf#phi",
                   HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
                   HFPFBinEdges::PhiEdges18WrapHigh().size() - 1, HFPFBinEdges::PhiEdges18WrapHigh().data()));

      fillAvgHistograms(posEnergyHists, H36Positive_meta[iMode][iID], H18PositiveMain_meta[iMode][iID],
                        H18PositiveWrapLow_meta[iMode][iID], H18PositiveWrapHigh_meta[iMode][iID], iMode, iID);
      fillAvgHistograms(negEnergyHists, H36Negative_meta[iMode][iID], H18NegativeMain_meta[iMode][iID],
                        H18NegativeWrapLow_meta[iMode][iID], H18NegativeWrapHigh_meta[iMode][iID], iMode, iID);
    }
  }

  Output.cd();
  for (TH2D *hist : posEnergyHists)
    hist->Write();
  for (TH2D *hist : negEnergyHists)
    hist->Write();
  for (int iMode = 0; iMode < modeLabels.size(); iMode++) {
    for (int iID = 0; iID < idLabels.size(); iID++) {
      H36Negative_meta[iMode][iID]->Write();
      H36Positive_meta[iMode][iID]->Write();
      H18NegativeMain_meta[iMode][iID]->Write();
      H18NegativeWrapLow_meta[iMode][iID]->Write();
      H18NegativeWrapHigh_meta[iMode][iID]->Write();
      H18PositiveMain_meta[iMode][iID]->Write();
      H18PositiveWrapLow_meta[iMode][iID]->Write();
      H18PositiveWrapHigh_meta[iMode][iID]->Write();
    }
  }

  return 0;
}
