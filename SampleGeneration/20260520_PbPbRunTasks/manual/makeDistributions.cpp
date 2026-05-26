#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

#include "Messenger.h"
#include "myhelper.h"
#include "utility.h"

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  std::string InputFileName = CL.Get("Input", "test.root");
  std::string OutputRootName = CL.Get("OutputRoot", "hf_pf_binning_check.root");
  const bool do18BinMerging = CL.GetBool("do18BinMerging", true);
  const int triggerChoice = CL.GetInt("TriggerChoice", -1); // 0: no trigger, 1: isNotBptxOR
  if (triggerChoice != 0 && triggerChoice != 1) {
    std::cerr << "Error: Invalid TriggerChoice parameter value. Expected 0 or 1." << std::endl;
    return -1;
  }
  std::cout << "Selected Trigger Choice: " << (triggerChoice == 0 ? "No Trigger" : "HLT_HIL1NotBptxOR_v") << std::endl;

  TFile *InputFile = TFile::Open(InputFileName.c_str());
  if (!InputFile || InputFile->IsZombie()) {
    std::cerr << "Error: Could not open input file " << InputFileName << std::endl;
    return -1;
  }

  PbPbUPCTrackTreeMessenger MTrackPbPbUPC(
      InputFile, myhelper::find_tree_from_list(InputFile, {"PbPbTracks/trackTree", "ppTracks/trackTree"}));
  PFTreeMessenger MPF(InputFile, myhelper::find_tree_from_list(InputFile, {"particleFlowAnalyser/pftree"}));
  TriggerTreeMessenger MTrigger(InputFile); // hltanalysis/HltTree

  TH2D H36Negative("h36Negative", "36-phi negative-eta occupancy;pf#eta;pf#phi",
                   HFPFBinEdges::EtaEdges36Negative().size() - 1, HFPFBinEdges::EtaEdges36Negative().data(),
                   HFPFBinEdges::PhiEdges36().size() - 1, HFPFBinEdges::PhiEdges36().data());
  TH2D H36Positive("h36Positive", "36-phi positive-eta occupancy;pf#eta;pf#phi",
                   HFPFBinEdges::EtaEdges36Positive().size() - 1, HFPFBinEdges::EtaEdges36Positive().data(),
                   HFPFBinEdges::PhiEdges36().size() - 1, HFPFBinEdges::PhiEdges36().data());
  TH2D H18NegativeMain("h18NegativeMain", "18-phi negative-eta occupancy main;pf#eta;pf#phi",
                       HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
                       HFPFBinEdges::PhiEdges18Main().size() - 1, HFPFBinEdges::PhiEdges18Main().data());
  TH2D H18NegativeWrapLow("h18NegativeWrapLow", "18-phi negative-eta occupancy wrap low;pf#eta;pf#phi",
                          HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
                          HFPFBinEdges::PhiEdges18WrapLow().size() - 1, HFPFBinEdges::PhiEdges18WrapLow().data());
  TH2D H18NegativeWrapHigh("h18NegativeWrapHigh", "18-phi negative-eta occupancy wrap high;pf#eta;pf#phi",
                           HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
                           HFPFBinEdges::PhiEdges18WrapHigh().size() - 1, HFPFBinEdges::PhiEdges18WrapHigh().data());
  TH2D H18PositiveMain("h18PositiveMain", "18-phi positive-eta occupancy main;pf#eta;pf#phi",
                       HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
                       HFPFBinEdges::PhiEdges18Main().size() - 1, HFPFBinEdges::PhiEdges18Main().data());
  TH2D H18PositiveWrapLow("h18PositiveWrapLow", "18-phi positive-eta occupancy wrap low;pf#eta;pf#phi",
                          HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
                          HFPFBinEdges::PhiEdges18WrapLow().size() - 1, HFPFBinEdges::PhiEdges18WrapLow().data());
  TH2D H18PositiveWrapHigh("h18PositiveWrapHigh", "18-phi positive-eta occupancy wrap high;pf#eta;pf#phi",
                           HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
                           HFPFBinEdges::PhiEdges18WrapHigh().size() - 1, HFPFBinEdges::PhiEdges18WrapHigh().data());
  std::vector<TH2D *> posHists = {&H36Positive, &H18PositiveMain, &H18PositiveWrapLow, &H18PositiveWrapHigh};
  std::vector<TH2D *> negHists = {&H36Negative, &H18NegativeMain, &H18NegativeWrapLow, &H18NegativeWrapHigh};
  std::vector<TH2D *> posEnergyHists = makeEnergyDistributionHistograms("pos", true);
  std::vector<TH2D *> negEnergyHists = makeEnergyDistributionHistograms("neg", false);

  // Make leading tower probability histograms
  TH2D H36Negative_leadingProb("h36Negative_leadingProb", "36-phi negative-eta leading tower prob. ;pf#eta;pf#phi",
                               HFPFBinEdges::EtaEdges36Negative().size() - 1, HFPFBinEdges::EtaEdges36Negative().data(),
                               HFPFBinEdges::PhiEdges36().size() - 1, HFPFBinEdges::PhiEdges36().data());
  TH2D H36Positive_leadingProb("h36Positive_leadingProb", "36-phi positive-eta leading tower prob. ;pf#eta;pf#phi",
                               HFPFBinEdges::EtaEdges36Positive().size() - 1, HFPFBinEdges::EtaEdges36Positive().data(),
                               HFPFBinEdges::PhiEdges36().size() - 1, HFPFBinEdges::PhiEdges36().data());
  TH2D H18NegativeMain_leadingProb(
      "h18NegativeMain_leadingProb", "18-phi negative-eta leading tower prob. main;pf#eta;pf#phi",
      HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
      HFPFBinEdges::PhiEdges18Main().size() - 1, HFPFBinEdges::PhiEdges18Main().data());
  TH2D H18NegativeWrapLow_leadingProb(
      "h18NegativeWrapLow_leadingProb", "18-phi negative-eta leading tower prob. wrap low;pf#eta;pf#phi",
      HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
      HFPFBinEdges::PhiEdges18WrapLow().size() - 1, HFPFBinEdges::PhiEdges18WrapLow().data());
  TH2D H18NegativeWrapHigh_leadingProb(
      "h18NegativeWrapHigh_leadingProb", "18-phi negative-eta leading tower prob. wrap high;pf#eta;pf#phi",
      HFPFBinEdges::EtaEdges18Negative().size() - 1, HFPFBinEdges::EtaEdges18Negative().data(),
      HFPFBinEdges::PhiEdges18WrapHigh().size() - 1, HFPFBinEdges::PhiEdges18WrapHigh().data());
  TH2D H18PositiveMain_leadingProb(
      "h18PositiveMain_leadingProb", "18-phi positive-eta leading tower prob. main;pf#eta;pf#phi",
      HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
      HFPFBinEdges::PhiEdges18Main().size() - 1, HFPFBinEdges::PhiEdges18Main().data());
  TH2D H18PositiveWrapLow_leadingProb(
      "h18PositiveWrapLow_leadingProb", "18-phi positive-eta leading tower prob. wrap low;pf#eta;pf#phi",
      HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
      HFPFBinEdges::PhiEdges18WrapLow().size() - 1, HFPFBinEdges::PhiEdges18WrapLow().data());
  TH2D H18PositiveWrapHigh_leadingProb(
      "h18PositiveWrapHigh_leadingProb", "18-phi positive-eta leading tower prob. wrap high;pf#eta;pf#phi",
      HFPFBinEdges::EtaEdges18Positive().size() - 1, HFPFBinEdges::EtaEdges18Positive().data(),
      HFPFBinEdges::PhiEdges18WrapHigh().size() - 1, HFPFBinEdges::PhiEdges18WrapHigh().data());
  std::vector<TH2D *> posHists_leadingProb = {&H36Positive_leadingProb, &H18PositiveMain_leadingProb,
                                              &H18PositiveWrapLow_leadingProb, &H18PositiveWrapHigh_leadingProb};
  std::vector<TH2D *> negHists_leadingProb = {&H36Negative_leadingProb, &H18NegativeMain_leadingProb,
                                              &H18NegativeWrapLow_leadingProb, &H18NegativeWrapHigh_leadingProb};

  int EntryCount = MTrigger.Tree->GetEntries();
  long long totalCandidates = 0;
  long long outOfRangeCandidates = 0;
  long long overflowCandidates = 0;
  long long selectedEvents = 0;

  for (long long iEvent = 0; iEvent < EntryCount; iEvent++) {
    MTrackPbPbUPC.GetEntry(iEvent);
    MPF.GetEntry(iEvent);
    MTrigger.GetEntry(iEvent);

    if (triggerChoice == 1 && !MTrigger.CheckTriggerStartWith("HLT_HIL1NotBptxOR_v")) {
      continue; // Trigger selection
    }
    selectedEvents++;

    double leading_Energy = -999.0;
    double leading_eta = -999.0;
    double leading_phi = -999.0;
    for (int iPf = 0; iPf < MPF.ID->size(); iPf++) {
      if (MPF.ID->at(iPf) != 6 && MPF.ID->at(iPf) != 7)
        continue;

      totalCandidates++;
      const double Eta = MPF.Eta->at(iPf);
      const double Phi = MPF.Phi->at(iPf);

      if (std::abs(Phi) > M_PI || std::abs(Eta) < 2.8 || std::abs(Eta) > 5.2)
        outOfRangeCandidates++;

      if (MPF.E->at(iPf) > leading_Energy) {
        leading_Energy = MPF.E->at(iPf);
        leading_eta = Eta;
        leading_phi = Phi;
      }

      int loc = determineHistLoc(Eta, Phi, do18BinMerging);
      int energyHistIndex = findEnergyHistogramIndex(Eta, Phi, do18BinMerging);
      if (energyHistIndex < 0) {
        std::cerr << "Error: Could not find energy histogram index for Eta: " << Eta << " Phi: " << Phi << std::endl;
        return -1;
      }
      if (Eta < 0) {
        negEnergyHists[energyHistIndex]->Fill(MPF.E->at(iPf), MPF.ID->at(iPf));
        negHists[loc]->Fill(Eta, Phi);
      } else {
        posEnergyHists[energyHistIndex]->Fill(MPF.E->at(iPf), MPF.ID->at(iPf));
        posHists[loc]->Fill(Eta, Phi);
      }
    }
    if (leading_Energy > 0) {
      int loc_leading = determineHistLoc(leading_eta, leading_phi, do18BinMerging);
      if (leading_phi >= HFPFBinEdges::PhiEdges18Main()[HFPFBinEdges::PhiEdges18Main().size() - 1] && do18BinMerging) {
        leading_phi *=
            -1; // Map wrap high region to wrap low region for leading tower probability histogram if merging is enabled
      }
      if (leading_eta < 0) {
        negHists_leadingProb[loc_leading]->Fill(leading_eta, leading_phi);
      } else {
        posHists_leadingProb[loc_leading]->Fill(leading_eta, leading_phi);
      }
    }
  }

  for (int iLoc = 0; iLoc < 4; iLoc++) {
    std::cout << "Overflow candidates in " << posHists[iLoc]->GetName() << ": "
              << posHists[iLoc]->GetEntries() - posHists[iLoc]->Integral() << std::endl;
    overflowCandidates += posHists[iLoc]->GetEntries() - posHists[iLoc]->Integral();
    std::cout << "Overflow candidates in " << negHists[iLoc]->GetName() << ": "
              << negHists[iLoc]->GetEntries() - negHists[iLoc]->Integral() << std::endl;
    overflowCandidates += negHists[iLoc]->GetEntries() - negHists[iLoc]->Integral();
  }

  std::cout << "Total candidates: " << totalCandidates << std::endl;
  std::cout << "Out of Range candidates" << outOfRangeCandidates << std::endl;
  std::cout << "Overflow candidates: " << overflowCandidates << std::endl;

  TH1D counts("counts", "Candidate counts;Category;Entries", 4, 0.0, 4.0);
  counts.SetBinContent(1, totalCandidates);
  counts.SetBinContent(2, outOfRangeCandidates);
  counts.SetBinContent(3, overflowCandidates);
  counts.SetBinContent(4, selectedEvents);
  counts.GetXaxis()->SetBinLabel(1, "Total");
  counts.GetXaxis()->SetBinLabel(2, "Out of Range");
  counts.GetXaxis()->SetBinLabel(3, "Overflow");
  counts.GetXaxis()->SetBinLabel(4, "Total Events");

  // Normalize by bin width to get density
  //  for (TH2D *hist : posHists) {
  //    for (int iBinX = 1; iBinX <= hist->GetNbinsX(); iBinX++) {
  //      for (int iBinY = 1; iBinY <= hist->GetNbinsY(); iBinY++) {
  //        const double BinContent = hist->GetBinContent(iBinX, iBinY);
  //        const double BinWidthX = hist->GetXaxis()->GetBinWidth(iBinX);
  //        const double BinWidthY = hist->GetYaxis()->GetBinWidth(iBinY);
  //        hist->SetBinContent(iBinX, iBinY, BinContent / (BinWidthX * BinWidthY));
  //      }
  //    }
  //  }
  //  for (TH2D *hist : negHists) {
  //    for (int iBinX = 1; iBinX <= hist->GetNbinsX(); iBinX++) {
  //      for (int iBinY = 1; iBinY <= hist->GetNbinsY(); iBinY++) {
  //        const double BinContent = hist->GetBinContent(iBinX, iBinY);
  //        const double BinWidthX = hist->GetXaxis()->GetBinWidth(iBinX);
  //        const double BinWidthY = hist->GetYaxis()->GetBinWidth(iBinY);
  //        hist->SetBinContent(iBinX, iBinY, BinContent / (BinWidthX * BinWidthY));
  //      }
  //    }
  //  }

  // Make average and standard deviation energy histograms
  std::vector<std::vector<TH2D *>> H36Negative_meta;
  std::vector<std::vector<TH2D *>> H36Positive_meta;
  std::vector<std::vector<TH2D *>> H18NegativeMain_meta;
  std::vector<std::vector<TH2D *>> H18NegativeWrapLow_meta;
  std::vector<std::vector<TH2D *>> H18NegativeWrapHigh_meta;
  std::vector<std::vector<TH2D *>> H18PositiveMain_meta;
  std::vector<std::vector<TH2D *>> H18PositiveWrapLow_meta;
  std::vector<std::vector<TH2D *>> H18PositiveWrapHigh_meta;
  std::vector<std::string> modeLabels = {"CountCrossCheck", "AvgEnergy", "StdDevEnergy"};
  std::vector<std::string> idLabels = {"All", "id6", "id7"};

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
                        H18PositiveWrapLow_meta[iMode][iID], H18PositiveWrapHigh_meta[iMode][iID], iMode, iID,
                        selectedEvents);
      fillAvgHistograms(negEnergyHists, H36Negative_meta[iMode][iID], H18NegativeMain_meta[iMode][iID],
                        H18NegativeWrapLow_meta[iMode][iID], H18NegativeWrapHigh_meta[iMode][iID], iMode, iID,
                        selectedEvents);
    }
  }

  TFile outputFile(OutputRootName.c_str(), "RECREATE");
  outputFile.cd();
  for (TH2D *hist : posHists)
    hist->Write();
  for (TH2D *hist : negHists)
    hist->Write();
  for (TH2D *hist : posHists_leadingProb) {
    hist->Write();
  }
  for (TH2D *hist : negHists_leadingProb) {
    hist->Write();
  }
  for (TH2D *hist : posEnergyHists)
    hist->Write();
  for (TH2D *hist : negEnergyHists)
    hist->Write();
  counts.Write();
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
