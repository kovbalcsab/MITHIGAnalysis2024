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
  const std::string InputSpec = CL.Get("Input", "staging_first100_files.txt");
  const std::string TreeName = CL.Get("Tree", "particleFlowAnalyser/pftree");
  const std::string OutputRootName = CL.Get("OutputRoot", "hf_pf_binning_check.root");
  const bool do18BinMerging = CL.GetBool("do18BinMerging", true);

  const std::vector<std::string> InputFiles = ExpandInputFiles(InputSpec);
  TChain Chain(TreeName.c_str());
  for (const std::string &FileName : InputFiles)
    Chain.Add(FileName.c_str());

  std::vector<int> *PFID = nullptr;
  std::vector<float> *PFEta = nullptr;
  std::vector<float> *PFPhi = nullptr;
  std::vector<float> *PFEnergy = nullptr;
  Chain.SetBranchAddress("pfId", &PFID);
  Chain.SetBranchAddress("pfEta", &PFEta);
  Chain.SetBranchAddress("pfPhi", &PFPhi);
  Chain.SetBranchAddress("pfE", &PFEnergy);

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

  long long totalCandidates = 0;
  long long outOfRangeCandidates = 0;
  long long overflowCandidates = 0;
  for (long long iEvent = 0; iEvent < Chain.GetEntries(); iEvent++) {
    Chain.GetEntry(iEvent);
    for (int iPf = 0; iPf < PFID->size(); iPf++) {
      if (PFID->at(iPf) != 6 && PFID->at(iPf) != 7)
        continue;

      totalCandidates++;
      const double Eta = PFEta->at(iPf);
      const double Phi = PFPhi->at(iPf);

      if (std::abs(Phi) > M_PI || std::abs(Eta) < 2.8 || std::abs(Eta) > 5.2)
        outOfRangeCandidates++;

      int loc = determineHistLoc(Eta, Phi, do18BinMerging);
      int energyHistIndex = findEnergyHistogramIndex(Eta, Phi, do18BinMerging);
      if (energyHistIndex < 0) {
        std::cerr << "Error: Could not find energy histogram index for Eta: " << Eta << " Phi: " << Phi << std::endl;
        return -1;
      }
      if (Eta < 0) {
        negEnergyHists[energyHistIndex]->Fill(PFEnergy->at(iPf), PFID->at(iPf));
        negHists[loc]->Fill(Eta, Phi);
      } else {
        posEnergyHists[energyHistIndex]->Fill(PFEnergy->at(iPf), PFID->at(iPf));
        posHists[loc]->Fill(Eta, Phi);
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

  TH1D counts("counts", "Candidate counts;Category;Entries", 3, 0.0, 3.0);
  counts.SetBinContent(1, totalCandidates);
  counts.SetBinContent(2, outOfRangeCandidates);
  counts.SetBinContent(3, overflowCandidates);
  counts.GetXaxis()->SetBinLabel(1, "Total");
  counts.GetXaxis()->SetBinLabel(2, "Out of Range");
  counts.GetXaxis()->SetBinLabel(3, "Overflow");

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
                        H18PositiveWrapLow_meta[iMode][iID], H18PositiveWrapHigh_meta[iMode][iID], iMode, iID);
      fillAvgHistograms(negEnergyHists, H36Negative_meta[iMode][iID], H18NegativeMain_meta[iMode][iID],
                        H18NegativeWrapLow_meta[iMode][iID], H18NegativeWrapHigh_meta[iMode][iID], iMode, iID);
    }
  }

  TFile outputFile(OutputRootName.c_str(), "RECREATE");
  outputFile.cd();
  for (TH2D *hist : posHists)
    hist->Write();
  for (TH2D *hist : negHists)
    hist->Write();
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
