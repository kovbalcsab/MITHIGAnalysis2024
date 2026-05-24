#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"
#include "TTree.h"

#include "CommandLine.h"
#include "InfoManager.h"
#include "Messenger.h"
#include "ProgressBar.h"
#include "myhelper.h"

int main(int argc, char *argv[]);
double GetMaxEnergyHF(PFTreeMessenger *M, double etaMin, double etaMax, double phiMin, double phiMax,
                      double ptMin = 0.);

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  CL.SetVerbose();
  string InputFileName = CL.Get("Input");
  string OutputFileName = CL.Get("Output");

  double Fraction = CL.GetDouble("Fraction", 1.00);

  float ZDCMinus1nThreshold = CL.GetDouble("ZDCMinus1nThreshold", 1000.);
  float ZDCPlus1nThreshold = CL.GetDouble("ZDCPlus1nThreshold", 1100.);
  int TriggerChoice =
      CL.GetInt("TriggerChoice", 0); // 0 = isNotBptxOR, -1 = isUnpairedBunchBptxMinus, 1 = isUnpairedBunchBptxPlus
  int nTrkFilter =
      CL.GetInt("nTrkFilter", 0); // 0 = no filter, 1 = nTrkInAcceptanceHP == 0, -1 = nTrkInAcceptanceHP > 0
  if (TriggerChoice != 0 && TriggerChoice != -1 && TriggerChoice != 1 && TriggerChoice != -2) {
    std::cerr << "Error: Invalid TriggerChoice parameter value. Expected 0, -1, -2 or 1." << std::endl;
    return -1;
  }
  if (nTrkFilter != 0 && nTrkFilter != 1 && nTrkFilter != -1) {
    std::cerr << "Error: Invalid nTrkFilter parameter value. Expected 0, 1, or -1." << std::endl;
    return -1;
  }

  TFile *InputFile = TFile::Open(InputFileName.c_str());

  HiEventTreeMessenger MEvent(InputFile); // hiEvtAnalyzer/HiTree
  PbPbUPCTrackTreeMessenger MTrackPbPbUPC(
      InputFile, myhelper::find_tree_from_list(InputFile, {"PbPbTracks/trackTree", "ppTracks/trackTree"}));
  GenParticleTreeMessenger MGen(InputFile); // HiGenParticleAna/hi
  PFTreeMessenger MPF(InputFile, myhelper::find_tree_from_list(InputFile, {"particleFlowAnalyser/pftree"}));
  SkimTreeMessenger MSkim(InputFile);                        // skimanalysis/HltTree
  TriggerTreeMessenger MTrigger(InputFile);                  // hltanalysis/HltTree
  DzeroTreeMessenger MDzero(InputFile);                      // Dfinder/ntDkpi
  DfinderGenTreeMessenger MDzeroGen(InputFile);              // Dfinder/ntGen
  ZDCTreeMessenger MZDC(InputFile, "zdcanalyzer/zdcrechit"); // zdcanalyzer/zdcrechit
  METFilterTreeMessenger MMETFilter(InputFile);              // l1MetFilterRecoTree/MetFilterRecoTree

  TFile *OutputFile = new TFile(OutputFileName.c_str(), "RECREATE");

  std::vector<float> etaBorders = {3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2};
  std::vector<float> phiBorders = {-M_PI, -M_PI / 3 * 2, -M_PI / 3, 0, M_PI / 3, M_PI / 3 * 2, M_PI};

  GeneralInfoManager man(OutputFile, "InfoDir", false);
  TTimeStamp *currentTime = new TTimeStamp();

  man.AddSourceFile(InputFileName, currentTime);
  man.AddCutParameter("TriggerChoice", TriggerChoice, currentTime);
  man.AddCutParameter("nTrkFilter", nTrkFilter, currentTime);

  std::cout << "Parameters used for this analysis:" << std::endl;
  man.PrintInfo();

  std::vector<std::vector<TH1D *>> hHFEMaxPlusMaps, hHFEMaxMinusMaps;
  for (size_t iEta = 0; iEta < etaBorders.size() - 1; iEta++) {
    std::vector<TH1D *> hHFEMaxPlusPhiBins, hHFEMaxMinusPhiBins;
    for (size_t iPhi = 0; iPhi < phiBorders.size() - 1; iPhi++) {
      hHFEMaxPlusPhiBins.push_back(new TH1D(Form("hHFEMaxPlus_eta%.1f_%.1f_phi%.1f_%.1f", etaBorders[iEta],
                                                 etaBorders[iEta + 1], phiBorders[iPhi], phiBorders[iPhi + 1]),
                                            "", 300, 0, 100));
      hHFEMaxMinusPhiBins.push_back(new TH1D(Form("hHFEMaxMinus_eta%.1f_%.1f_phi%.1f_%.1f", -etaBorders[iEta + 1],
                                                  -etaBorders[iEta], phiBorders[iPhi], phiBorders[iPhi + 1]),
                                             "", 300, 0, 100));
    }
    hHFEMaxPlusMaps.push_back(hHFEMaxPlusPhiBins);
    hHFEMaxMinusMaps.push_back(hHFEMaxMinusPhiBins);
  }
  TH1D *hHFEMaxPlusLeading = new TH1D("hHFEMaxPlusLeading", "", 300, 0, 100);
  TH1D *hHFEMaxMinusLeading = new TH1D("hHFEMaxMinusLeading", "", 300, 0, 100);

  TH1D *hNumberOfEventsAfterCuts = new TH1D("hNumberOfEventsAfterCuts", "", 4, -0.5, 3.5);
  hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(1, "NoCuts");
  hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(2, "Trigger");
  hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(3, "TrackFilter");
  hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(4, "ZDC");
  hNumberOfEventsAfterCuts->Sumw2();

  int EntryCount = MEvent.GetEntries() * Fraction;
  ProgressBar Bar(cout, EntryCount);

  // Event loop with event selection
  for (int iEntry = 0; iEntry < EntryCount; iEntry++) {
    if (iEntry % (EntryCount / 100) == 0) {
      Bar.Update(iEntry);
      Bar.Print();
    }

    MEvent.GetEntry(iEntry);
    MGen.GetEntry(iEntry);
    MTrackPbPbUPC.GetEntry(iEntry);
    MPF.GetEntry(iEntry);
    MSkim.GetEntry(iEntry);
    MTrigger.GetEntry(iEntry);
    MDzero.GetEntry(iEntry);
    MZDC.GetEntry(iEntry);
    MMETFilter.GetEntry(iEntry);
    hNumberOfEventsAfterCuts->Fill(0); // NoCuts

    // Apply event selection criteria here
    if ((TriggerChoice == 0 && !MTrigger.CheckTriggerStartWith("HLT_HIL1NotBptxOR_v")) ||
        (TriggerChoice == -1 && !MTrigger.CheckTriggerStartWith("HLT_HIL1UnpairedBunchBptxMinus_v")) ||
        (TriggerChoice == 1 && !MTrigger.CheckTriggerStartWith("HLT_HIL1UnpairedBunchBptxPlus_v"))) {
      continue; // Trigger selection
    }
    hNumberOfEventsAfterCuts->Fill(1); // Trigger passed

    // Count tracks in acceptance with high purity
    int nTrackInAcceptanceHP = 0;
    for (int iTrack = 0; iTrack < MTrackPbPbUPC.nTrk; iTrack++) {
      if (MTrackPbPbUPC.trkPt->at(iTrack) <= 0.3)
        continue;
      if (fabs(MTrackPbPbUPC.trkEta->at(iTrack)) >= 2.4)
        continue;
      if (MTrackPbPbUPC.highPurity->at(iTrack) == false)
        continue;
      nTrackInAcceptanceHP++;
    }
    if ((nTrkFilter == 1 && nTrackInAcceptanceHP != 0) || (nTrkFilter == -1 && nTrackInAcceptanceHP == 0)) {
      continue; // Track filter selection
    }
    hNumberOfEventsAfterCuts->Fill(2); // Track filter passed

    if (MZDC.sumPlus > ZDCPlus1nThreshold || MZDC.sumMinus > ZDCMinus1nThreshold)
      continue;                        // ZDC energy cut
    hNumberOfEventsAfterCuts->Fill(3); // ZDC cut passed

    // Fill HF E_max maps for empty BX events
    for (size_t iEta = 0; iEta < etaBorders.size() - 1; iEta++) {
      for (size_t iPhi = 0; iPhi < phiBorders.size() - 1; iPhi++) {
        double HFEMaxPlus =
            GetMaxEnergyHF(&MPF, etaBorders[iEta], etaBorders[iEta + 1], phiBorders[iPhi], phiBorders[iPhi + 1], 0.0);
        double HFEMaxMinus =
            GetMaxEnergyHF(&MPF, -etaBorders[iEta + 1], -etaBorders[iEta], phiBorders[iPhi], phiBorders[iPhi + 1], 0.0);
        hHFEMaxPlusMaps[iEta][iPhi]->Fill(HFEMaxPlus);
        hHFEMaxMinusMaps[iEta][iPhi]->Fill(HFEMaxMinus);
      }
    }
    hHFEMaxPlusLeading->Fill(GetMaxEnergyHF(&MPF, 3.0, 5.2, -M_PI, M_PI, 0.0));
    hHFEMaxMinusLeading->Fill(GetMaxEnergyHF(&MPF, -5.2, -3.0, -M_PI, M_PI, 0.0));
  }

  InputFile->Close();

  // Write histograms to output file
  OutputFile->cd();
  for (size_t iEta = 0; iEta < etaBorders.size() - 1; iEta++) {
    for (size_t iPhi = 0; iPhi < phiBorders.size() - 1; iPhi++) {
      hHFEMaxPlusMaps[iEta][iPhi]->Write();
      hHFEMaxMinusMaps[iEta][iPhi]->Write();
    }
  }
  hHFEMaxPlusLeading->Write();
  hHFEMaxMinusLeading->Write();
  hNumberOfEventsAfterCuts->Write();
  man.SaveToFile();
  OutputFile->Close();
  std::cout << "Analysis complete. Output written to " << OutputFileName << std::endl;
  return 0;
}

// ============================================================================ //
// Function to Retrieve Maximum Energy in HF Region within Specified Eta Range
// ============================================================================ //
double GetMaxEnergyHF(PFTreeMessenger *M, double etaMin, double etaMax, double phiMin, double phiMax, double ptMin) {
  if (M == nullptr)
    return -1;
  if (M->Tree == nullptr)
    return -1;

  double EMax = 0;
  for (int iPF = 0; iPF < M->ID->size(); iPF++) {
    if (M->PT->at(iPF) <= ptMin)
      continue;
    if ((M->ID->at(iPF) == 6 || M->ID->at(iPF) == 7) && M->Eta->at(iPF) > etaMin && M->Eta->at(iPF) < etaMax &&
        M->Phi->at(iPF) > phiMin && M->Phi->at(iPF) < phiMax) {
      if (M->E->at(iPF) > EMax)
        EMax = M->E->at(iPF);
    }
  }
  return EMax;
}
