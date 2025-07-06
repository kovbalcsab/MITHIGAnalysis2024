#include <algorithm>
#include <iostream>
#include <set>
#include <ctime>
#include <cstdlib>
using namespace std;

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "CommandLine.h"
#include "CommonFunctions.h"
#include "Messenger.h"
#include "ProgressBar.h"

#include "TrackResidualCorrector.h"
#include "tnp_weight.h"
#include "trackingEfficiency2017pp.h"
#include "trackingEfficiency2018PbPb.h"
#include "trackingEfficiency2023PbPb.h"
#include "trackingEfficiency2024ppref.h"
#include "trackingEfficiency2025OO.h"

#include "include/cent_OO_hijing_PF.h"
#include "include/skimSelectionBits_OO_PP.h"
#include "include/parseFSCandPPSInfo.h"

int main(int argc, char *argv[]) {
  string VersionString = "V8";
  srand(time(0));

  CommandLine CL(argc, argv);
  vector<string> InputFileNames = CL.GetStringVector("Input");
  string OutputFileName = CL.Get("Output");

  bool DoGenLevel = CL.GetBool("DoGenLevel", false);
  bool IsData = CL.GetBool("IsData", false);
  bool IsPP = CL.GetBool("IsPP", false);
  int Year = CL.GetInt("Year", 2024);

  vector<double> Fractions = CL.GetDoubleVector("Fractions");
  // bool ApplyZDCGapRejection = CL.GetBool("ApplyZDCGapRejection", false);
  // 0 for HIJING 00, 1 for Starlight SD, 2 for Starlight DD, 4 for HIJING alpha-O, -1 for data
  bool DebugMode = CL.GetBool("DebugMode", false);
  bool includeFSCandPPSMode = CL.GetBool("includeFSCandPPSMode", false);
  int saveTriggerBitsMode = CL.GetInt("saveTriggerBitsMode", 0);
  int nEvents = CL.GetInt("nEvents", 0);

  if (InputFileNames.size() > (Fractions.size()+1)) {
    cerr << "Error: Number of input files (" << InputFileNames.size() << ") exceeds number of fractions+1 (" << (Fractions.size()+1) << ")." << endl;
    return -1;
  }

  TFile OutputFile(OutputFileName.c_str(), "RECREATE");
  OutputFile.cd();

  TTree Tree("Tree", Form("Tree for UPC Dzero analysis (%s)", VersionString.c_str()));
  TTree InfoTree("InfoTree", "Information");
  vector<TFile*> InputFiles;
  ChargedHadronRAATreeMessenger MChargedHadronRAAOutput;
  vector<ChargedHadronRAATreeMessenger*> MChargedHadronRAAInputs;
  MChargedHadronRAAOutput.SetBranch(&Tree, saveTriggerBitsMode, DebugMode, includeFSCandPPSMode);

  float sumFraction = 0;
  for (size_t i = 0; i < Fractions.size(); i++) {
    if (Fractions[i] <= 0 || Fractions[i] > 1) {
      cerr << "Error: Fraction " << Fractions[i] << " is not in the range (0, 1]." << endl;
      return -1;
    }
    sumFraction += Fractions[i];
  }
  if (sumFraction > 1) {
    cerr << "Error: Sum of fractions exceeds 1.0." << endl;
    return -1;
  }
  if (Fractions.size() < InputFileNames.size()) {
    Fractions.push_back(1.0 - sumFraction); // Add the last fraction to make the total sum equal to 1
  }

  int sumOfEntries=0;
  bool abortShuffling = false;
  for (int iFile=0; iFile < InputFileNames.size(); iFile++) {
    InputFiles.push_back(new TFile(InputFileNames[iFile].c_str()));
    MChargedHadronRAAInputs.push_back(new ChargedHadronRAATreeMessenger(*InputFiles[iFile], string("Tree")));
    sumOfEntries += MChargedHadronRAAInputs[iFile]->GetEntries();
    if (MChargedHadronRAAInputs[iFile]->GetEntries() < nEvents*Fractions[iFile]) { abortShuffling = true; }
  }
  
  int runningSum = 0;
  if (!abortShuffling) {
    cout << "Total number of entries in all input files: " << sumOfEntries << endl;
    for (int iFile=0; iFile < InputFileNames.size(); iFile++) {
      int EntryCount = MChargedHadronRAAInputs[iFile]->GetEntries(); 
      int desiredEntries = nEvents * Fractions[iFile];
      
      if (runningSum >= nEvents) {
        break; // Stop processing if we have reached the desired number of events
      }

      /////////////////////////////////
      //////// Main Event Loop ////////
      /////////////////////////////////
      for (int iE = 0; iE < EntryCount; iE++) {
        if (runningSum >= nEvents) {
          break; // Stop processing if we have reached the desired number of events
        }
        if (rand() % EntryCount < desiredEntries) {
          runningSum++;
        } else {
          continue; // Skip this event if it is not selected
        }
        MChargedHadronRAAInputs[iFile]->GetEntry(iE);        

        ////////////////////////////////////////
        ////////// Global event stuff //////////
        ////////////////////////////////////////
        MChargedHadronRAAOutput.sampleType = MChargedHadronRAAInputs[iFile]->sampleType;
        MChargedHadronRAAOutput.Run = MChargedHadronRAAInputs[iFile]->Run;
        MChargedHadronRAAOutput.Lumi = MChargedHadronRAAInputs[iFile]->Lumi;
        MChargedHadronRAAOutput.Event = MChargedHadronRAAInputs[iFile]->Event;
        MChargedHadronRAAOutput.hiHF_pf = MChargedHadronRAAInputs[iFile]->hiHF_pf;
        MChargedHadronRAAOutput.hiHFPlus_pf = MChargedHadronRAAInputs[iFile]->hiHFPlus_pf;
        MChargedHadronRAAOutput.hiHFMinus_pf = MChargedHadronRAAInputs[iFile]->hiHFMinus_pf;
        MChargedHadronRAAOutput.hiBin = MChargedHadronRAAInputs[iFile]->hiBin;

        MChargedHadronRAAOutput.Ncoll = MChargedHadronRAAInputs[iFile]->Ncoll;
        MChargedHadronRAAOutput.Npart = MChargedHadronRAAInputs[iFile]->Npart;
        ////////////////////////////
        ////////// Vertex //////////
        ////////////////////////////
        MChargedHadronRAAOutput.VX = MChargedHadronRAAInputs[iFile]->VX;
        MChargedHadronRAAOutput.VY = MChargedHadronRAAInputs[iFile]->VY;
        MChargedHadronRAAOutput.VZ = MChargedHadronRAAInputs[iFile]->VZ;
        MChargedHadronRAAOutput.VXError = MChargedHadronRAAInputs[iFile]->VXError;
        MChargedHadronRAAOutput.VYError = MChargedHadronRAAInputs[iFile]->VYError;
        MChargedHadronRAAOutput.VZError = MChargedHadronRAAInputs[iFile]->VZError;
        MChargedHadronRAAOutput.isFakeVtx = MChargedHadronRAAInputs[iFile]->isFakeVtx;
        MChargedHadronRAAOutput.ptSumVtx = MChargedHadronRAAInputs[iFile]->ptSumVtx;
        MChargedHadronRAAOutput.nTracksVtx = MChargedHadronRAAInputs[iFile]->nTracksVtx;
        MChargedHadronRAAOutput.chi2Vtx = MChargedHadronRAAInputs[iFile]->chi2Vtx;
        MChargedHadronRAAOutput.ndofVtx = MChargedHadronRAAInputs[iFile]->ndofVtx;
        MChargedHadronRAAOutput.bestVtxIndx = MChargedHadronRAAInputs[iFile]->bestVtxIndx;
        
        MChargedHadronRAAOutput.nVtx = MChargedHadronRAAInputs[iFile]->nVtx;
        /////////////////////////////////////
        ////////// Event selection //////////
        /////////////////////////////////////
        MChargedHadronRAAOutput.ZDCsumPlus = MChargedHadronRAAInputs[iFile]->ZDCsumPlus;
        MChargedHadronRAAOutput.ZDCsumMinus = MChargedHadronRAAInputs[iFile]->ZDCsumMinus;
        MChargedHadronRAAOutput.ClusterCompatibilityFilter = MChargedHadronRAAInputs[iFile]->ClusterCompatibilityFilter;
        MChargedHadronRAAOutput.PVFilter = MChargedHadronRAAInputs[iFile]->PVFilter;
        MChargedHadronRAAOutput.mMaxL1HFAdcPlus = MChargedHadronRAAInputs[iFile]->mMaxL1HFAdcPlus;
        MChargedHadronRAAOutput.mMaxL1HFAdcMinus = MChargedHadronRAAInputs[iFile]->mMaxL1HFAdcMinus;
        MChargedHadronRAAOutput.VZ_pf = MChargedHadronRAAInputs[iFile]->VZ_pf;

        MChargedHadronRAAOutput.eventWeight = MChargedHadronRAAInputs[iFile]->eventWeight;
        if (IsData == true) {
          if (saveTriggerBitsMode == 1) { // OO triggers
              MChargedHadronRAAOutput.HLT_OxySingleJet16_ZDC1nAsymXOR_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxySingleJet16_ZDC1nAsymXOR_v1;
              MChargedHadronRAAOutput.HLT_OxySingleJet16_ZDC1nXOR_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxySingleJet16_ZDC1nXOR_v1;
              MChargedHadronRAAOutput.HLT_OxySingleJet24_ZDC1nAsymXOR_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxySingleJet24_ZDC1nAsymXOR_v1;
              MChargedHadronRAAOutput.HLT_OxySingleJet24_ZDC1nXOR_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxySingleJet24_ZDC1nXOR_v1;
              MChargedHadronRAAOutput.HLT_OxyZDC1nOR_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxyZDC1nOR_v1;
              MChargedHadronRAAOutput.HLT_OxyZeroBias_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxyZeroBias_v1;
              MChargedHadronRAAOutput.HLT_MinimumBiasHF_OR_BptxAND_v1 = MChargedHadronRAAInputs[iFile]->HLT_MinimumBiasHF_OR_BptxAND_v1;
              MChargedHadronRAAOutput.HLT_OxyL1SingleJet20_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxyL1SingleJet20_v1;
          } else if (saveTriggerBitsMode == 2) { // pO triggers
              MChargedHadronRAAOutput.HLT_OxyZeroBias_v1 = MChargedHadronRAAInputs[iFile]->HLT_OxyZeroBias_v1;
              MChargedHadronRAAOutput.HLT_OxyZDC1nOR_v1= MChargedHadronRAAInputs[iFile]->HLT_OxyZDC1nOR_v1;
              MChargedHadronRAAOutput.HLT_OxySingleMuOpen_NotMBHF2OR_v1= MChargedHadronRAAInputs[iFile]->HLT_OxySingleMuOpen_NotMBHF2OR_v1;
              MChargedHadronRAAOutput.HLT_OxySingleJet8_ZDC1nAsymXOR_v1= MChargedHadronRAAInputs[iFile]->HLT_OxySingleJet8_ZDC1nAsymXOR_v1;
              MChargedHadronRAAOutput.HLT_OxyNotMBHF2_v1= MChargedHadronRAAInputs[iFile]->HLT_OxyNotMBHF2_v1;
              MChargedHadronRAAOutput.HLT_OxyZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v1= MChargedHadronRAAInputs[iFile]->HLT_OxyZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v1;
              MChargedHadronRAAOutput.HLT_OxyZeroBias_MinPixelCluster400_v1= MChargedHadronRAAInputs[iFile]->HLT_OxyZeroBias_MinPixelCluster400_v1;
              MChargedHadronRAAOutput.HLT_MinimumBiasHF_OR_BptxAND_v1= MChargedHadronRAAInputs[iFile]->HLT_MinimumBiasHF_OR_BptxAND_v1;
              MChargedHadronRAAOutput.HLT_MinimumBiasHF_AND_BptxAND_v1= MChargedHadronRAAInputs[iFile]->HLT_MinimumBiasHF_AND_BptxAND_v1;
          }
        }
           
        
        MChargedHadronRAAOutput.HFEMaxPlus = MChargedHadronRAAInputs[iFile]->HFEMaxPlus;
        MChargedHadronRAAOutput.HFEMaxPlus2 = MChargedHadronRAAInputs[iFile]->HFEMaxPlus2;
        MChargedHadronRAAOutput.HFEMaxPlus3 = MChargedHadronRAAInputs[iFile]->HFEMaxPlus3;
        MChargedHadronRAAOutput.HFEMaxMinus = MChargedHadronRAAInputs[iFile]->HFEMaxMinus;
        MChargedHadronRAAOutput.HFEMaxMinus2 = MChargedHadronRAAInputs[iFile]->HFEMaxMinus2;
        MChargedHadronRAAOutput.HFEMaxMinus3 = MChargedHadronRAAInputs[iFile]->HFEMaxMinus3;

        int NTrack = MChargedHadronRAAInputs[iFile]->trkEta->size();
        MChargedHadronRAAOutput.nTrk = NTrack;
        for (int iTrack = 0; iTrack < NTrack; iTrack++) {
          MChargedHadronRAAOutput.trkEta->push_back(MChargedHadronRAAInputs[iFile]->trkEta->at(iTrack));
          MChargedHadronRAAOutput.trkPt->push_back(MChargedHadronRAAInputs[iFile]->trkPt->at(iTrack));
          MChargedHadronRAAOutput.trkPhi->push_back(MChargedHadronRAAInputs[iFile]->trkPhi->at(iTrack));
          MChargedHadronRAAOutput.trkPtError->push_back(MChargedHadronRAAInputs[iFile]->trkPtError->at(iTrack));
          MChargedHadronRAAOutput.highPurity->push_back(MChargedHadronRAAInputs[iFile]->highPurity->at(iTrack));
          MChargedHadronRAAOutput.trkDxyAssociatedVtx->push_back(MChargedHadronRAAInputs[iFile]->trkDxyAssociatedVtx->at(iTrack));
          MChargedHadronRAAOutput.trkDzAssociatedVtx->push_back(MChargedHadronRAAInputs[iFile]->trkDzAssociatedVtx->at(iTrack));
          MChargedHadronRAAOutput.trkDxyErrAssociatedVtx->push_back(MChargedHadronRAAInputs[iFile]->trkDxyErrAssociatedVtx->at(iTrack));
          MChargedHadronRAAOutput.trkDzErrAssociatedVtx->push_back(MChargedHadronRAAInputs[iFile]->trkDzErrAssociatedVtx->at(iTrack));
          MChargedHadronRAAOutput.trkAssociatedVtxIndx->push_back(MChargedHadronRAAInputs[iFile]->trkAssociatedVtxIndx->at(iTrack));
          MChargedHadronRAAOutput.trkCharge->push_back(MChargedHadronRAAInputs[iFile]->trkCharge->at(iTrack));
          MChargedHadronRAAOutput.trkNHits->push_back(MChargedHadronRAAInputs[iFile]->trkNHits->at(iTrack));
          MChargedHadronRAAOutput.trkNPixHits->push_back(MChargedHadronRAAInputs[iFile]->trkNPixHits->at(iTrack));
          MChargedHadronRAAOutput.trkNLayers->push_back(MChargedHadronRAAInputs[iFile]->trkNLayers->at(iTrack));
          MChargedHadronRAAOutput.trkNormChi2->push_back(MChargedHadronRAAInputs[iFile]->trkNormChi2->at(iTrack));
          MChargedHadronRAAOutput.pfEnergy->push_back(MChargedHadronRAAInputs[iFile]->pfEnergy->at(iTrack));

          MChargedHadronRAAOutput.trackWeight->push_back(MChargedHadronRAAInputs[iFile]->trackWeight->at(iTrack));
        } // end of loop over tracks (gen or reco)
        MChargedHadronRAAOutput.leadingPtEta1p0_sel = MChargedHadronRAAInputs[iFile]->leadingPtEta1p0_sel;
        MChargedHadronRAAOutput.multipicityEta1p0 = MChargedHadronRAAInputs[iFile]->multipicityEta1p0;
        MChargedHadronRAAOutput.multipicityEta2p4 = MChargedHadronRAAInputs[iFile]->multipicityEta2p4;

        ////////////////////////////
        ///// Debug variables //////
        ////////////////////////////

        if (DebugMode) {
          for (int iDebVtx = 0; iDebVtx < MChargedHadronRAAInputs[iFile]->nVtx; iDebVtx++) {
            MChargedHadronRAAOutput.AllxVtx->push_back(MChargedHadronRAAInputs[iFile]->AllxVtx->at(iDebVtx));
            MChargedHadronRAAOutput.AllyVtx->push_back(MChargedHadronRAAInputs[iFile]->AllyVtx->at(iDebVtx));
            MChargedHadronRAAOutput.AllzVtx->push_back(MChargedHadronRAAInputs[iFile]->AllzVtx->at(iDebVtx));
            MChargedHadronRAAOutput.AllxVtxError->push_back(MChargedHadronRAAInputs[iFile]->AllxVtxError->at(iDebVtx));
            MChargedHadronRAAOutput.AllyVtxError->push_back(MChargedHadronRAAInputs[iFile]->AllyVtxError->at(iDebVtx));
            MChargedHadronRAAOutput.AllzVtxError->push_back(MChargedHadronRAAInputs[iFile]->AllzVtxError->at(iDebVtx));
            MChargedHadronRAAOutput.AllisFakeVtx->push_back(MChargedHadronRAAInputs[iFile]->AllisFakeVtx->at(iDebVtx));
            MChargedHadronRAAOutput.AllnTracksVtx->push_back(MChargedHadronRAAInputs[iFile]->AllnTracksVtx->at(iDebVtx));
            MChargedHadronRAAOutput.Allchi2Vtx->push_back(MChargedHadronRAAInputs[iFile]->Allchi2Vtx->at(iDebVtx));
            MChargedHadronRAAOutput.AllndofVtx->push_back(MChargedHadronRAAInputs[iFile]->AllndofVtx->at(iDebVtx));
            MChargedHadronRAAOutput.AllptSumVtx->push_back(MChargedHadronRAAInputs[iFile]->AllptSumVtx->at(iDebVtx));
          }
        }

        ////////////////////////////
        /// PPS & FSC variables ////
        ////////////////////////////

        //currently not implemented

        ////////////////////////////////////////
        ///// Fill default selection bits //////
        ////////////////////////////////////////

        MChargedHadronRAAOutput.passBaselineEventSelection = MChargedHadronRAAInputs[iFile]->passBaselineEventSelection; 
        // Fill HF selection bits
        MChargedHadronRAAOutput.passL1HFAND_16_Offline = MChargedHadronRAAInputs[iFile]->passL1HFAND_16_Offline;  
        MChargedHadronRAAOutput.passL1HFOR_16_Offline = MChargedHadronRAAInputs[iFile]->passL1HFOR_16_Offline;  
        MChargedHadronRAAOutput.passL1HFAND_14_Offline = MChargedHadronRAAInputs[iFile]->passL1HFAND_14_Offline;  
        MChargedHadronRAAOutput.passL1HFOR_14_Offline = MChargedHadronRAAInputs[iFile]->passL1HFOR_14_Offline;

        MChargedHadronRAAOutput.passL1HFAND_16_Online = MChargedHadronRAAInputs[iFile]->passL1HFAND_16_Online;  
        MChargedHadronRAAOutput.passL1HFOR_16_Online = MChargedHadronRAAInputs[iFile]->passL1HFOR_16_Online; 
        MChargedHadronRAAOutput.passL1HFAND_14_Online = MChargedHadronRAAInputs[iFile]->passL1HFAND_14_Online;  
        MChargedHadronRAAOutput.passL1HFOR_14_Online = MChargedHadronRAAInputs[iFile]->passL1HFOR_14_Online; 
            
        MChargedHadronRAAOutput.FillEntry();
      }
        
    }
  }

  OutputFile.cd();
  Tree.Write();
  InfoTree.Write();

  OutputFile.Close();

  for (int iFile = 0; iFile < InputFiles.size(); iFile++) {
    InputFiles[iFile]->Close();
    delete MChargedHadronRAAInputs[iFile];
  }

  return 0;
}