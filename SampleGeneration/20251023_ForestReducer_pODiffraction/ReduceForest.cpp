#include <algorithm>
#include <iostream>
#include <set>
using namespace std;

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "CommandLine.h"
#include "CommonFunctions.h"
#include "Messenger.h"
#include "ProgressBar.h"

#include "include/cent_OO_hijing_PF.h"
#include "include/parseFSCandPPSInfo.h"
#include "include/skimSelectionBits_pO.h"

bool logical_or_vectBool(std::vector<bool> *vec) {
  return std::any_of(vec->begin(), vec->end(), [](bool b) { return b; });
}

// Helper function to convert a string to lowercase
std::string toLower(const std::string &str) {
  std::string lowerStr = str;
  std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
  return lowerStr;
}

int main(int argc, char *argv[]) {
  string VersionString = "V8";

  CommandLine CL(argc, argv);
  vector<string> InputFileNames = CL.GetStringVector("Input");
  string OutputFileName = CL.Get("Output");
  bool IsData = CL.GetBool("IsData", false);
  string CollisionSystem = CL.Get("CollisionSystem", "pO");
  double Fraction = CL.GetDouble("Fraction", 1.00);
  int ApplyTriggerRejection =
      CL.GetInteger("ApplyTriggerRejection", 0); // trigger = 0 for no rejection, 1 for ZeroBias, 2 for MinBias
  bool ApplyEventRejection = CL.GetBool("ApplyEventRejection", false);

  int sampleType = CL.GetInteger(
      "sampleType", -1); // -1 for data
  string PFTreeName = CL.Get("PFTree", "particleFlowAnalyser/pftree");
  string ZDCTreeName = CL.Get("ZDCTree", "zdcanalyzer/zdcrechit");
  string PPSTreeName = CL.Get("PPSTree", "ppsanalyzer/ppstracks");
  string FSCTreeName = CL.Get("FSCTree", "fscanalyzer/fscdigi");
  bool HideProgressBar = CL.GetBool("HideProgressBar", false);
  bool DebugMode = CL.GetBool("DebugMode", false);
  bool includeFSCandPPSMode = CL.GetBool("includeFSCandPPSMode", true);

  int saveTriggerBitsMode = 0;
  if (CollisionSystem != "OO" && CollisionSystem != "pO" && CollisionSystem != "NeNe") {
    std::cout << "ERROR: Collision system must be OO or pO. Exiting." << std::endl;
    return 1;
  }
  if (CollisionSystem == "OO" || CollisionSystem == "NeNe") {
    saveTriggerBitsMode = 1;
  } // save trigger bits for OO and NeNe
  else if (CollisionSystem == "pO") {
    saveTriggerBitsMode = 2;
  } // save trigger bits for pO

  TFile OutputFile(OutputFileName.c_str(), "RECREATE");
  TTree Tree("Tree", Form("Tree for pO diffractive analysis :) (%s)", VersionString.c_str()));
  TTree InfoTree("InfoTree", "Information");
  pODiffractiveTreeMessenger MpODiffractive;
  MpODiffractive.SetBranch(&Tree, saveTriggerBitsMode, DebugMode, includeFSCandPPSMode);

  for (string InputFileName : InputFileNames) {
    TFile InputFile(InputFileName.c_str());

    HiEventTreeMessenger MEvent(InputFile); // hiEvtAnalyzer/HiTree
    PPTrackTreeMessenger MTrack(InputFile, "ppTracks/trackTree");
    GenParticleTreeMessenger MGen(InputFile);      // HiGenParticleAna/hi
    SkimTreeMessenger MSkim(InputFile);            // skimanalysis/HltTree
    HFAdcMessenger MHFAdc(InputFile);              // HFAdcana/adc
    PFTreeMessenger MPF(InputFile, PFTreeName);    // particleFlowAnalyser/pftree
    ZDCTreeMessenger MZDC(InputFile, ZDCTreeName); // zdcanalyzer/zdcrechit
    TriggerTreeMessenger MTrigger(InputFile);      // hltanalysis/HltTree
    PPSTreeMessenger MPPS(InputFile, PPSTreeName); // ppsanalyzer/ppstracks
    FSCTreeMessenger MFSC(InputFile, FSCTreeName); // fscanalyzer/fscdigi
    // METFilterTreeMessenger MMETFilter(InputFile); // l1MetFilterRecoTree/MetFilterRecoTree

    int EntryCount = MEvent.GetEntries() * Fraction;
    ProgressBar Bar(cout, EntryCount);
    if (!HideProgressBar) {
      Bar.SetStyle(-1);
    }

    /////////////////////////////////
    //////// Main Event Loop ////////
    /////////////////////////////////
    for (int iE = 0; iE < EntryCount; iE++) {
      if (!HideProgressBar && (EntryCount < 300 || (iE % (EntryCount / 250)) == 0)) {
        Bar.Update(iE);
        Bar.Print();
      }
      MEvent.GetEntry(iE);
      MGen.GetEntry(iE);
      MTrack.GetEntry(iE);
      MSkim.GetEntry(iE);

      MHFAdc.GetEntry(iE);
      MPF.GetEntry(iE);
      MZDC.GetEntry(iE);
      MPPS.GetEntry(iE);
      MFSC.GetEntry(iE);
      MTrigger.GetEntry(iE);
      // MMETFilter.GetEntry(iE);

      ////////////////////////////////////////
      ////////// Global event stuff //////////
      ////////////////////////////////////////
      MpODiffractive.sampleType = sampleType;
      MpODiffractive.Run = MEvent.Run;
      MpODiffractive.Lumi = MEvent.Lumi;
      MpODiffractive.Event = MEvent.Event;
      MpODiffractive.hiHF_pf = MEvent.hiHF_pf;
      MpODiffractive.hiHFPlus_pf = MEvent.hiHFPlus_pf;
      MpODiffractive.hiHFMinus_pf = MEvent.hiHFMinus_pf;
      if (CollisionSystem == "OO" || CollisionSystem == "NeNe" ||
          CollisionSystem == "pO")
        MpODiffractive.hiBin = getHiBinFromhiHF(MEvent.hiHF_pf); // for pO this is not calibrated probably
      else
        MpODiffractive.hiBin = -9999;
      MpODiffractive.Ncoll = MEvent.Ncoll;
      MpODiffractive.Npart = MEvent.Npart;

      ////////////////////////////
      ////////// Vertex //////////
      ////////////////////////////

      int BestVertex = -1;
      for (int i = 0; i < MTrack.nVtx; i++) {
        if (BestVertex < 0 || MTrack.ptSumVtx->at(i) > MTrack.ptSumVtx->at(BestVertex))
          BestVertex = i;
      }
      if (BestVertex >= 0) {
        MpODiffractive.VX = MTrack.xVtx->at(BestVertex);
        MpODiffractive.VY = MTrack.yVtx->at(BestVertex);
        MpODiffractive.VZ = MTrack.zVtx->at(BestVertex);
        MpODiffractive.VXError = MTrack.xErrVtx->at(BestVertex);
        MpODiffractive.VYError = MTrack.yErrVtx->at(BestVertex);
        MpODiffractive.VZError = MTrack.zErrVtx->at(BestVertex);
        MpODiffractive.isFakeVtx = MTrack.isFakeVtx->at(BestVertex);
        MpODiffractive.ptSumVtx = MTrack.ptSumVtx->at(BestVertex);
        MpODiffractive.nTracksVtx = MTrack.nTracksVtx->at(BestVertex);
        MpODiffractive.chi2Vtx = MTrack.chi2Vtx->at(BestVertex);
        MpODiffractive.ndofVtx = MTrack.ndofVtx->at(BestVertex);
        MpODiffractive.bestVtxIndx = BestVertex;
      }
      MpODiffractive.nVtx = MTrack.nVtx;

      /////////////////////////////////////
      ////////// Event selection //////////
      /////////////////////////////////////
      MpODiffractive.ZDCsumPlus = MZDC.sumPlus;
      MpODiffractive.ZDCsumMinus = MZDC.sumMinus;
      MpODiffractive.ClusterCompatibilityFilter = MSkim.ClusterCompatibilityFilter;
      MpODiffractive.PVFilter = MSkim.PVFilter;
      MpODiffractive.mMaxL1HFAdcPlus = MHFAdc.mMaxL1HFAdcPlus;
      MpODiffractive.mMaxL1HFAdcMinus = MHFAdc.mMaxL1HFAdcMinus;
      MpODiffractive.hiHFPlus_pfle1 = MEvent.hiHFPlus_pfle1;
      MpODiffractive.hiHFPlus_pfle2 = MEvent.hiHFPlus_pfle2;
      MpODiffractive.hiHFPlus_pfle3 = MEvent.hiHFPlus_pfle3;
      MpODiffractive.hiHFMinus_pfle1 = MEvent.hiHFMinus_pfle1;
      MpODiffractive.hiHFMinus_pfle2 = MEvent.hiHFMinus_pfle2;
      MpODiffractive.hiHFMinus_pfle3 = MEvent.hiHFMinus_pfle3;

      if (CollisionSystem == "OO" || CollisionSystem == "NeNe") {
        if (IsData == true) {
          int HLT_OxyZeroBias_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyZeroBias_v1");
          int HLT_MinimumBiasHF_OR_BptxAND_v1 = MTrigger.CheckTriggerStartWith("HLT_MinimumBiasHF_OR_BptxAND_v1");
          MpODiffractive.HLT_OxySingleJet16_ZDC1nAsymXOR_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxySingleJet16_ZDC1nAsymXOR_v1");
          MpODiffractive.HLT_OxySingleJet16_ZDC1nXOR_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxySingleJet16_ZDC1nXOR_v1");
          MpODiffractive.HLT_OxySingleJet24_ZDC1nAsymXOR_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxySingleJet24_ZDC1nAsymXOR_v1");
          MpODiffractive.HLT_OxySingleJet24_ZDC1nXOR_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxySingleJet24_ZDC1nXOR_v1");
          MpODiffractive.HLT_OxyZDC1nOR_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyZDC1nOR_v1");
          MpODiffractive.HLT_OxyZeroBias_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyZeroBias_v1");
          MpODiffractive.HLT_MinimumBiasHF_OR_BptxAND_v1 =
              MTrigger.CheckTriggerStartWith("HLT_MinimumBiasHF_OR_BptxAND_v1");
          MpODiffractive.HLT_OxyL1SingleJet20_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyL1SingleJet20_v1");
          if (ApplyTriggerRejection == 1 && HLT_OxyZeroBias_v1 == 0) {
            continue;
          }
          if (ApplyTriggerRejection == 2 && HLT_MinimumBiasHF_OR_BptxAND_v1 == 0) {
            continue;
          }
        } // end of CollisionSystem == "OO" && IsData == true
      } // end of CollisionSystem == "OO"
      else if (CollisionSystem == "pO") {
        if (IsData == true) {
          MpODiffractive.HLT_OxyZeroBias_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyZeroBias_v1");
          MpODiffractive.HLT_OxyZDC1nOR_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyZDC1nOR_v1");
          MpODiffractive.HLT_OxySingleMuOpen_NotMBHF2OR_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxySingleMuOpen_NotMBHF2OR_v1");
          MpODiffractive.HLT_OxySingleJet8_ZDC1nAsymXOR_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxySingleJet8_ZDC1nAsymXOR_v1");
          MpODiffractive.HLT_OxyNotMBHF2_v1 = MTrigger.CheckTriggerStartWith("HLT_OxyNotMBHF2_v1");
          MpODiffractive.HLT_OxyZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxyZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v1");
          MpODiffractive.HLT_OxyZeroBias_MinPixelCluster400_v1 =
              MTrigger.CheckTriggerStartWith("HLT_OxyZeroBias_MinPixelCluster400_v1");
          MpODiffractive.HLT_MinimumBiasHF_OR_BptxAND_v1 =
              MTrigger.CheckTriggerStartWith("HLT_MinimumBiasHF_OR_BptxAND_v1");
          MpODiffractive.HLT_MinimumBiasHF_AND_BptxAND_v1 =
              MTrigger.CheckTriggerStartWith("HLT_MinimumBiasHF_AND_BptxAND_v1");
          if (ApplyTriggerRejection == 1 && MTrigger.CheckTriggerStartWith("HLT_OxyZeroBias_v1") == 0) {
            continue;
          }
          if (ApplyTriggerRejection == 2 && MTrigger.CheckTriggerStartWith("HLT_MinimumBiasHF_OR_BptxAND_v1") == 0) {
            continue;
          }
        } // end of CollisionSystem == "pO" && IsData == true
      } // end of CollisionSystem == "pO"

      ////////////////////////////////////////////////////
      //////// Fill Baseline evt. sel filters ////////////
      ////////////////////////////////////////////////////
      bool passBaselineEventSelection = false;
      if (CollisionSystem == "OO" || CollisionSystem == "NeNe" ||
                 CollisionSystem == "pO") { // PVfilter, PV position within 15 cm, ClusterCompatibilityFilter
        passBaselineEventSelection = getBaselinepOEventSel(MpODiffractive);
      }
      MpODiffractive.passBaselineEventSelection = passBaselineEventSelection;

      if (ApplyEventRejection && IsData == true) {
        if (passBaselineEventSelection == false) {
          continue;
        }
      }
      
      // loop over tracks
      int NTrack = MTrack.nTrk;
      MpODiffractive.nTrk = NTrack;
      int locMultiplicityEta2p4 = 0;

      for (int iTrack = 0; iTrack < NTrack; iTrack++) {
        if (MTrack.trkPt->at(iTrack) < 0.4 || abs(MTrack.trkEta->at(iTrack)) > 2.4)
          continue; // skip tracks with pT < 0.4 GeV/c and |eta| > 2.4

        if (MTrack.highPurity->at(iTrack)) {
          locMultiplicityEta2p4++;
        }

        float trkPt = MTrack.trkPt->at(iTrack);
        float trkEta = MTrack.trkEta->at(iTrack);
        float trkPhi = MTrack.trkPhi->at(iTrack);
        float trkPtError = MTrack.trkPtError->at(iTrack);
        bool highPurity = MTrack.highPurity->at(iTrack);
        float trkDxyAssociatedVtx = MTrack.trkDxyAssociatedVtx->at(iTrack);
        float trkDzAssociatedVtx = MTrack.trkDzAssociatedVtx->at(iTrack);
        float trkDxyErrAssociatedVtx = MTrack.trkDxyErrAssociatedVtx->at(iTrack);
        float trkDzErrAssociatedVtx = MTrack.trkDzErrAssociatedVtx->at(iTrack);
        int trkAssociatedVtxIndx = MTrack.trkAssociatedVtxIndx->at(iTrack);
        char trkCharge = MTrack.trkCharge->at(iTrack);
        char trkNHits = MTrack.trkNHits->at(iTrack);
        char trkNPixHits = MTrack.trkNPixHits->at(iTrack);
        char trkNLayers = MTrack.trkNLayers->at(iTrack);
        float trkNormChi2 = MTrack.trkNormChi2->at(iTrack);
        float pfEnergy = MTrack.pfEnergy->at(iTrack);

        MpODiffractive.trkEta->push_back(trkEta);
        MpODiffractive.trkPt->push_back(trkPt);
        MpODiffractive.trkPhi->push_back(trkPhi);
        MpODiffractive.trkPtError->push_back(trkPtError);
        MpODiffractive.highPurity->push_back(highPurity);
        MpODiffractive.trkDxyAssociatedVtx->push_back(trkDxyAssociatedVtx);
        MpODiffractive.trkDzAssociatedVtx->push_back(trkDzAssociatedVtx);
        MpODiffractive.trkDxyErrAssociatedVtx->push_back(trkDxyErrAssociatedVtx);
        MpODiffractive.trkDzErrAssociatedVtx->push_back(trkDzErrAssociatedVtx);
        MpODiffractive.trkAssociatedVtxIndx->push_back(trkAssociatedVtxIndx);
        MpODiffractive.trkCharge->push_back(trkCharge);
        MpODiffractive.trkNHits->push_back(trkNHits);
        MpODiffractive.trkNPixHits->push_back(trkNPixHits);
        MpODiffractive.trkNLayers->push_back(trkNLayers);
        MpODiffractive.trkNormChi2->push_back(trkNormChi2);
        MpODiffractive.pfEnergy->push_back(pfEnergy);
        // not all forests have dE/dx info
        if (MTrack.dedxAllLikelihood) MpODiffractive.dedxAllLikelihood->push_back(MTrack.dedxAllLikelihood->at(iTrack));
        if (MTrack.dedxPixelLikelihood) MpODiffractive.dedxPixelLikelihood->push_back(MTrack.dedxPixelLikelihood->at(iTrack));
        if (MTrack.dedxStripLikelihood) MpODiffractive.dedxStripLikelihood->push_back(MTrack.dedxStripLikelihood->at(iTrack));

      } // end of loop over tracks (gen or reco)

      MpODiffractive.multiplicityEta2p4 = locMultiplicityEta2p4;

      ////////////////////////////
      ////// PF variables ////////
      ////////////////////////////

      int N_pf = MPF.PT->size();
      for (int iPF = 0; iPF < N_pf; iPF++) {
        MpODiffractive.pfID->push_back(MPF.ID->at(iPF));
        MpODiffractive.pfPT->push_back(MPF.PT->at(iPF));
        MpODiffractive.pfE->push_back(MPF.E->at(iPF));
        MpODiffractive.pfEta->push_back(MPF.Eta->at(iPF));
        MpODiffractive.pfPhi->push_back(MPF.Phi->at(iPF));
        MpODiffractive.pfM->push_back(MPF.M->at(iPF));
      }

      ////////////////////////////
      ///// Debug variables //////
      ////////////////////////////

      if (DebugMode) {
        for (int iDebVtx = 0; iDebVtx < MTrack.nVtx; iDebVtx++) {
          MpODiffractive.AllxVtx->push_back(MTrack.xVtx->at(iDebVtx));
          MpODiffractive.AllyVtx->push_back(MTrack.yVtx->at(iDebVtx));
          MpODiffractive.AllzVtx->push_back(MTrack.zVtx->at(iDebVtx));
          MpODiffractive.AllxVtxError->push_back(MTrack.xErrVtx->at(iDebVtx));
          MpODiffractive.AllyVtxError->push_back(MTrack.yErrVtx->at(iDebVtx));
          MpODiffractive.AllzVtxError->push_back(MTrack.zErrVtx->at(iDebVtx));
          MpODiffractive.AllisFakeVtx->push_back(MTrack.isFakeVtx->at(iDebVtx));
          MpODiffractive.AllnTracksVtx->push_back(MTrack.nTracksVtx->at(iDebVtx));
          MpODiffractive.Allchi2Vtx->push_back(MTrack.chi2Vtx->at(iDebVtx));
          MpODiffractive.AllndofVtx->push_back(MTrack.ndofVtx->at(iDebVtx));
          MpODiffractive.AllptSumVtx->push_back(MTrack.ptSumVtx->at(iDebVtx));
        }
      }

      ////////////////////////////
      /// PPS & FSC variables ////
      ////////////////////////////

      if (includeFSCandPPSMode) {
        // PPS variables
        if (MPPS.n > PPSMAXN) {
          std::cout << "ERROR: in the PPS tree of the forest n > PPSMAXN; skipping PPS information filling"
                    << std::endl;
        } else {
          for (int iPPS = 0; iPPS < MPPS.n; iPPS++) {
            fillPPSInfo(MpODiffractive, MPPS, iPPS);
          }
        }

        // FSC variables
        if (MFSC.n > FSCMAXN) {
          std::cout << "ERROR: in the FSC tree of the forest n > FSCMAXN; skipping FSC information filling"
                    << std::endl;
        } else {
          for (int iFSC = 0; iFSC < MFSC.n; iFSC++) {
            fillFSCInfo(MpODiffractive, MFSC, iFSC);
          }
        }
      }

      MpODiffractive.FillEntry();
    }
    if (!HideProgressBar) {
      Bar.Update(EntryCount);
      Bar.Print();
      Bar.PrintLine();
    } // if (!HideProgressBar)
    InputFile.Close();
  }
  OutputFile.cd();
  Tree.Write();
  InfoTree.Write();

  OutputFile.Close();

  return 0;
}
