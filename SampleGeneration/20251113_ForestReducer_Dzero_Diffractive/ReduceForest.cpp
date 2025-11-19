// ===================================================
// Script Based on Yi Chen Reducer script
// Yi Chen (Vanderbilt) Yi.Chen@cern.ch
// https://github.com/FHead/PhysicsZHadronEEC
// ===================================================

#include <iostream>
#include <set>
#include <algorithm>
#include <vector>
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

#include "include/DmesonSelection.h"
#include "include/PIDScoring.h"
#include "include/parseFSCInfo.h"

const std::vector<std::string> triggerPaths2023 = {
  "HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v",
  "HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v",
  "HLT_HIUPC_SingleJet8_ZDC1nXOR_MaxPixelCluster50000_v",
  "HLT_HIUPC_SingleJet8_ZDC1nAsymXOR_MaxPixelCluster50000_v",
  "HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v",
  "HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v",
  "HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v"
};

const std::vector<std::string> triggerPaths2025 = {
  "HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v",
  "HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v",
  "HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v",
  "HLT_HIUPC_ZeroBias_MaxPixelCluster10000_v",
  "HLT_HIUPC_SingleJet8_ZDC1nXOR_MaxPixelCluster10000_v",
  "HLT_HIUPC_SingleJet8_ZDC1nAsymXOR_MaxPixelCluster10000_v",
  "HLT_HIUPC_ZDC1nOR_MaxPixelCluster10000_v",
  "HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v",
  "HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v",
  "HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v",
  "HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v",
  "HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v",
  "HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster1000_v"
};

bool logical_or_vectBool(std::vector<bool>* vec) {
    return std::any_of(vec->begin(), vec->end(), [](bool b) { return b; });
}

// Helper function to convert a string to lowercase
std::string toLower(const std::string& str) {
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    return lowerStr;
}

int main(int argc, char *argv[]);
double GetMaxEnergyHF(PFTreeMessenger *M, double etaMin, double etaMax);

int main(int argc, char *argv[]) {
  string VersionString = "V8";

  CommandLine CL(argc, argv);

  vector<string> InputFileNames = CL.GetStringVector("Input");
  string OutputFileName = CL.Get("Output");

  // bool DoGenLevel                    = CL.GetBool("DoGenLevel", true);
  bool IsData = CL.GetBool("IsData", false);
  bool IsGammaNMCtype = CL.GetBool("IsGammaNMCtype", true); // This is only meaningful when IsData==false. gammaN: BeamA, Ngamma: BeamB
  int Year = CL.GetInt("Year", 2023);
  bool DoPID = CL.GetBool("DoPID", true);

  double Fraction = CL.GetDouble("Fraction", 1.00);
  float ZDCMinus1nThreshold = CL.GetDouble("ZDCMinus1nThreshold", 1000.);
  float ZDCPlus1nThreshold = CL.GetDouble("ZDCPlus1nThreshold", 1100.);
  int ApplyTriggerRejection = CL.GetInteger("ApplyTriggerRejection", 0);
  bool ApplyEventRejection = CL.GetBool("ApplyEventRejection", false);
  bool ApplyZDCGapRejection = CL.GetBool("ApplyZDCGapRejection", false);

  // options for how to reject non-selected D candidates (case-insensitive)
  // "NO":              keep all D's
  // "OR":              keep D's that pass the OR logic of all D selections
  // "Nominal":         reject !DpassCutNominal
  // "Loose":           reject !DpassCutLoose
  // "SystDsvpvSig":    reject !DpassCutSystDsvpvSig
  // "SystDtrkPt":      reject !DpassCutSystDtrkPt
  // "SystDalpha":      reject !DpassCutSystDalpha
  // "SystDdtheta":     reject !DpassCutSystDdtheta
  // "SystDchi2cl":     reject !DpassCutSystDchi2cl
  // "PASOR":           keep D's that pass the OR logic of all PAS selections
  //                    (this is a tighter selection than "OR")
  // "PAS":             reject !DpassCut23PAS
  // "PASSystDsvpvSig": reject !DpassCut23PASSystDsvpvSig
  // "PASSystDtrkPt":   reject !DpassCut23PASSystDtrkPt
  // "PASSystDalpha":   reject !DpassCut23PASSystDalpha
  // "PASSystDchi2cl":  reject !DpassCut23PASSystDchi2cl
  string ApplyDRejection = toLower(CL.Get("ApplyDRejection", "no"));
  string TrackTreeName = CL.Get("TrackTree", "PbPbTracks/trackTree");
  string PFTreeName = CL.Get("PFTree", "particleFlowAnalyser/pftree");
  string ZDCTreeName = CL.Get("ZDCTree", "zdcanalyzer/zdcrechit");
  string FSCTreeName = CL.Get("FSCTree", "fscanalyzer/fscdigi");
  bool HideProgressBar = CL.GetBool("HideProgressBar", false);
  TFile OutputFile(OutputFileName.c_str(), "RECREATE");
  TTree Tree("Tree", Form("Tree for UPC Dzero analysis (%s)", VersionString.c_str()));
  TTree InfoTree("InfoTree", "Information");
  DzeroDiffractivePbPbTreeMessenger MDzeroDiffractive;
  MDzeroDiffractive.SetBranch(&Tree);

  // Allowed ApplyDRejection in lowercase
  std::set<std::string> allowedApplyDRejection = {
    "no", "or", "nominal", "loose",
    "systdsvpvsig", "systdtrkpt", "systdalpha", "systddtheta", "systdchi2cl",
    "pasor", "pas",
    "passystdsvpvsig", "passystdtrkpt", "passystdalpha", "passystdchi2cl"
  };
  // Validate the argument
  if (allowedApplyDRejection.find(ApplyDRejection) != allowedApplyDRejection.end()) {
    std::cout << "D filtering criterion: " << ApplyDRejection << std::endl;
    // Add your program logic here
  } else {
    std::cerr << "[Warning] Invalid ApplyDRejection. Set to NO rejection. Allowed ApplyDRejection are: ";
    for (const auto& ele : allowedApplyDRejection) {
        std::cout << ele << ", ";
    }
    std::cout << std::endl;
    ApplyDRejection = "no";
  }

  for (string InputFileName : InputFileNames) {
    TFile InputFile(InputFileNames[0].c_str());

    HiEventTreeMessenger MEvent(InputFile); // hiEvtAnalyzer/HiTree
    PbPbUPCTrackTreeMessenger MTrackPbPbUPC(InputFile, TrackTreeName); // PbPbTracks/trackTree
    GenParticleTreeMessenger MGen(InputFile); // HiGenParticleAna/hi
    PFTreeMessenger MPF(InputFile, PFTreeName); // particleFlowAnalyser/pftree
    SkimTreeMessenger MSkim(InputFile); // skimanalysis/HltTree
    TriggerTreeMessenger MTrigger(InputFile); // hltanalysis/HltTree
    DzeroTreeMessenger MDzero(InputFile); // Dfinder/ntDkpi
    DzeroGenTreeMessenger MDzeroGen(InputFile); // Dfinder/ntGen
    ZDCTreeMessenger MZDC(InputFile, ZDCTreeName); // zdcanalyzer/zdcrechit
    METFilterTreeMessenger MMETFilter(InputFile); // l1MetFilterRecoTree/MetFilterRecoTree
    FSCTreeMessenger MFSC(InputFile, FSCTreeName); // fscanalyzer/fscdigi

    int EntryCount = MEvent.GetEntries() * Fraction;
    ProgressBar Bar(cout, EntryCount);
    if (!HideProgressBar) {
      Bar.SetStyle(-1);
    }
    
    vector<TF1*> dedxFunctions = ImportPIDRoot("../../CommonCode/root/DzeroUPC_dedxMap.root");
    TF1* fdedxPionCenter  = dedxFunctions[0];
    TF1* fdedxPionSigmaLo = dedxFunctions[1];
    TF1* fdedxPionSigmaHi = dedxFunctions[2];
    TF1* fdedxKaonCenter  = dedxFunctions[3];
    TF1* fdedxKaonSigmaLo = dedxFunctions[4];
    TF1* fdedxKaonSigmaHi = dedxFunctions[5];
    TF1* fdedxProtCenter  = dedxFunctions[6];
    TF1* fdedxProtSigmaLo = dedxFunctions[7];
    TF1* fdedxProtSigmaHi = dedxFunctions[8];
    
    /////////////////////////////////
    //////// Main Event Loop ////////
    /////////////////////////////////
    for (int iE = 0; iE < EntryCount; iE++) {
      if (!HideProgressBar &&
        (EntryCount < 300 || (iE % (EntryCount / 250)) == 0)
      ) {
        Bar.Update(iE);
        Bar.Print();
      }

      MEvent.GetEntry(iE);
      MGen.GetEntry(iE);
      MTrackPbPbUPC.GetEntry(iE);
      MPF.GetEntry(iE);
      MSkim.GetEntry(iE);
      MTrigger.GetEntry(iE);
      MDzero.GetEntry(iE);
      if (IsData == false)
        MDzeroGen.GetEntry(iE);
      MZDC.GetEntry(iE);
      MDzeroDiffractive.Clear();
      MMETFilter.GetEntry(iE);
      if (MFSC.Tree != nullptr) MFSC.GetEntry(iE);
      
      ////////////////////////////////////////
      ////////// Global event stuff //////////
      ////////////////////////////////////////

      MDzeroDiffractive.Run = MEvent.Run;
      MDzeroDiffractive.Lumi = MEvent.Lumi;
      MDzeroDiffractive.Event = MEvent.Event;

      ////////////////////////////
      ////////// Vertex //////////
      ////////////////////////////

      int BestVertex = -1;

      for (int i = 0; i < MTrackPbPbUPC.nVtx; i++) {
        if (BestVertex < 0 || MTrackPbPbUPC.ptSumVtx->at(i) > MTrackPbPbUPC.ptSumVtx->at(BestVertex))
          BestVertex = i;
      }
      if (BestVertex >= 0) {
        MDzeroDiffractive.VX = MTrackPbPbUPC.xVtx->at(BestVertex);
        MDzeroDiffractive.VY = MTrackPbPbUPC.yVtx->at(BestVertex);
        MDzeroDiffractive.VZ = MTrackPbPbUPC.zVtx->at(BestVertex);
        MDzeroDiffractive.VXError = MTrackPbPbUPC.xErrVtx->at(BestVertex);
        MDzeroDiffractive.VYError = MTrackPbPbUPC.yErrVtx->at(BestVertex);
        MDzeroDiffractive.VZError = MTrackPbPbUPC.zErrVtx->at(BestVertex);
      }
      MDzeroDiffractive.nVtx = MTrackPbPbUPC.nVtx;
      
      /////////////////////////////////////
      ////////// Event selection //////////
      /////////////////////////////////////
      if (IsData == false) {
        MDzeroDiffractive.Gsize = MDzeroGen.Gsize;
        for (int iDGen = 0; iDGen < MDzeroGen.Gsize; iDGen++) {
          MDzeroDiffractive.Gpt->push_back(MDzeroGen.Gpt[iDGen]);
          MDzeroDiffractive.Gy->push_back(MDzeroGen.Gy[iDGen]);
          bool isSignalGen =
              (MDzeroGen.GisSignal[iDGen] == 1 || MDzeroGen.GisSignal[iDGen] == 2) && MDzeroGen.Gpt[iDGen] > 0.;
          bool isPromptGen = MDzeroGen.GBAncestorpdgId[iDGen] == 0;
          bool isFeeddownGen = (MDzeroGen.GBAncestorpdgId[iDGen] >= 500 && MDzeroGen.GBAncestorpdgId[iDGen] < 600) ||
                                (MDzeroGen.GBAncestorpdgId[iDGen] > -600 && MDzeroGen.GBAncestorpdgId[iDGen] <= -500);
          MDzeroDiffractive.GisSignalCalc->push_back(isSignalGen);
          MDzeroDiffractive.GisSignalCalcPrompt->push_back(isSignalGen && isPromptGen);
          MDzeroDiffractive.GisSignalCalcFeeddown->push_back(isSignalGen && isFeeddownGen);
        }
      }
      if (IsData == true) {
        std::vector<std::string> triggerPaths;
        if (Year == 2023) {
          triggerPaths = triggerPaths2023;
        } else if (Year == 2025) {
          triggerPaths = triggerPaths2025;
        }
        for (int iPath = 0; iPath < triggerPaths.size(); iPath++) {
          int triggerResult = MTrigger.CheckTriggerStartWith(triggerPaths[iPath]);
          if (triggerResult != -1) {
            MDzeroDiffractive.triggerNames->push_back(triggerPaths[iPath]);
            MDzeroDiffractive.triggerBits->push_back(triggerResult);
          }
        }
      }
      
      if (IsData == true) {
        MDzeroDiffractive.ZDCsumPlus = MZDC.sumPlus;
        MDzeroDiffractive.ZDCsumMinus = MZDC.sumMinus;
        bool selectedBkgFilter = MSkim.ClusterCompatibilityFilter == 1 && MMETFilter.cscTightHalo2015Filter;
        bool selectedVtxFilter = (MSkim.PVFilter == 1 && MTrackPbPbUPC.zVtx != nullptr && MTrackPbPbUPC.zVtx->size() > 0 && fabs(MTrackPbPbUPC.zVtx->at(0)) < 15.);
        if (ApplyEventRejection && IsData && (selectedBkgFilter == false || selectedVtxFilter == false)) continue;
        MDzeroDiffractive.selectedBkgFilter = selectedBkgFilter;
        MDzeroDiffractive.selectedVtxFilter = selectedVtxFilter;
        bool ZDCgammagamma = (MZDC.sumMinus < ZDCMinus1nThreshold && MZDC.sumPlus < ZDCPlus1nThreshold);
        bool ZDCgammaN = (MZDC.sumMinus > ZDCMinus1nThreshold && MZDC.sumPlus < ZDCPlus1nThreshold);
        bool ZDCNgamma = (MZDC.sumMinus < ZDCMinus1nThreshold && MZDC.sumPlus > ZDCPlus1nThreshold);
        MDzeroDiffractive.ZDCgammagamma = ZDCgammagamma;
        MDzeroDiffractive.ZDCgammaN = ZDCgammaN;
        MDzeroDiffractive.ZDCNgamma = ZDCNgamma;
        if (!ZDCgammagamma && !ZDCgammaN && !ZDCNgamma) {
          continue; // filter out events with both ZDC sides having neutrons
        }
      } // end of if (IsData == true)
      else { // if (IsData == false)
        // MDzeroDiffractive.ZDCsumPlus = MZDC.sumPlus;
        // MDzeroDiffractive.ZDCsumMinus = MZDC.sumMinus;
        bool selectedBkgFilter = MSkim.ClusterCompatibilityFilter == 1; // METFilter always true for MC
        bool selectedVtxFilter = (MSkim.PVFilter == 1 && MTrackPbPbUPC.zVtx != nullptr && MTrackPbPbUPC.zVtx->size() > 0 && fabs(MTrackPbPbUPC.zVtx->at(0)) < 15.);
        MDzeroDiffractive.selectedBkgFilter = selectedBkgFilter;
        MDzeroDiffractive.selectedVtxFilter = selectedVtxFilter;
        bool ZDCgammaN =  IsGammaNMCtype;
        bool ZDCNgamma = !IsGammaNMCtype;
        MDzeroDiffractive.ZDCgammaN = ZDCgammaN;
        MDzeroDiffractive.ZDCNgamma = ZDCNgamma;
      } // end of if (IsData == false)
      
      // Loop through the specified ranges for gapgammaN and gapNgamma
      // gammaN[4] and Ngamma[4] are nominal selection criteria
      //float EMaxHFPlus = GetMaxEnergyHF(&MPF, 3., 5.2);
      //float EMaxHFMinus = GetMaxEnergyHF(&MPF, -5.2, -3.);
      
      float EMaxHFPlus = MEvent.hiHFPlus_pfle1;
      float EMaxHFMinus = MEvent.hiHFMinus_pfle1;
      MDzeroDiffractive.HFEMaxPlus = EMaxHFPlus;
      MDzeroDiffractive.HFEMaxMinus = EMaxHFMinus;
      bool gapgammaN = (EMaxHFPlus < 9.2 && EMaxHFMinus > 8.6);
      bool gapNgamma = (EMaxHFMinus < 8.6 && EMaxHFPlus > 9.2);
      bool gapgammagamma = (EMaxHFMinus < 8.6 && EMaxHFPlus < 9.2);

      MDzeroDiffractive.gapgammaN = gapgammaN;
      MDzeroDiffractive.gapNgamma = gapNgamma;
      MDzeroDiffractive.gapgammagamma = gapgammagamma;
      
      int nTrackInAcceptanceHP = 0;
      for (int iTrack = 0; iTrack < MTrackPbPbUPC.nTrk; iTrack++) {
        if (MTrackPbPbUPC.trkPt->at(iTrack) <= 0.2)
          continue;
        if (fabs(MTrackPbPbUPC.trkEta->at(iTrack)) >= 2.4)
          continue;
        if (MTrackPbPbUPC.highPurity->at(iTrack) == false)
          continue;
        nTrackInAcceptanceHP++;
      }
      MDzeroDiffractive.nTrackInAcceptanceHP = nTrackInAcceptanceHP;

      // ----------------- FSC information ---------------------------

      if (MFSC.Tree != nullptr) {
        // FSC variables
        if (MFSC.n > FSCMAXN) {
          std::cout << "ERROR: in the FSC tree of the forest n > FSCMAXN; skipping FSC information filling"
                    << std::endl;
        } else {
          for (int iFSC = 0; iFSC < MFSC.n; iFSC++) {
            fillFSCInfo(MDzeroDiffractive, MFSC, iFSC);
          }
        }
      }

      // -------------------------------------------------------------
      
      if (MDzero.Dsize == 20000) cout << "[WARNING] Event " << iE << " has max Dsize!" << endl;
      int countSelDzero = 0;
      for (int iD = 0; iD < MDzero.Dsize; iD++) {
        bool DpassCutNominal_           = DpassCutNominal(MDzero, iD);
        bool DpassCutSystDsvpvSig_      = DpassCutSystDsvpvSig(MDzero, iD);
        bool DpassCutSystDtrkPt_        = DpassCutSystDtrkPt(MDzero, iD);
        bool DpassCutSystDalpha_        = DpassCutSystDalpha(MDzero, iD);
        bool DpassCutSystDalphaDdtheta_ = DpassCutSystDalphaDdtheta(MDzero, iD);
        bool DpassCutSystDdtheta_       = DpassCutSystDdtheta(MDzero, iD);
        bool DpassCutSystDchi2cl_       = DpassCutSystDchi2cl(MDzero, iD);
        bool DpassCutLoose_             = DpassCutLoose(MDzero, iD);
        bool DpassCut23PAS_             = DpassCut23PAS(MDzero, iD);
        bool DpassCut23LowPt_           = DpassCut23LowPt(MDzero, iD);
        bool DpassCut23PASSystDsvpvSig_ = DpassCut23PASSystDsvpvSig(MDzero, iD);
        bool DpassCut23PASSystDtrkPt_   = DpassCut23PASSystDtrkPt(MDzero, iD);
        bool DpassCut23PASSystDalpha_   = DpassCut23PASSystDalpha(MDzero, iD);
        bool DpassCut23PASSystDchi2cl_  = DpassCut23PASSystDchi2cl(MDzero, iD);
        if (IsData)
        {
          if (ApplyDRejection=="or") { // OR for new cuts only
            if (!DpassCutNominal_           &&
                !DpassCutSystDsvpvSig_      &&
                !DpassCutSystDtrkPt_        &&
                !DpassCutSystDalpha_        &&
                !DpassCutSystDdtheta_       &&
                !DpassCutSystDalphaDdtheta_ &&
                !DpassCutSystDchi2cl_       &&
                !DpassCutLoose_ ) continue;
          }
          else if (ApplyDRejection=="nominal"         && !DpassCutNominal_) continue;
          else if (ApplyDRejection=="systdsvpvsig"    && !DpassCutSystDsvpvSig_) continue;
          else if (ApplyDRejection=="systdtrkpt"      && !DpassCutSystDtrkPt_) continue;
          else if (ApplyDRejection=="systdalpha"      && !DpassCutSystDalpha_) continue;
          else if (ApplyDRejection=="systddtheta"     && !DpassCutSystDdtheta_) continue;
          else if (ApplyDRejection=="systdchi2cl"     && !DpassCutSystDchi2cl_) continue;
          else if (ApplyDRejection=="loose"           && !DpassCutLoose_) continue;
          else if (ApplyDRejection=="pasor") { // OR for PAS23 cuts only
            if (!DpassCut23PAS_ &&
                !DpassCut23PASSystDsvpvSig_ &&
                !DpassCut23PASSystDtrkPt_ &&
                !DpassCut23PASSystDalpha_ &&
                !DpassCut23PASSystDchi2cl_ ) continue;
          }
          else if (ApplyDRejection=="pas"             && !DpassCut23PAS_) continue;
          else if (ApplyDRejection=="passystdsvpvsig" && !DpassCut23PASSystDsvpvSig_) continue;
          else if (ApplyDRejection=="passystdtrkpt"   && !DpassCut23PASSystDtrkPt_) continue;
          else if (ApplyDRejection=="passystdalpha"   && !DpassCut23PASSystDalpha_) continue;
          else if (ApplyDRejection=="passystdchi2cl"  && !DpassCut23PASSystDchi2cl_) continue;
        }
        countSelDzero++;
        MDzeroDiffractive.Dpt->push_back(             MDzero.Dpt[iD]);
        MDzeroDiffractive.Dy->push_back(              MDzero.Dy[iD]);
        MDzeroDiffractive.Dmass->push_back(           MDzero.Dmass[iD]);
        MDzeroDiffractive.Dchi2cl->push_back(         MDzero.Dchi2cl[iD]);
        MDzeroDiffractive.DsvpvDistance->push_back(   MDzero.DsvpvDistance[iD]);
        MDzeroDiffractive.DsvpvDisErr->push_back(     MDzero.DsvpvDisErr[iD]);
        MDzeroDiffractive.DsvpvDistance_2D->push_back(MDzero.DsvpvDistance_2D[iD]);
        MDzeroDiffractive.DsvpvDisErr_2D->push_back(  MDzero.DsvpvDisErr_2D[iD]);
  //        MDzeroDiffractive.Dip3d->push_back(           MDzero.Dip3d[iD]);
  //        MDzeroDiffractive.Dip3derr->push_back(        MDzero.Dip3derr[iD]);
        MDzeroDiffractive.Dalpha->push_back(          MDzero.Dalpha[iD]);
        MDzeroDiffractive.Ddtheta->push_back(         MDzero.Ddtheta[iD]);
        
        MDzeroDiffractive.Dtrk1P->push_back(          MDzero.Dtrk1P[iD]);
        MDzeroDiffractive.Dtrk1Pt->push_back(         MDzero.Dtrk1Pt[iD]);
        MDzeroDiffractive.Dtrk1PtErr->push_back(      MDzero.Dtrk1PtErr[iD]);
        MDzeroDiffractive.Dtrk1Eta->push_back(        MDzero.Dtrk1Eta[iD]);
        MDzeroDiffractive.Dtrk1dedx->push_back(       MDzero.Dtrk1dedx[iD]);
        MDzeroDiffractive.Dtrk1MassHypo->push_back(   MDzero.Dtrk1MassHypo[iD]);
        MDzeroDiffractive.Dtrk1PixelHit->push_back(   MDzero.Dtrk1PixelHit[iD]);
        MDzeroDiffractive.Dtrk1StripHit->push_back(   MDzero.Dtrk1StripHit[iD]);
        if (DoPID && MDzero.Dtrk1P[iD] < 2.5) // Only give valid PID for p_track < 2 GeV
        {
          MDzeroDiffractive.Dtrk1PionScore->push_back(GetPIDScore(
            MDzero.Dtrk1P[iD], MDzero.Dtrk1dedx[iD],
            fdedxPionCenter, fdedxPionSigmaLo, fdedxPionSigmaHi));
          MDzeroDiffractive.Dtrk1KaonScore->push_back(GetPIDScore(
            MDzero.Dtrk1P[iD], MDzero.Dtrk1dedx[iD],
            fdedxKaonCenter, fdedxKaonSigmaLo, fdedxKaonSigmaHi));
          MDzeroDiffractive.Dtrk1ProtScore->push_back(GetPIDScore(
            MDzero.Dtrk1P[iD], MDzero.Dtrk1dedx[iD],
            fdedxProtCenter, fdedxProtSigmaLo, fdedxProtSigmaHi));
        } else {
          MDzeroDiffractive.Dtrk1PionScore->push_back(-999.);
          MDzeroDiffractive.Dtrk1KaonScore->push_back(-999.);
          MDzeroDiffractive.Dtrk1ProtScore->push_back(-999.);
        }
        
        MDzeroDiffractive.Dtrk2P->push_back(          MDzero.Dtrk2P[iD]);
        MDzeroDiffractive.Dtrk2Pt->push_back(         MDzero.Dtrk2Pt[iD]);
        MDzeroDiffractive.Dtrk2PtErr->push_back(      MDzero.Dtrk2PtErr[iD]);
        MDzeroDiffractive.Dtrk2Eta->push_back(        MDzero.Dtrk2Eta[iD]);
        MDzeroDiffractive.Dtrk2dedx->push_back(       MDzero.Dtrk2dedx[iD]);
        MDzeroDiffractive.Dtrk2MassHypo->push_back(   MDzero.Dtrk2MassHypo[iD]);
        MDzeroDiffractive.Dtrk2PixelHit->push_back(   MDzero.Dtrk2PixelHit[iD]);
        MDzeroDiffractive.Dtrk2StripHit->push_back(   MDzero.Dtrk2StripHit[iD]);
        if (DoPID && MDzero.Dtrk2P[iD] < 2.5) // Only do PID for p_track < 2 GeV
        {
          MDzeroDiffractive.Dtrk2PionScore->push_back(GetPIDScore(
            MDzero.Dtrk2P[iD], MDzero.Dtrk2dedx[iD],
            fdedxPionCenter, fdedxPionSigmaLo, fdedxPionSigmaHi));
          MDzeroDiffractive.Dtrk2KaonScore->push_back(GetPIDScore(
            MDzero.Dtrk2P[iD], MDzero.Dtrk2dedx[iD],
            fdedxKaonCenter, fdedxKaonSigmaLo, fdedxKaonSigmaHi));
          MDzeroDiffractive.Dtrk2ProtScore->push_back(GetPIDScore(
            MDzero.Dtrk2P[iD], MDzero.Dtrk2dedx[iD],
            fdedxProtCenter, fdedxProtSigmaLo, fdedxProtSigmaHi));
        } else {
          MDzeroDiffractive.Dtrk2PionScore->push_back(-999.);
          MDzeroDiffractive.Dtrk2KaonScore->push_back(-999.);
          MDzeroDiffractive.Dtrk2ProtScore->push_back(-999.);
        }
        
        MDzeroDiffractive.DpassCutNominal->push_back(           DpassCutNominal_);
        MDzeroDiffractive.DpassCutSystDsvpvSig->push_back(      DpassCutSystDsvpvSig_);
        MDzeroDiffractive.DpassCutSystDtrkPt->push_back(        DpassCutSystDtrkPt_);
        MDzeroDiffractive.DpassCutSystDalpha->push_back(        DpassCutSystDalpha_);
        MDzeroDiffractive.DpassCutSystDdtheta->push_back(       DpassCutSystDdtheta_);
        MDzeroDiffractive.DpassCutSystDalphaDdtheta->push_back( DpassCutSystDalphaDdtheta_);
        MDzeroDiffractive.DpassCutSystDchi2cl->push_back(       DpassCutSystDchi2cl_);
        MDzeroDiffractive.DpassCutLoose->push_back(             DpassCutLoose_);
        MDzeroDiffractive.DpassCut23PAS->push_back(             DpassCut23PAS_);
        MDzeroDiffractive.DpassCut23PASSystDsvpvSig->push_back( DpassCut23PASSystDsvpvSig_);
        MDzeroDiffractive.DpassCut23PASSystDtrkPt->push_back(   DpassCut23PASSystDtrkPt_);
        MDzeroDiffractive.DpassCut23PASSystDalpha->push_back(   DpassCut23PASSystDalpha_);
        MDzeroDiffractive.DpassCut23PASSystDchi2cl->push_back(  DpassCut23PASSystDchi2cl_);
        if (IsData == false) {
          MDzeroDiffractive.Dgen->push_back(MDzero.Dgen[iD]);
          bool isSignalGenMatched = (MDzero.Dgen[iD] == 23333 || MDzero.Dgen[iD] == 41022 || MDzero.Dgen[iD] == 41044) && MDzero.Dgenpt[iD] > 0.;
          bool isPromptGenMatched = MDzero.DgenBAncestorpdgId[iD] == 0;
          bool isFeeddownGenMatched = (MDzero.DgenBAncestorpdgId[iD] >= 500 && MDzero.DgenBAncestorpdgId[iD] < 600) ||
                                      (MDzero.DgenBAncestorpdgId[iD] > -600 && MDzero.DgenBAncestorpdgId[iD] <= -500);
          MDzeroDiffractive.DisSignalCalc->push_back(isSignalGenMatched);
          MDzeroDiffractive.DisSignalCalcPrompt->push_back(isSignalGenMatched && isPromptGenMatched);
          MDzeroDiffractive.DisSignalCalcFeeddown->push_back(isSignalGenMatched && isFeeddownGenMatched);
        }
      }
      MDzeroDiffractive.Dsize = countSelDzero;
      MDzeroDiffractive.FillEntry();
    }
    if (!HideProgressBar) {
      Bar.Update(EntryCount);
      Bar.Print();
      Bar.PrintLine();
    }

    InputFile.Close();
  }
  
  OutputFile.cd();
  Tree.Write();
  InfoTree.Write();

  OutputFile.Close();

  return 0;
}

// ============================================================================ //
// Function to Retrieve Maximum Energy in HF Region within Specified Eta Range
// ============================================================================ //
double GetMaxEnergyHF(PFTreeMessenger *M, double etaMin = 3., double etaMax = 5.) {
  if (M == nullptr)
    return -1;
  if (M->Tree == nullptr)
    return -1;

  double EMax = 0;
  for (int iPF = 0; iPF < M->ID->size(); iPF++) {
    if ((M->ID->at(iPF) == 6 || M->ID->at(iPF) == 7) && M->Eta->at(iPF) > etaMin && M->Eta->at(iPF) < etaMax) {
      if (M->E->at(iPF) > EMax)
        EMax = M->E->at(iPF);
    }
  }
  return EMax;
}
