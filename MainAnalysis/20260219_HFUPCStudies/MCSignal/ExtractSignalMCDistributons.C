#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>
#include <string>

using namespace std;
#include "CommandLine.h" // Yi's Commandline bundle
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "utilities.h"   // Yen-Jie's random utility functions
#include "InfoManager.h"

#define DMASS 1.86484
#define DMASSMIN 1.66
#define DMASSMAX 2.26
#define DMASSNBINS 48

#define CheckAndSetBranch(TREE, BRANCH)                                 \
  if (TREE->GetBranch(#BRANCH)) {                                       \
    TREE->SetBranchStatus(#BRANCH, 1);                                  \
    TREE->SetBranchAddress(#BRANCH, &BRANCH);                           \
  } else {                                                              \
    std::cout<<__FUNCTION__<<": " #BRANCH " is not found in the input tree."<<std::endl; \
  }

bool isResolved(int ProcessID) {return ProcessID == 113; }; // 113 is resolved photon process ID
bool isDirect(int ProcessID) {return (ProcessID == 281 || ProcessID == 282 || ProcessID == 283 || ProcessID == 284); }; // 281-284 is direct photon process ID

bool eventSelection(DzeroUPCTreeMessenger *b, int IsResolved, int useClusterCompatibilityFilter) {  
  if (IsResolved == 1 && isResolved(b->ProcessID) == false)
    return false;
  else if (IsResolved == 0 && isDirect(b->ProcessID) == false)
    return false;

  if ((useClusterCompatibilityFilter == 1 && b->ClusterCompatibilityFilter == false) || b->selectedVtxFilter == false)
    return false;

  if (b->nVtx >= 3) return false;
  return true;
}

class DataAnalyzer {
public:
  TFile *inf, *outf;
  string title;
  DzeroUPCTreeMessenger *MDzeroUPC;
  TH1D *hDmass;
  TH1D *hNumberOfEventsAfterCuts;
  
  TTree *OutputTree;

  DataAnalyzer(const char *filename, const char *outFilename, const char *mytitle = "")
      : inf(new TFile(filename)), title(mytitle), MDzeroUPC(new DzeroUPCTreeMessenger(*inf, string("Tree"))),
        outf(new TFile(outFilename, "recreate")) {
    outf->cd();
    OutputTree = new TTree(Form("OutputTree%s", title.c_str()), "");
  }

  ~DataAnalyzer() {
    deleteHistograms();
    inf->Close();
    outf->Close();
    delete MDzeroUPC;
  }

  void analyze(GeneralInfoManager &man) {
    outf->cd();
    hDmass = new TH1D(Form("hDmass%s", title.c_str()), "", DMASSNBINS, DMASSMIN, DMASSMAX);
    hNumberOfEventsAfterCuts = new TH1D(Form("hNumberOfEventsAfterCuts%s", title.c_str()), "", 4, -0.5, 3.5);
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(1, "Total");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(2, "Event Selection");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(3, "Gen D Signal Event Selection");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(4, "Reco D Signal Event Selection");

    float HFEMaxPlus, HFEMaxMinus;
    float HFEMaxPlus_forest, HFEMaxMinus_forest;
    float HFEMaxPlus_eta5, HFEMaxMinus_eta5;
    float HFEMaxPlus_pt0p1, HFEMaxMinus_pt0p1;
    
    OutputTree->Branch("HFEMaxPlus", &HFEMaxPlus, "HFEMaxPlus/F");
    OutputTree->Branch("HFEMaxMinus", &HFEMaxMinus, "HFEMaxMinus/F");
    OutputTree->Branch("HFEMaxPlus_forest", &HFEMaxPlus_forest, "HFEMaxPlus_forest/F");
    OutputTree->Branch("HFEMaxMinus_forest", &HFEMaxMinus_forest, "HFEMaxMinus_forest/F");
    OutputTree->Branch("HFEMaxPlus_eta5", &HFEMaxPlus_eta5, "HFEMaxPlus_eta5/F");
    OutputTree->Branch("HFEMaxMinus_eta5", &HFEMaxMinus_eta5, "HFEMaxMinus_eta5/F");
    OutputTree->Branch("HFEMaxPlus_pt0p1", &HFEMaxPlus_pt0p1, "HFEMaxPlus_pt0p1/F");
    OutputTree->Branch("HFEMaxMinus_pt0p1", &HFEMaxMinus_pt0p1, "HFEMaxMinus_pt0p1/F");

    hDmass->Sumw2();
    hNumberOfEventsAfterCuts->Sumw2();

    double DzeroPTMin = man.GetCutParameterValue("DzeroPTMin");
    double DzeroPTMax = man.GetCutParameterValue("DzeroPTMax");
    double DzeroYMin = man.GetCutParameterValue("DzeroYMin");
    double DzeroYMax = man.GetCutParameterValue("DzeroYMax");
    int IsResolved = man.GetCutParameterValue("IsResolved");
    if (IsResolved != -1 && IsResolved != 0 && IsResolved != 1) {
      std::cout << "Invalid value for IsResolved parameter. It should be -1 (no selection), 0 (only direct), or 1 (only resolved)." << std::endl;
      return;
    }

    int useClusterCompatibilityFilter = man.GetCutParameterValue("useClusterCompatibilityFilter");
    if (useClusterCompatibilityFilter != 0 && useClusterCompatibilityFilter != 1)    {
      std::cout << "Invalid value for useClusterCompatibilityFilter parameter. It should be 0 (no filter) or 1 (apply cluster compatibility filter)." << std::endl;
      return;
    }

    int useGenDSignalEventSelection = man.GetCutParameterValue("useGenDSignalEventSelection");
    if (useGenDSignalEventSelection != 0 && useGenDSignalEventSelection != 1)    {
      std::cout << "Invalid value for useGenDSignalEventSelection parameter. It should be 0 (no D signal event selection) or 1 (only consider events with at least one gen-level D candidate that falls into the (pt,y) bin of interest as signal MC events)." << std::endl;
      return;
    }

    int useRecoDSignalEventSelection = man.GetCutParameterValue("useRecoDSignalEventSelection");
    if (useRecoDSignalEventSelection != 0 && useRecoDSignalEventSelection != 1)    {
      std::cout << "Invalid value for useRecoDSignalEventSelection parameter. It should be 0 (no reco D signal event selection) or 1 (only consider events with at least one reco-level D candidate that falls into the (pt,y) bin of interest as signal MC events)." << std::endl;
      return;
    }

    int DoSystD = man.GetCutParameterValue("DoSystD");
    if (DoSystD < 0 || DoSystD > 4)    {
      std::cout << "Invalid value for DoSystD parameter. It should be 0 (nominal), 1 (Dsvpv variation), 2 (DtrkPt variation), 3 (Dalpha variation), or 4 (Dchi2cl variation)." << std::endl;
      return;
    }

    unsigned long nEntry = MDzeroUPC->GetEntries();
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      if (i % 1000 == 0) {
        Bar.Update(i);
        Bar.Print();
      }
      hNumberOfEventsAfterCuts->Fill(0); // Total events

      // Event selection
      if (!eventSelection(MDzeroUPC, IsResolved, useClusterCompatibilityFilter))
        continue;
      hNumberOfEventsAfterCuts->Fill(1); // After event selection

      // Check if the event is a signal MC event, a.k.a., having at least one gen-level D candidate that falls into the (pt,y) bin of interest
      if (useGenDSignalEventSelection == 1) {
        bool isSigMCEvt = false;
        for (unsigned long j = 0; j < MDzeroUPC->Gpt->size(); j++) {
          if (MDzeroUPC->GisSignalCalc->at(j) == false)
            continue;

          if (MDzeroUPC->Gpt->at(j) < DzeroPTMin)
            continue;
          if (MDzeroUPC->Gpt->at(j) > DzeroPTMax)
            continue;
          if (MDzeroUPC->Gy->at(j) < DzeroYMin)
            continue;
          if (MDzeroUPC->Gy->at(j) > DzeroYMax)
            continue;
          isSigMCEvt = true;
        }
        if (!isSigMCEvt)
          continue;
      }
      hNumberOfEventsAfterCuts->Fill(2); // After gen-level D signal event selection

      // Do reco-level D signal event selection if needed
      if (useRecoDSignalEventSelection == 1)      {
        bool hasRecoDSignalCandidate = false;

        bool doTrkFilter = false;
        if (MDzeroUPC->Dtrk1PtErr != nullptr &&
            MDzeroUPC->Dtrk2PtErr != nullptr &&
            MDzeroUPC->Dtrk1PixelHit != nullptr &&
            MDzeroUPC->Dtrk1StripHit != nullptr &&
            MDzeroUPC->Dtrk2PixelHit != nullptr &&
            MDzeroUPC->Dtrk2StripHit != nullptr) doTrkFilter = true;

        for (unsigned long j = 0; j < MDzeroUPC->Dalpha->size(); j++) {
          if (MDzeroUPC->Dpt->at(j) < DzeroPTMin)
            continue;
          if (MDzeroUPC->Dpt->at(j) > DzeroPTMax)
            continue;
          if (MDzeroUPC->Dy->at(j) < DzeroYMin)
            continue;
          if (MDzeroUPC->Dy->at(j) > DzeroYMax)
            continue;
          if (MDzeroUPC->Dgen->at(j) != 23333) // Gen matching condition
            continue;
          if (DoSystD==0 && MDzeroUPC->DpassCut23PAS->at(j) == false) continue;
          // if (DoSystD==0 && MDzeroUPC->DpassCut23LowPt->at(j) == false) continue;
          if (DoSystD==1 && MDzeroUPC->DpassCut23PASSystDsvpvSig->at(j) == false) continue;
          if (DoSystD==2 && MDzeroUPC->DpassCut23PASSystDtrkPt->at(j) == false) continue;
          if (DoSystD==3 && MDzeroUPC->DpassCut23PASSystDalpha->at(j) == false) continue;
          if (DoSystD==4 && MDzeroUPC->DpassCut23PASSystDchi2cl->at(j) == false) continue;
//          if (DoSystD==0 && MDzeroUPC->DpassCutDefault->at(j) == false) continue;
//          if (DoSystD==1 && MDzeroUPC->DpassCutSystDsvpvSig->at(j) == false) continue;
//          if (DoSystD==2 && MDzeroUPC->DpassCutSystDtrkPt->at(j) == false) continue;
//          if (DoSystD==3 && MDzeroUPC->DpassCutSystDalpha->at(j) == false) continue;
//          if (DoSystD==4 && MDzeroUPC->DpassCutSystDchi2cl->at(j) == false) continue;
          if (doTrkFilter) {
            if (
              (MDzeroUPC->Dtrk1PtErr->at(j) / MDzeroUPC->Dtrk1Pt->at(j)) > 0.1 ||
              (MDzeroUPC->Dtrk2PtErr->at(j) / MDzeroUPC->Dtrk2Pt->at(j)) > 0.1
            ) continue;
           if (
             (MDzeroUPC->Dtrk1PixelHit->at(j) + MDzeroUPC->Dtrk1StripHit->at(j)) < 11 ||
             (MDzeroUPC->Dtrk2PixelHit->at(j) + MDzeroUPC->Dtrk2StripHit->at(j)) < 11
           ) continue;
          }

          hDmass->Fill((*MDzeroUPC->Dmass)[j]);
          hasRecoDSignalCandidate = true;
        } // end of reco-level Dzero loop

        if (!hasRecoDSignalCandidate)
          continue;
      }
      hNumberOfEventsAfterCuts->Fill(3); // After reco-level D signal event selection
      
      // Fill HF E_max distributions for all events
      HFEMaxPlus=MDzeroUPC->HFEMaxPlus;
      HFEMaxMinus=MDzeroUPC->HFEMaxMinus;
      HFEMaxPlus_forest=MDzeroUPC->HFEMaxPlus_forest;
      HFEMaxMinus_forest=MDzeroUPC->HFEMaxMinus_forest;
      HFEMaxPlus_eta5=MDzeroUPC->HFEMaxPlus_eta5;
      HFEMaxMinus_eta5=MDzeroUPC->HFEMaxMinus_eta5;
      HFEMaxPlus_pt0p1=MDzeroUPC->HFEMaxPlus_pt0p1;
      HFEMaxMinus_pt0p1=MDzeroUPC->HFEMaxMinus_pt0p1;
      OutputTree->Fill(); 
    }     // end of event loop
  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
    smartWrite(hDmass);
    smartWrite(hNumberOfEventsAfterCuts);
    OutputTree->Write();
  }

private:
  void deleteHistograms() {
    delete hDmass;
    delete hNumberOfEventsAfterCuts;
  }
};

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  string inputFile = CL.Get("Input", "mergedSample.root"); // Input file
  string outputFile = CL.Get("Output", "output.root");     // Output file

  double DzeroPTMin = CL.GetDouble("DzeroPTMin", 0.); // Minimum pt of gen-level D0 candidates to be considered as signal MC events
  double DzeroPTMax = CL.GetDouble("DzeroPTMax", 100.); // Maximum pt of gen-level D0 candidates to be considered as signal MC events
  double DzeroYMin = CL.GetDouble("DzeroYMin", -10.); // Minimum rapidity of gen-level D0 candidates to be considered as signal MC events
  double DzeroYMax = CL.GetDouble("DzeroYMax", 10.); // Maximum rapidity of gen-level D0 candidates to be considered as signal MC events
  int IsResolved = CL.GetInt("isResolved", -1); // -1 = no selection, 0 = only direct, 1 = only resolved 
  int useClusterCompatibilityFilter = CL.GetInt("useClusterCompatibilityFilter", 0); // 0 = no filter, 1 = apply cluster compatibility filter
  int useGenDSignalEventSelection = CL.GetInt("useGenDSignalEventSelection", 1); // 0 = no D signal event selection, 1 = only consider events with at least one gen-level D candidate that falls into the (pt,y) bin of interest as signal MC events
  int useRecoDSignalEventSelection = CL.GetInt("useRecoDSignalEventSelection", 0); // 0 = no reco D signal event selection, 1 = only consider events with at least one reco-level, gen matched D candidate that falls into the (pt,y) bin of interest as signal MC events
  int DoSystD = CL.GetInt("DoSystD", 0);             // Systematic study: apply the alternative D selections
                                                     // 0 = nominal, 1 = Dsvpv variation, 2: DtrkPt variation
                                                     // 3 = Dalpha variation, 4: Dchi2cl variation
  int IsGammaN = CL.GetInt("IsGammaN", 0); // 0 = N-Gamma, 1 = Gamma-N

  DataAnalyzer analyzer(inputFile.c_str(), outputFile.c_str());

  GeneralInfoManager man(analyzer.outf, "InfoDir", false);
  TTimeStamp *currentTime = new TTimeStamp();
  
  man.AddSourceFile(inputFile, currentTime);
  man.AddCutParameter("DzeroPTMin", DzeroPTMin, currentTime);
  man.AddCutParameter("DzeroPTMax", DzeroPTMax, currentTime);
  man.AddCutParameter("DzeroYMin", DzeroYMin, currentTime);
  man.AddCutParameter("DzeroYMax", DzeroYMax, currentTime);
  man.AddCutParameter("IsResolved", IsResolved, currentTime);
  man.AddCutParameter("useClusterCompatibilityFilter", useClusterCompatibilityFilter, currentTime);
  man.AddCutParameter("useGenDSignalEventSelection", useGenDSignalEventSelection, currentTime);
  man.AddCutParameter("useRecoDSignalEventSelection", useRecoDSignalEventSelection, currentTime);
  man.AddCutParameter("DoSystD", DoSystD, currentTime);
  man.AddCutParameter("IsGammaN", IsGammaN, currentTime);

  std::cout << "Parameters used for this analysis:" << std::endl;
  man.PrintInfo();

  analyzer.analyze(man);
  analyzer.writeHistograms(analyzer.outf);
  man.SaveToFile();
  std::cout << "Analysis complete. Output saved to " << outputFile << std::endl;

  return 0;
}