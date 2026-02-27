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

bool eventSelection(DzeroUPCTreeMessenger *b) {  
  // [FIXME] Cluster comp? 
  if (b->selectedBkgFilter == false || b->selectedVtxFilter == false)
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
  TH1D *hHFEMaxPlus, *hHFEMaxMinus;
  TH2D *hHFEMaxMinus_vs_hHFEMaxPlus;

  TH1D *hHFEMaxPlus_forest, *hHFEMaxMinus_forest;
  TH2D *hHFEMaxMinus_vs_hHFEMaxPlus_forest;

  TH1D *hHFEMaxPlus_eta5, *hHFEMaxMinus_eta5;
  TH2D *hHFEMaxMinus_vs_hHFEMaxPlus_eta5;

  TH1D *hHFEMaxPlus_pt0p1, *hHFEMaxMinus_pt0p1;
  TH2D *hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1;

  DataAnalyzer(const char *filename, const char *outFilename, const char *mytitle = "")
      : inf(new TFile(filename)), title(mytitle), MDzeroUPC(new DzeroUPCTreeMessenger(*inf, string("Tree"))),
        outf(new TFile(outFilename, "recreate")) {
    outf->cd();
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

    hHFEMaxPlus = new TH1D(Form("hHFEMaxPlus%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus = new TH1D(Form("hHFEMaxMinus%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_vs_hHFEMaxPlus = new TH2D(Form("hHFEMaxMinus_vs_hHFEMaxPlus%s", title.c_str()), "", 400, 0, 80, 400, 0, 80);

    hHFEMaxPlus_forest = new TH1D(Form("hHFEMaxPlus_forest%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_forest = new TH1D(Form("hHFEMaxMinus_forest%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_vs_hHFEMaxPlus_forest = new TH2D(Form("hHFEMaxMinus_vs_hHFEMaxPlus_forest%s", title.c_str()), "", 400, 0, 80, 400, 0, 80);

    hHFEMaxPlus_eta5 = new TH1D(Form("hHFEMaxPlus_eta5%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_eta5 = new TH1D(Form("hHFEMaxMinus_eta5%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_vs_hHFEMaxPlus_eta5 = new TH2D(Form("hHFEMaxMinus_vs_hHFEMaxPlus_eta5%s", title.c_str()), "", 400, 0, 80, 400, 0, 80);

    hHFEMaxPlus_pt0p1 = new TH1D(Form("hHFEMaxPlus_pt0p1%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_pt0p1 = new TH1D(Form("hHFEMaxMinus_pt0p1%s", title.c_str()), "", 400, 0, 80);
    hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1 = new TH2D(Form("hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1%s", title.c_str()), "", 400, 0, 80, 400, 0, 80);

    hDmass->Sumw2();
    hHFEMaxPlus->Sumw2();
    hHFEMaxMinus->Sumw2();
    hHFEMaxMinus_vs_hHFEMaxPlus->Sumw2();
    hHFEMaxPlus_forest->Sumw2();
    hHFEMaxMinus_forest->Sumw2();
    hHFEMaxMinus_vs_hHFEMaxPlus_forest->Sumw2();
    hHFEMaxPlus_eta5->Sumw2();
    hHFEMaxMinus_eta5->Sumw2();
    hHFEMaxMinus_vs_hHFEMaxPlus_eta5->Sumw2();
    hHFEMaxPlus_pt0p1->Sumw2();
    hHFEMaxMinus_pt0p1->Sumw2();
    hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1->Sumw2();

    unsigned long nEntry = MDzeroUPC->GetEntries();
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    int triggerChoice = man.GetCutParameterValue("TriggerChoice");
    if (triggerChoice != 0 && triggerChoice != -1 && triggerChoice != 1) {
      std::cerr << "Error: Invalid TriggerChoice parameter value. Expected 0, -1, or 1." << std::endl;
      return;
    }

    int nTrkFilter = man.GetCutParameterValue("nTrkFilter");
    if (nTrkFilter != 0 && nTrkFilter != 1 && nTrkFilter != -1) {
      std::cerr << "Error: Invalid nTrkFilter parameter value. Expected 0, 1, or -1." << std::endl;
      return;
    }

    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      if (i % 1000 == 0) {
        Bar.Update(i);
        Bar.Print();
      }

      // Trigger selection  
      if (triggerChoice == 0 && !isNotBptxOR) continue; // isNotBptxOR
      else if (triggerChoice == -1 && !isUnpairedBunchBptxMinus) continue; // isUnpairedBunchBptxMinus
      else if (triggerChoice == 1 && !isUnpairedBunchBptxPlus) continue; // isUnpairedBunchBptxPlus

      // Track filter selection
      if (nTrkFilter == 1 && MDzeroUPC->nTrackInAcceptanceHP != 0) continue; // nTrkInAcceptanceHP == 0
      else if (nTrkFilter == -1 && !(MDzeroUPC->nTrackInAcceptanceHP > 0)) continue; // nTrkInAcceptanceHP > 0
      
      // Fill HF E_max distributions for all events
      hHFEMaxMinus->Fill(MDzeroUPC->HFEMaxMinus);
      hHFEMaxPlus->Fill(MDzeroUPC->HFEMaxPlus);
      hHFEMaxMinus_vs_hHFEMaxPlus->Fill(MDzeroUPC->HFEMaxMinus, MDzeroUPC->HFEMaxPlus);

      hHFEMaxMinus_forest->Fill(MDzeroUPC->HFEMaxMinus_forest);
      hHFEMaxPlus_forest->Fill(MDzeroUPC->HFEMaxPlus_forest);
      hHFEMaxMinus_vs_hHFEMaxPlus_forest->Fill(MDzeroUPC->HFEMaxMinus_forest, MDzeroUPC->HFEMaxPlus_forest);

      hHFEMaxPlus_eta5->Fill(MDzeroUPC->HFEMaxPlus_eta5);
      hHFEMaxMinus_eta5->Fill(MDzeroUPC->HFEMaxMinus_eta5);
      hHFEMaxMinus_vs_hHFEMaxPlus_eta5->Fill(MDzeroUPC->HFEMaxMinus_eta5, MDzeroUPC->HFEMaxPlus_eta5);

      hHFEMaxPlus_pt0p1->Fill(MDzeroUPC->HFEMaxPlus_pt0p1);
      hHFEMaxMinus_pt0p1->Fill(MDzeroUPC->HFEMaxMinus_pt0p1);
      hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1->Fill(MDzeroUPC->HFEMaxMinus_pt0p1, MDzeroUPC->HFEMaxPlus_pt0p1);
    }     // end of event loop
  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
    smartWrite(hHFEMaxPlus);
    smartWrite(hHFEMaxMinus);
    smartWrite(hHFEMaxMinus_vs_hHFEMaxPlus);

    smartWrite(hHFEMaxPlus_forest);
    smartWrite(hHFEMaxMinus_forest);
    smartWrite(hHFEMaxMinus_vs_hHFEMaxPlus_forest);

    smartWrite(hHFEMaxPlus_eta5);
    smartWrite(hHFEMaxMinus_eta5);
    smartWrite(hHFEMaxMinus_vs_hHFEMaxPlus_eta5);

    smartWrite(hHFEMaxPlus_pt0p1);
    smartWrite(hHFEMaxMinus_pt0p1);
    smartWrite(hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1);
  }

private:
  void deleteHistograms() {
    delete hHFEMaxPlus;
    delete hHFEMaxMinus;
    delete hHFEMaxMinus_vs_hHFEMaxPlus;

    delete hHFEMaxPlus_forest;
    delete hHFEMaxMinus_forest;
    delete hHFEMaxMinus_vs_hHFEMaxPlus_forest;

    delete hHFEMaxPlus_eta5;
    delete hHFEMaxMinus_eta5;
    delete hHFEMaxMinus_vs_hHFEMaxPlus_eta5;

    delete hHFEMaxPlus_pt0p1;
    delete hHFEMaxMinus_pt0p1;
    delete hHFEMaxMinus_vs_hHFEMaxPlus_pt0p1;
  }
};

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  int TriggerChoice = CL.GetInt("TriggerChoice", 0); // 0 = isNotBptxOR, -1 = isUnpairedBunchBptxMinus, 1 = isUnpairedBunchBptxPlus
  int nTrkFilter = CL.GetInt("nTrkFilter", 0); // 0 = no filter, 1 = nTrkInAcceptanceHP == 0, -1 = nTrkInAcceptanceHP > 0

  string inputFile = CL.Get("Input", "mergedSample.root"); // Input file
  string outputFile = CL.Get("Output", "output.root");     // Output file

  DataAnalyzer analyzer(inputFile.c_str(), outputFile.c_str());

  GeneralInfoManager man(analyzer.outf, "InfoDir", false);
  TTimeStamp *currentTime = new TTimeStamp();
  
  man.AddSourceFile(inputFile, currentTime);
  man.AddCutParameter("TriggerChoice", TriggerChoice, currentTime); // 0 = no trigger selection, 1 = isL1ZDCOr, 2 = isL1ZDCXORJet8
  man.AddCutParameter("nTrkFilter", nTrkFilter, currentTime);

  std::cout << "Parameters used for this analysis:" << std::endl;
  man.PrintInfo();

  analyzer.analyze(man);
  analyzer.writeHistograms(analyzer.outf);
  man.SaveToFile();
  std::cout << "Analysis complete. Output saved to " << outputFile << std::endl;

  return 0;
}