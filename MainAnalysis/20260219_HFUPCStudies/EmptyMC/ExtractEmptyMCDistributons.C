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

class DataAnalyzer {
public:
  TFile *inf, *outf;
  string title;
  DzeroUPCTreeMessenger *MDzeroUPC;
  
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

    hNumberOfEventsAfterCuts = new TH1D(Form("hNumberOfEventsAfterCuts%s", title.c_str()), "", 1, -0.5, 0.5);
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(1, "NoCuts");
    hNumberOfEventsAfterCuts->Sumw2();

    unsigned long nEntry = MDzeroUPC->GetEntries();
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      if (i % 1000 == 0) {
        Bar.Update(i);
        Bar.Print();
      }
      hNumberOfEventsAfterCuts->Fill(0);

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
    OutputTree->Write();
    smartWrite(hNumberOfEventsAfterCuts);
  }

private:
  void deleteHistograms() {
    delete hNumberOfEventsAfterCuts;
  }
};

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  string inputFile = CL.Get("Input", "mergedSample.root"); // Input file
  string outputFile = CL.Get("Output", "output.root");     // Output file

  DataAnalyzer analyzer(inputFile.c_str(), outputFile.c_str());

  GeneralInfoManager man(analyzer.outf, "InfoDir", false);
  TTimeStamp *currentTime = new TTimeStamp();
  
  man.AddSourceFile(inputFile, currentTime);

  std::cout << "Parameters used for this analysis:" << std::endl;
  man.PrintInfo();

  analyzer.analyze(man);
  analyzer.writeHistograms(analyzer.outf);
  man.SaveToFile();
  std::cout << "Analysis complete. Output saved to " << outputFile << std::endl;

  return 0;
}