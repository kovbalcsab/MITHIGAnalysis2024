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

using namespace std;
#include "CommandLine.h" // Yi's Commandline bundle
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "helpMessage.h" // Print out help message
#include "parameter.h"   // The parameters used in the analysis
#include "utilities.h"   // Yen-Jie's random utility functions

//============================================================//
// Function to check for configuration errors
//============================================================//
bool checkError(const Parameters &par) { return false; }

//======= trackSelection =====================================//
// Check if the track pass selection criteria
//============================================================//
bool dzeroSelection(DzeroUPCTreeMessenger *b, Parameters par, int j) { return true; }

//======= eventSelection =====================================//
// Check if the event pass eventSelection criteria
//============================================================//
bool eventSelection(DzeroUPCTreeMessenger *b, const Parameters &par) {
  if (par.IsData == false) return true;
  if (b->selectedBkgFilter == false || b->selectedVtxFilter == false) return false;
  if (par.IsGammaN && (b->ZDCgammaN && b->gapgammaN) == false) return false;
  if (!par.IsGammaN && (b->ZDCNgamma && b->gapNgamma) == false) return false;
  if (par.TriggerChoice == 1 && b->isL1ZDCOr == false) return false;
  if (par.TriggerChoice == 2 && b->isL1ZDCXORJet8 == false) return false;
  if (b->nVtx >= 3) return false;
  return true;
}

class DataAnalyzer {
public:
  TFile *inf, *outf;
  TH1D *hDmass;
  DzeroUPCTreeMessenger *MDzeroUPC;
  TNtuple *nt;
  string title;
  TH1D *hDpt;
  TH1D *hDptGen;


  DataAnalyzer(const char *filename, const char *outFilename, const char *mytitle = "")
      : inf(new TFile(filename)), MDzeroUPC(new DzeroUPCTreeMessenger(*inf, string("Tree"))), title(mytitle),
        outf(new TFile(outFilename, "recreate")) {
    outf->cd();
    nt = new TNtuple("nt", "D0 mass tree", "Dmass:Dgen");
  }

  ~DataAnalyzer() {
    deleteHistograms();
    delete nt;
    inf->Close();
    outf->Close();
    delete MDzeroUPC;
  }

  void analyze(Parameters &par) {
    outf->cd();
    hDmass = new TH1D(Form("hDmass%s", title.c_str()), "", 60, 1.7, 2.0);
    hDpt = new TH1D(Form("hDpt%s", title.c_str()), "", 100, 0, 20);
    hDptGen = new TH1D(Form("hDptGen%s", title.c_str()), "", 100, 0, 20);

    hDmass->Sumw2();
    hDpt->Sumw2();
    hDptGen->Sumw2();

    par.printParameters();
    unsigned long nEntry = MDzeroUPC->GetEntries() * par.scaleFactor;
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    //===============================================================
    // data loop
    //===============================================================
    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      if (i % 1000 == 0) {
        Bar.Update(i);
        Bar.Print();
      }
      // Check if the event passes the selection criteria
      if (eventSelection(MDzeroUPC, par)) {
        for (unsigned long j = 0; j < MDzeroUPC->Dalpha->size(); j++) {
          // pt cuts
          if (MDzeroUPC->Dpt->at(j) < par.MinDzeroPT) continue;
          if (MDzeroUPC->Dpt->at(j) > par.MaxDzeroPT) continue;
          // rapidity cuts
          if (MDzeroUPC->Dy->at(j) < par.MinDzeroY) continue;
          if (MDzeroUPC->Dy->at(j) > par.MaxDzeroY) continue;
          if (MDzeroUPC->DpassCut->at(j) == false) continue;
          if(!par.DoInclusive && par.DoPrompt && !(*MDzeroUPC->DisSignalCalcPrompt)[j] ) continue;
          if(!par.DoInclusive && par.DoNonPrompt && !(*MDzeroUPC->DisSignalCalcFeeddown)[j] )continue;
          hDmass->Fill((*MDzeroUPC->Dmass)[j]);
          if (!par.IsData) {
            nt->Fill((*MDzeroUPC->Dmass)[j], (*MDzeroUPC->Dgen)[j]);
            if (MDzeroUPC->Dgen->at(j) == 23333) {
              hDpt->Fill((*MDzeroUPC->Dpt)[j]);
            }
          } else
            nt->Fill((*MDzeroUPC->Dmass)[j], 0);
        } // end of reco-level Dzero loop
      }   // end of event selection
      //===============================================================

      //===============================================================
      // mc loop
      //===============================================================
      if (!par.IsData) {
        for (unsigned long j = 0; j < MDzeroUPC->Gpt->size(); j++) {
          if (MDzeroUPC->Gpt->at(j) < par.MinDzeroPT) continue;
          if (MDzeroUPC->Gpt->at(j) > par.MaxDzeroPT) continue;
          if (MDzeroUPC->Gy->at(j) < par.MinDzeroY) continue;
          if (MDzeroUPC->Gy->at(j) > par.MaxDzeroY)continue;
          if (MDzeroUPC->GisSignalCalc->at(j) == false) continue;
          if(!par.DoInclusive && par.DoPrompt && !(MDzeroUPC->GisSignalCalcPrompt->at(j)) ) continue;
          if(!par.DoInclusive && par.DoNonPrompt && !(MDzeroUPC->GisSignalCalcFeeddown->at(j)) )continue;
          hDptGen->Fill(MDzeroUPC->Gpt->at(j));
        } // end of gen-level Dzero loop
      }   // end of gen-level Dzero loop
    }     // end of event loop
    //===============================================================

  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
    smartWrite(hDmass);
    hDpt->Write();
    hDptGen->Write();
    smartWrite(nt);
  }

private:
  void deleteHistograms() {
    delete hDmass;
    delete hDpt;
    delete hDptGen;
  }
};

//============================================================//
// Main analysis
//============================================================//
int main(int argc, char *argv[]) {
  if (printHelpMessage(argc, argv))
    return 0;
  CommandLine CL(argc, argv);
  float MinDzeroPT = CL.GetDouble("MinDzeroPT", 1);  // Minimum Dzero transverse momentum threshold for Dzero selection.
  float MaxDzeroPT = CL.GetDouble("MaxDzeroPT", 2);  // Maximum Dzero transverse momentum threshold for Dzero selection.
  float MinDzeroY = CL.GetDouble("MinDzeroY", -2);   // Minimum Dzero rapidity threshold for Dzero selection.
  float MaxDzeroY = CL.GetDouble("MaxDzeroY", +2);   // Maximum Dzero rapidity threshold for Dzero selection.
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  int TriggerChoice = CL.GetInt("TriggerChoice", 2); // 0 = no trigger sel, 1 = isL1ZDCOr, 2 = isL1ZDCXORJet8
  float scaleFactor = CL.GetDouble("scaleFactor", 1); // Scale factor for the number of events to be processed.
  bool IsData = CL.GetBool("IsData", 0);              // Data or MC
  bool DoInclusive = CL.GetBool("DoInclusive", 1);    // Whether or not we want there to be an inclusive selection, used mostly as a protection
  bool DoPrompt    = CL.GetBool("DoPrompt", 0);       // Whether to apply the prompt selection
  bool DoNonPrompt = CL.GetBool("DoNonPrompt", 0);    // Whether to apply the non-prompt selection
  Parameters par(MinDzeroPT, MaxDzeroPT, MinDzeroY, MaxDzeroY, IsGammaN, TriggerChoice, IsData, scaleFactor, DoInclusive, DoPrompt, DoNonPrompt);
  par.input = CL.Get("Input", "mergedSample.root"); // Input file
  par.output = CL.Get("Output", "output.root");     // Output file
  par.nThread = CL.GetInt("nThread", 1);            // The number of threads to be used for parallel processing.
  par.nChunk =
      CL.GetInt("nChunk", 1); // Specifies which chunk (segment) of the data to process, used in parallel processing.
  if (checkError(par))
    return -1;
  std::cout << "Parameters are set" << std::endl;
  // Analyze Data
  DataAnalyzer analyzer(par.input.c_str(), par.output.c_str(), "");
  analyzer.analyze(par);
  analyzer.writeHistograms(analyzer.outf);
  saveParametersToHistograms(par, analyzer.outf);
  cout << "done!" << analyzer.outf->GetName() << endl;
}
