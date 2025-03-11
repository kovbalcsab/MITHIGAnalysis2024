#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>
#include <vector>

using namespace std;
#include "CommandLine.h" // Yi's Commandline bundle
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "helpMessage.h" // Print out help message
#include "parameter.h"   // The parameters used in the analysis
#include "utilities.h"   // Yen-Jie's random utility functions

// How to run: ./ExtractGenlevelHists --Input "your_input_mc_file.root" \
										--IsGammaN ($IsGammaN: 0 or 1) \
										--Output "genLevelHists.root" &

class DataAnalyzer {
public:
  TFile *inf, *outf;
  DzeroUPCTreeMessenger *MDzeroUPC;
  string title;

  TH2D *hGDpt_GDy; 

  DataAnalyzer(const char *filename, const char *outFilename, const char *mytitle = "")
      : inf(new TFile(filename)), MDzeroUPC(new DzeroUPCTreeMessenger(*inf, string("Tree"))), title(mytitle),
        outf(new TFile(outFilename, "recreate")) {
    outf->cd();
  }

  ~DataAnalyzer() {
    deleteHistograms();
    inf->Close();
    outf->Close();
    delete MDzeroUPC;
  }

  void analyze(Parameters &par) {
    outf->cd();

    initiateHistograms(par);

    par.printParameters();
    unsigned long nEntry = MDzeroUPC->GetEntries() * par.scaleFactor;
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      if (i % 1000 == 0) {
        Bar.Update(i);
        Bar.Print();
      }

      for (unsigned long j = 0; j < MDzeroUPC->Gpt->size(); j++) {
        if (MDzeroUPC->Gpt->at(j) < par.MinDzeroPT)
          continue;
        if (MDzeroUPC->Gpt->at(j) > par.MaxDzeroPT)
          continue;
        if (MDzeroUPC->Gy->at(j) < par.MinDzeroY)
          continue;
        if (MDzeroUPC->Gy->at(j) > par.MaxDzeroY)
          continue;
        if (MDzeroUPC->GisSignalCalc->at(j) == false)
          continue;
        hGDpt_GDy->Fill(MDzeroUPC->Gpt->at(j), MDzeroUPC->Gy->at(j));
      }
    }     // end of event loop
  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();

    smartWrite(hGDpt_GDy);
  }

private:
  void deleteHistograms() {
    delete hGDpt_GDy;
  }

  void initiateHistograms(Parameters &par) {
    hGDpt_GDy = new TH2D(Form("hGDpt_GDy%s", title.c_str()), "", 30, 0, 6, 30, -3, 3);

    hGDpt_GDy->Sumw2();
  }
};

//============================================================//
// Main analysis
//============================================================//
int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  float scaleFactor = CL.GetDouble("scaleFactor", 1); // Scale factor for the number of events to be processed.

  Parameters par(0, 6, -3, 3, IsGammaN, 0, 0, scaleFactor,
                 0, 0);
  par.input = CL.Get("Input", "mergedSample.root"); // Input file
  par.output = CL.Get("Output", "output.root");     // Output file
  par.nThread = CL.GetInt("nThread", 1);            // The number of threads to be used for parallel processing.
  par.nChunk =
      CL.GetInt("nChunk", 1); // Specifies which chunk (segment) of the data to process, used in parallel processing.
  // Analyze Data
  DataAnalyzer analyzer(par.input.c_str(), par.output.c_str(), "");
  analyzer.analyze(par);
  analyzer.writeHistograms(analyzer.outf);
  saveParametersToHistograms(par, analyzer.outf);
  cout << "done!" << analyzer.outf->GetName() << endl;
}
