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
#include "InfoManager.h"
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "utilities.h"   // Yen-Jie's random utility functions

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  string inputFile = CL.Get("Input", "mergedSample.root"); // Input file
  string outputFile = CL.Get("Output", "output.txt");      // Output file
  double quantile = CL.GetDouble("Quantile", 0.99);        // Quantile to compute

  TFile *inf = new TFile(inputFile.c_str());
  if (!inf || inf->IsZombie()) {
    cerr << "Error: Could not open input file " << inputFile << endl;
    return 1;
  }
  TTree *tree = (TTree *)inf->Get("OutputTree");

  TH1D *HFEMaxPlusHist =
      new TH1D("HFEMaxPlusHist", "HFEMaxPlus Distribution", 2000, 0, 200); // Adjust bins and range as needed

  TH1D *HFEMaxMinusHist =
      new TH1D("HFEMaxMinusHist", "HFEMaxMinus Distribution", 2000, 0, 200); // Adjust bins and range as needed
                                                                             //
  // fill histograms
  tree->Draw("HFEMaxPlus_forest>>HFEMaxPlusHist");
  tree->Draw("HFEMaxMinus_forest>>HFEMaxMinusHist");

  double xp_minus, xp_plus;
  double p = quantile;

  HFEMaxMinusHist->GetQuantiles(1, &xp_minus, &p);
  HFEMaxPlusHist->GetQuantiles(1, &xp_plus, &p);

  std::ofstream outFile(outputFile);
  outFile << "Quantile: " << quantile << endl;
  outFile << "HFEMaxPlus_forest pos: " << xp_plus << endl;
  outFile << "HFEMaxMinus_forest pos: " << xp_minus << endl;
  outFile.close();

  delete inf;

  return 0;
}
