#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>
using namespace std;

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"
#include "TTree.h"

#include "CommandLine.h"

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  string InputFileName = CL.Get("Input");
  string OutputFileName = CL.Get("Output");
  int ADCThreshold = CL.GetInt("ADCThreshold", 19);

  TFile *InputFile = TFile::Open(InputFileName.c_str());
  if (!InputFile || InputFile->IsZombie()) {
    cerr << "Error: Could not open input file " << InputFileName << endl;
    return -1;
  }

  TTree *HiEvtTree = (TTree *)InputFile->Get("hiEvtAnalyzer/HiTree");
  TTree *HFADCTree = (TTree *)InputFile->Get("HFAdcana/adc");

  int HFPlusADC, HFMinusADC;
  float HFPlusLeading, HFMinusLeading;
  HiEvtTree->SetBranchAddress("hiHFPlus_pfle1", &HFPlusLeading);
  HiEvtTree->SetBranchAddress("hiHFMinus_pfle1", &HFMinusLeading);
  HFADCTree->SetBranchAddress("mMaxL1HFAdcPlus", &HFPlusADC);
  HFADCTree->SetBranchAddress("mMaxL1HFAdcMinus", &HFMinusADC);

  TH1D *hHFPlusLeading_Xn0nAll = new TH1D("hHFPlusLeading_Xn0nAll",
                                          "HF+ Leading Energy (ZB + at least one side ADC <=19);Energy (GeV);Events",
                                          400, 0, 50);
  TH1D *hHFMinusLeading_Xn0nAll = new TH1D("hHFMinusLeading_Xn0nAll",
                                           "HF- Leading Energy (ZB + at least one side ADC <=19);Energy (GeV);Events",
                                           400, 0, 50);
  TH1D *hHFPlusLeading_Xn0nOneSide =
      new TH1D("hHFPlusLeading_Xn0nOneSide",
               "HF+ Leading Energy (ZB + side with ADC<=19 filled);Energy (GeV);Events",
               400, 0, 50);
  TH1D *hHFMinusLeading_Xn0nOneSide =
      new TH1D("hHFMinusLeading_Xn0nOneSide",
               "HF- Leading Energy (ZB + side with ADC<=19 filled);Energy (GeV);Events",
               400, 0, 50);
  TH1D *hHFPlusLeading_0n0n =
      new TH1D("hHFPlusLeading_0n0n", "HF+ Leading Energy (ZB+ both ADC<=19);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0n0n =
      new TH1D("hHFMinusLeading_0n0n", "HF- Leading Energy (ZB+ both ADC<=19);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_NoCut =
      new TH1D("hHFPlusLeading_NoCut", "HF+ Leading Energy (ZB events);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_NoCut =
      new TH1D("hHFMinusLeading_NoCut", "HF- Leading Energy (ZB events);Energy (GeV);Events", 400, 0, 50);

  int nEntries = HiEvtTree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    HiEvtTree->GetEntry(i);
    HFADCTree->GetEntry(i);

    if (HFPlusADC <= ADCThreshold || HFMinusADC <= ADCThreshold) {
      hHFMinusLeading_Xn0nAll->Fill(HFMinusLeading);
      hHFPlusLeading_Xn0nAll->Fill(HFPlusLeading);
    }
    if (HFPlusADC <= ADCThreshold && HFMinusADC <= ADCThreshold) {
      hHFMinusLeading_0n0n->Fill(HFMinusLeading);
      hHFPlusLeading_0n0n->Fill(HFPlusLeading);
    }
    if (HFPlusADC <= ADCThreshold) {
      hHFPlusLeading_Xn0nOneSide->Fill(HFPlusLeading);
    }
    if (HFMinusADC <= ADCThreshold) {
      hHFMinusLeading_Xn0nOneSide->Fill(HFMinusLeading);
    }
    hHFPlusLeading_NoCut->Fill(HFPlusLeading);
    hHFMinusLeading_NoCut->Fill(HFMinusLeading);
  }

  TFile *OutputFile = new TFile(OutputFileName.c_str(), "RECREATE");
  OutputFile->cd();
  hHFPlusLeading_Xn0nAll->Write();
  hHFMinusLeading_Xn0nAll->Write();
  hHFPlusLeading_Xn0nOneSide->Write();
  hHFMinusLeading_Xn0nOneSide->Write();
  hHFPlusLeading_0n0n->Write();
  hHFMinusLeading_0n0n->Write();
  hHFPlusLeading_NoCut->Write();
  hHFMinusLeading_NoCut->Write();

  OutputFile->Close();

  return 0;
}
