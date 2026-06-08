#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>
using namespace std;

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TTimeStamp.h"
#include "TTree.h"

#include "CommandLine.h"

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  string InputFileName = CL.Get("Input");
  string OutputFileName = CL.Get("Output");
  int ADCThreshold = CL.GetInt("ADCThreshold", 19);
  int TriggerChoice =
      CL.GetInt("TriggerChoice", 0); // 0: no trigger cut, 1: HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v1
  double ZDCThreshold = CL.GetDouble("ZDCThreshold", 600.0);

  TFile *InputFile = TFile::Open(InputFileName.c_str());
  if (!InputFile || InputFile->IsZombie()) {
    cerr << "Error: Could not open input file " << InputFileName << endl;
    return -1;
  }

  TTree *HiEvtTree = (TTree *)InputFile->Get("hiEvtAnalyzer/HiTree");
  TTree *HFADCTree = (TTree *)InputFile->Get("HFAdcana/adc");
  TTree *ZDCTree = (TTree *)InputFile->Get("zdcanalyzer/zdcrechit");
  TTree *TriggerTree = (TTree *)InputFile->Get("hltanalysis/HltTree");
  TTree *SkimTree = (TTree *)InputFile->Get("skimanalysis/HltTree");
  TTree *TrackTree = (TTree *)InputFile->Get("PbPbTracks/trackTree");

  int HFPlusADC, HFMinusADC;
  int HFPlusADC_phi, HFPlusADC_eta, HFMinusADC_phi, HFMinusADC_eta;
  float HFPlusLeading, HFMinusLeading;
  float ZDCSumPlus, ZDCSumMinus;
  int HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v1;
  int pprimaryVertexFilter;
  std::vector<float> *TrackPt = nullptr;
  std::vector<float> *TrackEta = nullptr;
  std::vector<bool> *TrackHP = nullptr;

  HiEvtTree->SetBranchAddress("hiHFPlus_pfle1", &HFPlusLeading);
  HiEvtTree->SetBranchAddress("hiHFMinus_pfle1", &HFMinusLeading);
  HFADCTree->SetBranchAddress("mMaxL1HFAdcPlus", &HFPlusADC);
  HFADCTree->SetBranchAddress("mMaxL1HFAdcMinus", &HFMinusADC);
  HFADCTree->SetBranchAddress("mMaxietaPlus", &HFPlusADC_eta);
  HFADCTree->SetBranchAddress("mMaxiphiPlus", &HFPlusADC_phi);
  HFADCTree->SetBranchAddress("mMaxietaMinus", &HFMinusADC_eta);
  HFADCTree->SetBranchAddress("mMaxiphiMinus", &HFMinusADC_phi);

  ZDCTree->SetBranchAddress("sumPlus", &ZDCSumPlus);
  ZDCTree->SetBranchAddress("sumMinus", &ZDCSumMinus);
  TriggerTree->SetBranchAddress("HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v1",
                                &HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v1);
  SkimTree->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  TrackTree->SetBranchAddress("trkPt", &TrackPt);
  TrackTree->SetBranchAddress("trkEta", &TrackEta);
  TrackTree->SetBranchAddress("highPurity", &TrackHP);

  TH1D *ZDCSumPlusHist = new TH1D("ZDCSumPlusHist", "ZDC Sum Plus;Energy (GeV);Events", 400, 0, 4000);
  TH1D *ZDCSumPlus_0nHist = new TH1D("ZDCSumPlus_0nHist", "ZDC Sum Plus;Energy (GeV);Events", 400, 0, 4000);
  TH1D *ZDCSumMinusHist = new TH1D("ZDCSumMinusHist", "ZDC Sum Minus;Energy (GeV);Events", 400, 0, 4000);
  TH1D *ZDCSumMinus_0nHist = new TH1D("ZDCSumMinus_0nHist", "ZDC Sum Minus;Energy (GeV);Events", 400, 0, 4000);

  TH2D *HFPlusEtaPhi = new TH2D("HFPlusEtaPhi", "HF+ Leading ADC Eta-Phi;Eta;Phi", 100, 0, 100, 100, 0, 100);
  TH2D *HFMinusEtaPhi = new TH2D("HFMinusEtaPhi", "HF- Leading ADC Eta-Phi;Eta;Phi", 100, -100, 0, 100, 0, 100);

  TH1D *hHFPlusLeading_Xn0nAll = new TH1D(
      "hHFPlusLeading_Xn0nAll", "HF+ Leading Energy (ZB + at least one side ADC <=19);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_Xn0nAll =
      new TH1D("hHFMinusLeading_Xn0nAll", "HF- Leading Energy (ZB + at least one side ADC <=19);Energy (GeV);Events",
               400, 0, 50);
  TH1D *hHFPlusLeading_Xn0nOneSide =
      new TH1D("hHFPlusLeading_Xn0nOneSide", "HF+ Leading Energy (ZB + side with ADC<=19 filled);Energy (GeV);Events",
               400, 0, 50);
  TH1D *hHFMinusLeading_Xn0nOneSide =
      new TH1D("hHFMinusLeading_Xn0nOneSide", "HF- Leading Energy (ZB + side with ADC<=19 filled);Energy (GeV);Events",
               400, 0, 50);
  TH1D *hHFPlusLeading_0n0n =
      new TH1D("hHFPlusLeading_0n0n", "HF+ Leading Energy (ZB+ both ADC<=19);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0n0n =
      new TH1D("hHFMinusLeading_0n0n", "HF- Leading Energy (ZB+ both ADC<=19);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_NoCut =
      new TH1D("hHFPlusLeading_NoCut", "HF+ Leading Energy (ZB events);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_NoCut =
      new TH1D("hHFMinusLeading_NoCut", "HF- Leading Energy (ZB events);Energy (GeV);Events", 400, 0, 50);

  TH1D *hHFPlusLeading_0n0n_ZB_OfflineSel =
      new TH1D("hHFPlusLeading_0n0n_ZB_OfflineSel",
               "HF+ Leading Energy (ZB + 0nA offline selection);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0n0n_ZB_OfflineSel =
      new TH1D("hHFMinusLeading_0n0n_ZB_OfflineSel",
               "HF- Leading Energy (ZB + 0nA offline selection);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_0nAn_ZB_OfflineSel =
      new TH1D("hHFPlusLeading_0nAn_ZB_OfflineSel",
               "HF+ Leading Energy (ZB + 0nA offline selection);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0nAn_ZB_OfflineSel =
      new TH1D("hHFMinusLeading_0nAn_ZB_OfflineSel",
               "HF- Leading Energy (ZB + 0nA offline selection);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_0nXn_ZB_OfflineSel =
      new TH1D("hHFPlusLeading_0nXn_ZB_OfflineSel",
               "HF+ Leading Energy (ZB + 0nA offline selection);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0nXn_ZB_OfflineSel =
      new TH1D("hHFMinusLeading_0nXn_ZB_OfflineSel",
               "HF- Leading Energy (ZB + 0nA offline selection);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_0n0n_ZB_OfflineSel_Trig =
      new TH1D("hHFPlusLeading_0n0n_ZB_OfflineSel_Trig",
               "HF+ Leading Energy (ZB + 0nA offline selection + ADC cut);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0n0n_ZB_OfflineSel_Trig =
      new TH1D("hHFMinusLeading_0n0n_ZB_OfflineSel_Trig",
               "HF- Leading Energy (ZB + 0nA offline selection + ADC cut);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_0nAn_ZB_OfflineSel_Trig =
      new TH1D("hHFPlusLeading_0nAn_ZB_OfflineSel_Trig",
               "HF+ Leading Energy (ZB + 0nA offline selection + ADC cut);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0nAn_ZB_OfflineSel_Trig =
      new TH1D("hHFMinusLeading_0nAn_ZB_OfflineSel_Trig",
               "HF- Leading Energy (ZB + 0nA offline selection + ADC cut);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFPlusLeading_0nXn_ZB_OfflineSel_Trig =
      new TH1D("hHFPlusLeading_0nXn_ZB_OfflineSel_Trig",
               "HF+ Leading Energy (ZB + 0nA offline selection + ADC cut);Energy (GeV);Events", 400, 0, 50);
  TH1D *hHFMinusLeading_0nXn_ZB_OfflineSel_Trig =
      new TH1D("hHFMinusLeading_0nXn_ZB_OfflineSel_Trig",
               "HF- Leading Energy (ZB + 0nA offline selection + ADC cut);Energy (GeV);Events", 400, 0, 50);
  TH1D *hMult_NoCut = new TH1D("hMult_NoCut", "Track Multiplicity (ZB events);Multiplicity;Events", 400, 0, 800);
  TH1D *hMult_0nAn = new TH1D("hMult_0nAn", "Track Multiplicity (ZB events);Multiplicity;Events", 400, 0, 800);
  TH1D *hMult_0nXn = new TH1D("hMult_0nXn", "Track Multiplicity (ZB events);Multiplicity;Events", 400, 0, 800);
  TH1D *hMult_0n0n = new TH1D("hMult_0n0n", "Track Multiplicity (ZB events);Multiplicity;Events", 400, 0, 800);

  int nEntries = HiEvtTree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    HiEvtTree->GetEntry(i);
    HFADCTree->GetEntry(i);
    ZDCTree->GetEntry(i);
    TriggerTree->GetEntry(i);
    SkimTree->GetEntry(i);
    TrackTree->GetEntry(i);

    if (TriggerChoice == 1 && HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v1 != 1) {
      continue; // Apply trigger cut if specified
    }
    // if (pprimaryVertexFilter != 1)
    //   continue; // Apply primary vertex filter
    int nHPTrk = 0;
    for (int iTrk = 0; iTrk < TrackPt->size(); iTrk++) {
      if (TrackHP->at(iTrk) && abs(TrackEta->at(iTrk)) < 2.4 && TrackPt->at(iTrk) > 0.4) {
        nHPTrk++;
      }
    }
    hMult_NoCut->Fill(nHPTrk);

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
    ZDCSumPlusHist->Fill(ZDCSumPlus);
    ZDCSumMinusHist->Fill(ZDCSumMinus);

    if (ZDCSumPlus < ZDCThreshold || ZDCSumMinus < ZDCThreshold) {
      hMult_0nAn->Fill(nHPTrk);
    }
    if (ZDCSumPlus < ZDCThreshold) {
      ZDCSumPlus_0nHist->Fill(ZDCSumPlus);
      hHFPlusLeading_0nAn_ZB_OfflineSel->Fill(HFPlusLeading);
      if (HFPlusADC <= ADCThreshold) {
        hHFPlusLeading_0nAn_ZB_OfflineSel_Trig->Fill(HFPlusLeading);
      } else {
        HFPlusEtaPhi->Fill(HFPlusADC_eta, HFPlusADC_phi);
      }
    }
    if (ZDCSumMinus < ZDCThreshold) {
      ZDCSumMinus_0nHist->Fill(ZDCSumMinus);
      hHFMinusLeading_0nAn_ZB_OfflineSel->Fill(HFMinusLeading);
      if (HFMinusADC <= ADCThreshold) {
        hHFMinusLeading_0nAn_ZB_OfflineSel_Trig->Fill(HFMinusLeading);
      } else {
        HFMinusEtaPhi->Fill(HFMinusADC_eta, HFMinusADC_phi);
      }
    }
    if (ZDCSumPlus < ZDCThreshold && ZDCSumMinus >= ZDCThreshold) {
      hMult_0nXn->Fill(nHPTrk);
      hHFPlusLeading_0nXn_ZB_OfflineSel->Fill(HFPlusLeading);
      if (HFPlusADC <= ADCThreshold) {
        hHFPlusLeading_0nXn_ZB_OfflineSel_Trig->Fill(HFPlusLeading);
      }
    }
    if (ZDCSumMinus < ZDCThreshold && ZDCSumPlus >= ZDCThreshold) {
      hMult_0nXn->Fill(nHPTrk);
      hHFMinusLeading_0nXn_ZB_OfflineSel->Fill(HFMinusLeading);
      if (HFMinusADC <= ADCThreshold) {
        hHFMinusLeading_0nXn_ZB_OfflineSel_Trig->Fill(HFMinusLeading);
      }
    }

    if (ZDCSumPlus < ZDCThreshold && ZDCSumMinus < ZDCThreshold) {
      hMult_0n0n->Fill(nHPTrk);
      hHFPlusLeading_0n0n_ZB_OfflineSel->Fill(HFPlusLeading);
      hHFMinusLeading_0n0n_ZB_OfflineSel->Fill(HFMinusLeading);
      if (HFPlusADC <= ADCThreshold) {
        hHFPlusLeading_0n0n_ZB_OfflineSel_Trig->Fill(HFPlusLeading);
      }
      if (HFMinusADC <= ADCThreshold) {
        hHFMinusLeading_0n0n_ZB_OfflineSel_Trig->Fill(HFMinusLeading);
      }
    }
  }

  TFile *OutputFile = new TFile(OutputFileName.c_str(), "RECREATE");
  OutputFile->cd();
  ZDCSumPlusHist->Write();
  ZDCSumMinusHist->Write();
  ZDCSumPlus_0nHist->Write();
  ZDCSumMinus_0nHist->Write();
  hHFPlusLeading_Xn0nAll->Write();
  hHFMinusLeading_Xn0nAll->Write();
  hHFPlusLeading_Xn0nOneSide->Write();
  hHFMinusLeading_Xn0nOneSide->Write();
  hHFPlusLeading_0n0n->Write();
  hHFMinusLeading_0n0n->Write();
  hHFPlusLeading_NoCut->Write();
  hHFMinusLeading_NoCut->Write();
  hHFMinusLeading_0n0n_ZB_OfflineSel->Write();
  hHFPlusLeading_0n0n_ZB_OfflineSel->Write();
  hHFPlusLeading_0nAn_ZB_OfflineSel->Write();
  hHFMinusLeading_0nAn_ZB_OfflineSel->Write();
  hHFPlusLeading_0nXn_ZB_OfflineSel->Write();
  hHFMinusLeading_0nXn_ZB_OfflineSel->Write();
  hHFPlusLeading_0n0n_ZB_OfflineSel_Trig->Write();
  hHFMinusLeading_0n0n_ZB_OfflineSel_Trig->Write();
  hHFPlusLeading_0nAn_ZB_OfflineSel_Trig->Write();
  hHFMinusLeading_0nAn_ZB_OfflineSel_Trig->Write();
  hHFPlusLeading_0nXn_ZB_OfflineSel_Trig->Write();
  hHFMinusLeading_0nXn_ZB_OfflineSel_Trig->Write();
  hMult_NoCut->Write();
  hMult_0nAn->Write();
  hMult_0nXn->Write();
  hMult_0n0n->Write();
  HFPlusEtaPhi->Write();
  HFMinusEtaPhi->Write();

  OutputFile->Close();

  return 0;
}
