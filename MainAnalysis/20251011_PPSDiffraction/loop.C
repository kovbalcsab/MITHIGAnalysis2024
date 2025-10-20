#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm> // for std::min
#include <iostream>
#include <vector>
void loop() {

  gStyle->SetOptStat(0);
  // Open the ROOT file
  TFile *f = TFile::Open("/data00/bakovacs/OOsamples/Skims/20250703_pO_PhysicsIonPhysics1_393952_full.root");
  // TFile *f = TFile::Open("20250703_pO_PhysicsIonPhysics0_393952_full.root");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open file!\n";
    return;
  }

  TTree *t = nullptr;
  f->GetObject("Tree", t); // ← update if your tree is named differently
  if (!t) {
    std::cerr << "Tree not found!\n";
    f->Close();
    return;
  }

  // x axis = number of tracks, y axis = number of events
  TH1F *hNtrNoCond = new TH1F("hNtrNoCond", "; n. HP tracks |#eta|<2.4, p_{T}>0.3 GeV; Entries", 1000, -0.5, 999.5);
  TH1F *hNtrPPS = new TH1F("hNtrPPS", "; n. HP tracks |#eta|<2.4, p_{T}>0.3 GeV; Entries", 1000, -0.5, 999.5);
  TH1F *htrkPhiNoCond = new TH1F("htrkPhiNoCond", "; track #phi; Entries", 200, -3.15, 3.15);
  TH1F *htrkPhiPPS = new TH1F("htrkPhiPPS", "; track #phi; Entries", 200, -3.15, 3.15);
  TH2F *hNtrvsFSC3TOT_chargefCNoCond =
      new TH2F("hNtrvsFSC3TOT_chargefCNoCond", "; n. HP tracks |#eta|<2.4, p_{T}>0.3 GeV; FSC3 sum [fC]", 100, 0, 1000, 4000, 0, 400000 * 4);
  TH2F *hNtrvsFSC3TOT_chargefCPPS =
      new TH2F("hNtrvsFSC3TOT_chargefCPPS", "; n. HP tracks |#eta|<2.4, p_{T}>0.3 GeV; FSC3 sum [fC]", 100, 0, 1000, 4000, 0, 400000 * 4);
  TH2F *hZDCpvsZDCmNoCond = new TH2F("hZDCpvsZDCmNoCond", "; ZDC+; ZDC-", 100, 0, 8000, 100, 0, 8000);
  TH2F *hZDCpvsZDCmPPS = new TH2F("hZDCpvsZDCmPPS", "; ZDC+; ZDC-", 100, 0, 8000, 100, 0, 8000);
  TH1F *hZDCmNoCond = new TH1F("hZDCmNoCond", "; ZDC-; Entries", 100, 0, 8000);
  TH1F *hZDCmPPS = new TH1F("hZDCmPPS", "; ZDC-; Entries", 100, 0, 8000);
  TH1F *hZDCpNoCond = new TH1F("hZDCpNoCond", "; ZDC+; Entries", 100, 0, 8000);
  TH1F *hZDCpPPS = new TH1F("hZDCpPPS", "; ZDC+; Entries", 100, 0, 8000);
  TH1F *hFSC2bottomM_adcNoCond = new TH1F("hFSC2bottomM_adcNoCond", "; FSC2bottomM_adcNoCond; Entries", 100, 0, 256);
  TH1F *hFSC2bottomM_adcPPS = new TH1F("hFSC2bottomM_adcPPS", "; FSC2bottomM_adcPPS; Entries", 100, 0, 256);
  TH1F *hFSC3bottomleftM_adcNoCond =
      new TH1F("hFSC3bottomleftM_adcNoCond", "; FSC3bottomleftM_adcNoCond; Entries", 100, 0, 256);
  TH1F *hFSC3bottomleftM_adcPPS = new TH1F("hFSC3bottomleftM_adcPPS", "; FSC3bottomleftM_adcPPS; Entries", 100, 0, 256);
  TH1F *hFSC2bottomM_chargefCNoCond =
      new TH1F("hFSC2bottomM_chargefCNoCond", "; FSC2bottomM_chargefCNoCond; Entries", 10000, 0, 400000);
  TH1F *hFSC2bottomM_chargefCPPS =
      new TH1F("hFSC2bottomM_chargefCPPS", "; FSC2bottomM_chargefCPPS; Entries", 10000, 0, 400000);
  TH1F *hFSC3bottomleftM_chargefCNoCond =
      new TH1F("hFSC3bottomleftM_chargefCNoCond", "; FSC3bottomleftM_chargefCNoCond; Entries", 10000, 0, 400000);
  TH1F *hFSC3bottomleftM_chargefCPPS =
      new TH1F("hFSC3bottomleftM_chargefCPPS", "; FSC3bottomleftM_chargefCPPS; Entries", 10000, 0, 400000);
  TH2F *hFSC3TOT_adcvsZDCmNoCond = new TH2F("hFSC3TOT_adcvsZDCmNoCond", "; ZDC-; FSC3 sum [ADC]", 100, 0, 8000, 400, 0, 256*4);
  TH2F *hFSC3TOT_adcvsZDCmPPS = new TH2F("hFSC3TOT_adcvsZDCmPPS", "; ZDC-; FSC3 sum [ADC]", 100, 0, 8000, 400, 0, 256*4);
  TH2F *hFSC3TOT_chargefCvsZDCmNoCond = new TH2F("hFSC3TOT_chargefCvsZDCmNoCond", "; ZDC-; FSC3 sum [fC]", 100, 0, 8000, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_chargefCvsZDCmPPS = new TH2F("hFSC3TOT_chargefCvsZDCmPPS", "; ZDC-; FSC3 sum [fC]", 100, 0, 8000, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_adcvsFSC2TOT_adcNoCond = new TH2F("hFSC3TOT_adcvsFSC2TOT_adcNoCond", "; FSC2 sum [ADC]; FSC3 sum [ADC]", 200, 0, 256*2, 400, 0, 256*4);
  TH2F *hFSC3TOT_adcvsFSC2TOT_adcPPS = new TH2F("hFSC3TOT_adcvsFSC2TOT_adcPPS", "; FSC2 sum [ADC]; FSC3 sum [ADC]", 200, 0, 256*2, 400, 0, 256*4);
  TH2F *hFSC3TOT_chargefCvsFSC2TOT_chargefCNoCond = new TH2F("hFSC3TOT_chargefCvsFSC2TOT_chargefCNoCond", "; FSC2 sum [fC]; FSC3 sum [fC]", 2000, 0, 400000*2, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_chargefCvsFSC2TOT_chargefCPPS = new TH2F("hFSC3TOT_chargefCvsFSC2TOT_chargefCPPS", "; FSC2 sum [fC]; FSC3 sum [fC]", 2000, 0, 400000*2, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_chargefCvsHFEMaxMinusNoCond = new TH2F("hFSC3TOT_chargefCvsHFEMaxMinusNoCond", "; HFE Max-; FSC3 sum [fC]", 400, 0, 100, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_chargefCvsHFEMaxMinusPPS = new TH2F("hFSC3TOT_chargefCvsHFEMaxMinusPPS", "; HFE Max-; FSC3 sum [fC]", 400, 0, 100, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_chargefCvshiHFMinusNoCond = new TH2F("hFSC3TOT_chargefCvshiHFMinusNoCond", "; hiHF-; FSC3 sum [fC]", 200, 0, 200, 4000, 0, 400000*4);
  TH2F *hFSC3TOT_chargefCvshiHFMinusPPS = new TH2F("hFSC3TOT_chargefCvshiHFMinusPPS", ";; hiHF-; FSC3 sum [fC]", 200, 0, 200, 4000, 0, 400000*4);

  TH1F *hMissingptAngleNoCond = new TH1F("hMissingptAngleNoCond", "; Angle of missing p_{T}; Entries", 100, -3.15, 3.15);
  TH1F *hMissingptAnglePPS = new TH1F("hMissingptAnglePPS", "; Angle of missing p_{T}; Entries", 100, -3.15, 3.15);
  TH2F *PPS0xvsy = new TH2F("PPS0xvsy", "; PPSStation0M_x; PPSStation0M_y", 150, -10, 20, 150, -10, 20);
  TH2F *PPS2xvsy = new TH2F("PPS2xvsy", "; PPSStation2M_x; PPSStation2M_y", 150, -10, 20, 150, -10, 20);
  TH1F *h2Ddistance = new TH1F("h2Ddistance", "; Distance; Entries", 100, 0, 20);
  TH2I *hPPS0vsPPS2Size = new TH2I("hPPS0vsPPS2Size", "; PPSStation0M size; PPSStation2M size", 20, 0, 20, 20, 0, 20);
  bool passBaselineEventSelection = false;
  float VZ;
  float VZError;
  float ZDCsumMinus, ZDCsumPlus;
  float HFEMaxMinus, HFEMaxPlus, hiHFMinus_pf;
  bool HLT_OxyZeroBias_v1 = false;

  std::vector<float> *trkPt = nullptr;
  std::vector<float> *trkEta = nullptr;
  std::vector<float> *trkPhi = nullptr;
  std::vector<bool> *highPurity = nullptr;

  std::vector<float> *PPSStation0M_x = nullptr;
  std::vector<float> *PPSStation0M_y = nullptr;
  std::vector<float> *PPSStation2M_x = nullptr;
  std::vector<float> *PPSStation2M_y = nullptr;
  std::vector<int> *FSC2topM_adc = nullptr;
  std::vector<int> *FSC2bottomM_adc = nullptr;
  std::vector<int> *FSC3bottomleftM_adc = nullptr;
  std::vector<int> *FSC3bottomrightM_adc = nullptr;
  std::vector<int> *FSC3topleftM_adc = nullptr;
  std::vector<int> *FSC3toprightM_adc = nullptr;
  std::vector<float> *FSC2topM_chargefC = nullptr;
  std::vector<float> *FSC2bottomM_chargefC = nullptr;
  std::vector<float> *FSC3bottomleftM_chargefC = nullptr;
  std::vector<float> *FSC3bottomrightM_chargefC = nullptr;
  std::vector<float> *FSC3topleftM_chargefC = nullptr;
  std::vector<float> *FSC3toprightM_chargefC = nullptr;

  t->SetBranchAddress("VZ", &VZ);
  t->SetBranchAddress("VZError", &VZError);
  t->SetBranchAddress("ZDCsumMinus", &ZDCsumMinus);
  t->SetBranchAddress("ZDCsumPlus", &ZDCsumPlus);
  t->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
  t->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
  t->SetBranchAddress("hiHFMinus_pf", &hiHFMinus_pf);

  t->SetBranchAddress("HLT_OxyZeroBias_v1", &HLT_OxyZeroBias_v1);

  t->SetBranchAddress("passBaselineEventSelection", &passBaselineEventSelection);
  t->SetBranchAddress("trkPt", &trkPt);
  t->SetBranchAddress("trkEta", &trkEta);
  t->SetBranchAddress("trkPhi", &trkPhi);
  t->SetBranchAddress("highPurity", &highPurity);

  t->SetBranchAddress("PPSStation0M_x", &PPSStation0M_x);
  t->SetBranchAddress("PPSStation0M_y", &PPSStation0M_y);
  t->SetBranchAddress("PPSStation2M_x", &PPSStation2M_x);
  t->SetBranchAddress("PPSStation2M_y", &PPSStation2M_y);
  t->SetBranchAddress("FSC2topM_adc", &FSC2topM_adc);
  t->SetBranchAddress("FSC2bottomM_adc", &FSC2bottomM_adc);
  t->SetBranchAddress("FSC3bottomleftM_adc", &FSC3bottomleftM_adc);
  t->SetBranchAddress("FSC3bottomrightM_adc", &FSC3bottomrightM_adc);
  t->SetBranchAddress("FSC3topleftM_adc", &FSC3topleftM_adc);
  t->SetBranchAddress("FSC3toprightM_adc", &FSC3toprightM_adc);
  t->SetBranchAddress("FSC2topM_chargefC", &FSC2topM_chargefC);
  t->SetBranchAddress("FSC2bottomM_chargefC", &FSC2bottomM_chargefC);
  t->SetBranchAddress("FSC3bottomleftM_chargefC", &FSC3bottomleftM_chargefC);
  t->SetBranchAddress("FSC3bottomrightM_chargefC", &FSC3bottomrightM_chargefC);
  t->SetBranchAddress("FSC3topleftM_chargefC", &FSC3topleftM_chargefC);
  t->SetBranchAddress("FSC3toprightM_chargefC", &FSC3toprightM_chargefC);

  const Long64_t nEntries = t->GetEntries();
  const Long64_t nMax = std::min<Long64_t>(nEntries, 1000000000000); // cap to 1000 as in your code
  int multHighPurity2p4Pt0p3 = 0;
  std::cerr << "Processing " << nMax << " entries...\n";

  for (Long64_t i = 0; i < nMax; ++i) {
    t->GetEntry(i);
    multHighPurity2p4Pt0p3 = 0;
    if (!passBaselineEventSelection) {
      continue;
    }
    // if (HFEMaxPlus > 8.0 || HFEMaxMinus > 8.0) {
    //   continue;
    // }

    if (!HLT_OxyZeroBias_v1) {
      continue;
    } // skip events not passing the trigger

    if (fabs(VZ) > 15.0)
      continue;

    // std::cout << "Event " << i << ": passBaselineEventSelection = " << passBaselineEventSelection << "\n";

    int size_PPSStation0M_x = PPSStation0M_x->size();
    int size_PPSStation0M_y = PPSStation0M_y->size();
    int size_PPSStation2M_x = PPSStation2M_x->size();
    int size_PPSStation2M_y = PPSStation2M_y->size();
    if (size_PPSStation0M_x != size_PPSStation0M_y) {
      std::cerr << "Warning: PPSStation0M_x and PPSStation0M_y have different sizes!\n";
      continue;
    }
    if (size_PPSStation2M_x != size_PPSStation2M_y) {
      std::cerr << "Warning: PPSStation2M_x and PPSStation2M_y have different sizes!\n";
      continue;
    }
    hPPS0vsPPS2Size->Fill(size_PPSStation0M_x, size_PPSStation2M_x);
    double missingpt_x = 0.0;
    double missingpt_y = 0.0;

    for (size_t j = 0; j < trkPt->size(); ++j) {
      if (highPurity->at(j)) {
        if (trkPt->at(j) > 0.3 && fabs(trkEta->at(j)) < 2.4) {
          multHighPurity2p4Pt0p3++;
          missingpt_x += trkPt->at(j) * cos(trkPhi->at(j));
          missingpt_y += trkPt->at(j) * sin(trkPhi->at(j));
          htrkPhiNoCond->Fill(trkPhi->at(j), trkPt->at(j));
          if (size_PPSStation0M_x > 0 && size_PPSStation0M_y > 0 && size_PPSStation2M_x > 0 && size_PPSStation2M_y > 0) {
            htrkPhiPPS->Fill(trkPhi->at(j), trkPt->at(j));
          }
        }
      }

      // std::cout << "  Track " << j << ": Pt = " << trkPt->at(j) << ", Eta = " << trkEta->at(j) << ", Phi = " <<
      // trkPhi->at(j) << "\n";
    }
    
    int FSC3TOT_adc = FSC3bottomleftM_adc->at(2) + FSC3bottomrightM_adc->at(2) + FSC3topleftM_adc->at(2) +
                      FSC3toprightM_adc->at(2);
    float FSC3TOT_chargefC = FSC3bottomleftM_chargefC->at(2) + FSC3bottomrightM_chargefC->at(2) +
                              FSC3topleftM_chargefC->at(2) + FSC3toprightM_chargefC->at(2);
    int FSC2TOT_adc = FSC2topM_adc->at(2) + FSC2bottomM_adc->at(2);
    float FSC2TOT_chargefC = FSC2topM_chargefC->at(2) + FSC2bottomM_chargefC->at(2);
    if (size_PPSStation0M_x > 0 && size_PPSStation0M_y > 0 && size_PPSStation2M_x > 0 && size_PPSStation2M_y > 0) {
      hNtrPPS->Fill(multHighPurity2p4Pt0p3);
      hZDCpvsZDCmPPS->Fill(ZDCsumPlus, ZDCsumMinus);
      hZDCmPPS->Fill(ZDCsumMinus);
      hZDCpPPS->Fill(ZDCsumPlus);
      hFSC2bottomM_adcPPS->Fill(FSC2bottomM_adc->at(2));
      hFSC2bottomM_chargefCPPS->Fill(FSC2bottomM_chargefC->at(2));
      hFSC3bottomleftM_adcPPS->Fill(FSC3bottomleftM_adc->at(2));
      hFSC3bottomleftM_chargefCPPS->Fill(FSC3bottomleftM_chargefC->at(2));
      hFSC3TOT_adcvsZDCmPPS->Fill(ZDCsumMinus, FSC3TOT_adc);
      hFSC3TOT_chargefCvsZDCmPPS->Fill(ZDCsumMinus, FSC3TOT_chargefC);
      hFSC3TOT_adcvsFSC2TOT_adcPPS->Fill(FSC2TOT_adc, FSC3TOT_adc);
      hFSC3TOT_chargefCvsFSC2TOT_chargefCPPS->Fill(FSC2TOT_chargefC, FSC3TOT_chargefC);
      hNtrvsFSC3TOT_chargefCPPS->Fill(multHighPurity2p4Pt0p3, FSC3TOT_chargefC);
      PPS0xvsy->Fill(PPSStation0M_x->at(0), PPSStation0M_y->at(0));
      PPS2xvsy->Fill(PPSStation2M_x->at(0), PPSStation2M_y->at(0));
      float distance =
          sqrt(pow(PPSStation0M_x->at(0) - PPSStation2M_x->at(0), 2) + pow(PPSStation0M_y->at(0) - PPSStation2M_y->at(0), 2));
      h2Ddistance->Fill(distance);
      hMissingptAnglePPS->Fill(atan2(missingpt_y, missingpt_x));
      hFSC3TOT_chargefCvsHFEMaxMinusPPS->Fill(HFEMaxMinus, FSC3TOT_chargefC);
      hFSC3TOT_chargefCvshiHFMinusPPS->Fill(hiHFMinus_pf, FSC3TOT_chargefC);
    }
    hNtrNoCond->Fill(multHighPurity2p4Pt0p3);
    hZDCpvsZDCmNoCond->Fill(ZDCsumPlus, ZDCsumMinus);
    hZDCmNoCond->Fill(ZDCsumMinus);
    hZDCpNoCond->Fill(ZDCsumPlus);
    hFSC2bottomM_adcNoCond->Fill(FSC2bottomM_adc->at(2));
    hFSC2bottomM_chargefCNoCond->Fill(FSC2bottomM_chargefC->at(2));
    hFSC3bottomleftM_adcNoCond->Fill(FSC3bottomleftM_adc->at(2));
    hFSC3bottomleftM_chargefCNoCond->Fill(FSC3bottomleftM_chargefC->at(2));
    hFSC3TOT_adcvsZDCmNoCond->Fill(ZDCsumMinus, FSC3TOT_adc);
    hFSC3TOT_chargefCvsZDCmNoCond->Fill(ZDCsumMinus, FSC3TOT_chargefC);
    hFSC3TOT_adcvsFSC2TOT_adcNoCond->Fill(FSC2TOT_adc, FSC3TOT_adc);
    hFSC3TOT_chargefCvsFSC2TOT_chargefCNoCond->Fill(FSC2TOT_chargefC, FSC3TOT_chargefC);
    hNtrvsFSC3TOT_chargefCNoCond->Fill(multHighPurity2p4Pt0p3, FSC3TOT_chargefC);
    hMissingptAngleNoCond->Fill(atan2(missingpt_y, missingpt_x));
    hFSC3TOT_chargefCvsHFEMaxMinusNoCond->Fill(HFEMaxMinus, FSC3TOT_chargefC);
    hFSC3TOT_chargefCvshiHFMinusNoCond->Fill(hiHFMinus_pf, FSC3TOT_chargefC);
  }

  TFile *outFile = TFile::Open("output.root", "RECREATE");
  hNtrNoCond->Write();
  hNtrPPS->Write();
  hZDCpvsZDCmNoCond->Write();
  hZDCpvsZDCmPPS->Write();
  hZDCmNoCond->Write();
  hZDCmPPS->Write();
  hZDCpNoCond->Write();
  hZDCpPPS->Write();
  hFSC2bottomM_adcNoCond->Write();
  hFSC2bottomM_adcPPS->Write();
  hFSC3bottomleftM_adcNoCond->Write();
  hFSC3bottomleftM_adcPPS->Write();
  hFSC2bottomM_chargefCNoCond->Write();
  hFSC2bottomM_chargefCPPS->Write();
  hFSC3bottomleftM_chargefCNoCond->Write();
  hFSC3bottomleftM_chargefCPPS->Write();
  hFSC3TOT_adcvsZDCmNoCond->Write();
  hFSC3TOT_adcvsZDCmPPS->Write();
  hFSC3TOT_chargefCvsZDCmNoCond->Write();
  hFSC3TOT_chargefCvsZDCmPPS->Write();
  hFSC3TOT_adcvsFSC2TOT_adcNoCond->Write();
  hFSC3TOT_adcvsFSC2TOT_adcPPS->Write();
  hFSC3TOT_chargefCvsFSC2TOT_chargefCNoCond->Write();
  hFSC3TOT_chargefCvsFSC2TOT_chargefCPPS->Write();
  hNtrvsFSC3TOT_chargefCNoCond->Write();
  hNtrvsFSC3TOT_chargefCPPS->Write();
  PPS0xvsy->Write();
  PPS2xvsy->Write();
  h2Ddistance->Write();
  htrkPhiNoCond->Write();
  htrkPhiPPS->Write();
  hMissingptAngleNoCond->Write();
  hMissingptAnglePPS->Write();
  hPPS0vsPPS2Size->Write();
  hFSC3TOT_chargefCvsHFEMaxMinusNoCond->Write();
  hFSC3TOT_chargefCvsHFEMaxMinusPPS->Write();
  hFSC3TOT_chargefCvshiHFMinusNoCond->Write();
  hFSC3TOT_chargefCvshiHFMinusPPS->Write();
  outFile->Close();
  f->Close();
}

int main() {
  loop();
  return 0;
}
