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

bool eventSelection(DzeroUPCTreeMessenger *b, TH1D* hNumberOfEventsAfterCuts, bool TriggerChoice, bool IsGammaN) {
  if (TriggerChoice == 1 && b->isL1ZDCOr == false)
    return false;
  if (TriggerChoice == 2 && b->isL1ZDCXORJet8 == false)
    return false;
  hNumberOfEventsAfterCuts->Fill(1); // After trigger selection

  if (b->cscTightHalo2015Filter == false || b->selectedVtxFilter == false)
    return false;
  if (b->nVtx >= 3) return false;
  hNumberOfEventsAfterCuts->Fill(2); // After vertex selection

  if (IsGammaN && b->ZDCgammaN == false) {
    return false; // ZDC selection for GammaN
  } else if (!IsGammaN && b->ZDCNgamma == false) {
    return false; // ZDC selection for NGamma
  }
  hNumberOfEventsAfterCuts->Fill(3); // After ZDC selection

  return true;
}

class DataAnalyzer {
public:
  TFile *inf, *outf;
  string title;
  DzeroUPCTreeMessenger *MDzeroUPC;
  TH1D *hNumberOfEventsAfterCuts;
  TH1D* hDmass;
  TH1D* hDpt;
  TH1D* hDeta;
  
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
    hNumberOfEventsAfterCuts = new TH1D(Form("hNumberOfEventsAfterCuts%s", title.c_str()), "", 5, -0.5, 4.5);
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(1, "Total");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(2, "Trigger");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(3, "Vertex");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(4, "ZDC");
    hNumberOfEventsAfterCuts->GetXaxis()->SetBinLabel(5, "D selection");

    hDmass = new TH1D(Form("hDmass%s", title.c_str()), "", 100, 1.7, 2.0);
    hDpt = new TH1D(Form("hDpt%s", title.c_str()), "", 100, 0, 10);
    hDeta = new TH1D(Form("hDeta%s", title.c_str()), "", 100, -2.4, 2.4);

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

    hNumberOfEventsAfterCuts->Sumw2();
    hDmass->Sumw2();
    hDpt->Sumw2();
    hDeta->Sumw2();

    unsigned long nEntry = MDzeroUPC->GetEntries();
    ProgressBar Bar(cout, nEntry);
    Bar.SetStyle(1);

    int triggerChoice = man.GetCutParameterValue("TriggerChoice");
    if (triggerChoice != 0 && triggerChoice != 1 && triggerChoice != 2) {
      std::cerr << "Error: Invalid TriggerChoice parameter value. Expected 0, 1, or 2." << std::endl;
      return;
    }

    bool IsGammaN = man.GetCutParameterValue("IsGammaN");

    double MinDzeroPT = man.GetCutParameterValue("MinDzeroPT");
    double MaxDzeroPT = man.GetCutParameterValue("MaxDzeroPT");
    double MinDzeroY = man.GetCutParameterValue("MinDzeroY");
    double MaxDzeroY = man.GetCutParameterValue("MaxDzeroY");
    int DoSystD = man.GetCutParameterValue("DoSystD");
    if (DoSystD < 0 || DoSystD > 4) {
      std::cerr << "Error: Invalid DoSystD parameter value. Expected 0, 1, 2, 3, or 4." << std::endl;
      return;
    }

    for (unsigned long i = 0; i < nEntry; i++) {
      MDzeroUPC->GetEntry(i);
      //if (i % 1000 == 0) {
      //  Bar.Update(i);
      //  Bar.Print();
      //}
      hNumberOfEventsAfterCuts->Fill(0); // Total events

      if (!eventSelection(MDzeroUPC, hNumberOfEventsAfterCuts, triggerChoice, IsGammaN)) {
        continue;
      }

      // D selections
      bool foundRecoD = false;
      
      bool doTrkFilter = false;
      if (MDzeroUPC->Dtrk1PtErr != nullptr &&
          MDzeroUPC->Dtrk2PtErr != nullptr &&
          MDzeroUPC->Dtrk1PixelHit != nullptr &&
          MDzeroUPC->Dtrk1StripHit != nullptr &&
          MDzeroUPC->Dtrk2PixelHit != nullptr &&
          MDzeroUPC->Dtrk2StripHit != nullptr) doTrkFilter = true;
      for (unsigned long j = 0; j < MDzeroUPC->Dsize; j++) {
        if (MDzeroUPC->Dpt->at(j) < MinDzeroPT)
          continue;
        if (MDzeroUPC->Dpt->at(j) > MaxDzeroPT)
          continue;
        if (MDzeroUPC->Dy->at(j) < MinDzeroY)
          continue;
        if (MDzeroUPC->Dy->at(j) > MaxDzeroY)
          continue;
        if (DoSystD==0 && MDzeroUPC->DpassCut23PAS->at(j) == false) continue;
        if (DoSystD==1 && MDzeroUPC->DpassCut23PASSystDsvpvSig->at(j) == false) continue;
        if (DoSystD==2 && MDzeroUPC->DpassCut23PASSystDtrkPt->at(j) == false) continue;
        if (DoSystD==3 && MDzeroUPC->DpassCut23PASSystDalpha->at(j) == false) continue;
        if (DoSystD==4 && MDzeroUPC->DpassCut23PASSystDchi2cl->at(j) == false) continue;
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
        hDpt->Fill((*MDzeroUPC->Dpt)[j]);
        hDeta->Fill((*MDzeroUPC->Dy)[j]);
        foundRecoD = true;
      } // end of reco-level Dzero loop

      if (!foundRecoD) continue;
      hNumberOfEventsAfterCuts->Fill(4); // After D selection

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
    smartWrite(hDmass);
    smartWrite(hDpt);
    smartWrite(hDeta);
  }

private:
  void deleteHistograms() {
    delete hNumberOfEventsAfterCuts;
    delete hDmass;
    delete hDpt;
    delete hDeta;
  }
};

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  int TriggerChoice = CL.GetInt("TriggerChoice", 1); // 0 = no selection, 1 = isL1ZDCOr, 2 = isL1ZDCXORJet8
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  double MinDzeroPT = CL.GetDouble("MinDzeroPT", 2);  // Minimum Dzero transverse momentum threshold for Dzero selection.
  double MaxDzeroPT = CL.GetDouble("MaxDzeroPT", 5);
  double MinDzeroY = CL.GetDouble("MinDzeroY", -2);   // Minimum Dzero rapidity threshold for Dzero selection.
  double MaxDzeroY = CL.GetDouble("MaxDzeroY", +2);
  int DoSystD = CL.GetInt("DoSystD", 0); // 0 = no selection, 1 = DsvpvSig, 2 = DtrkPt, 3 = Dalpha, 4 = Dchi2cl

  string inputFile = CL.Get("Input", "mergedSample.root"); // Input file
  string outputFile = CL.Get("Output", "output.root");     // Output file

  DataAnalyzer analyzer(inputFile.c_str(), outputFile.c_str());

  GeneralInfoManager man(analyzer.outf, "InfoDir", false);
  TTimeStamp *currentTime = new TTimeStamp();

  man.AddSourceFile(inputFile, currentTime);
  man.AddCutParameter("TriggerChoice", TriggerChoice, currentTime);
  man.AddCutParameter("IsGammaN", IsGammaN, currentTime);
  man.AddCutParameter("MinDzeroPT", MinDzeroPT, currentTime);
  man.AddCutParameter("MaxDzeroPT", MaxDzeroPT, currentTime);
  man.AddCutParameter("MinDzeroY", MinDzeroY, currentTime);
  man.AddCutParameter("MaxDzeroY", MaxDzeroY, currentTime);
  man.AddCutParameter("DoSystD", DoSystD, currentTime);

  std::cout << "Parameters used for this analysis:" << std::endl;
  man.PrintInfo();

  analyzer.analyze(man);
  analyzer.writeHistograms(analyzer.outf);
  man.SaveToFile();
  std::cout << "Analysis complete. Output saved to " << outputFile << std::endl;

  return 0;
}