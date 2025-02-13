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
  if (par.IsData)
  {
    if (par.TriggerChoice == 1 && b->isL1ZDCOr == false)
      return false;
    if (par.TriggerChoice == 2 && b->isL1ZDCXORJet8 == false)
      return false;
  }

  if (b->selectedBkgFilter == false || b->selectedVtxFilter == false)
    return false;

  if (par.DoSystRapGap==-1)
  {
    // alternative (loose) rapidity gap selection
    if (par.IsGammaN && b->gammaN_EThreshSyst15() == false)
      return false;
    if (!par.IsGammaN && b->Ngamma_EThreshSyst15() == false)
      return false;
  }
  else if (par.DoSystRapGap==1)
  {
    // alternative (tight) rapidity gap selection
    if (par.IsGammaN && b->gammaN_EThreshSyst5p5() == false)
      return false;
    if (!par.IsGammaN && b->Ngamma_EThreshSyst5p5() == false)
      return false;
  } 
  else if (par.DoSystRapGap > 9) {
    // Custom rapidity gap threshold decision
    if (par.IsGammaN && b->gammaN_EThreshCustom(((float)par.DoSystRapGap)/10.) == false)
      return false;
    if (!par.IsGammaN && b->Ngamma_EThreshCustom(((float)par.DoSystRapGap)/10.) == false)
      return false;
  } 
  else
  {
    // nominal rapidity gap selection
    if (par.IsGammaN && (b->ZDCgammaN && b->gapgammaN) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && b->gapNgamma) == false)
      return false;
  }

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
  TH1D *hDenEvtEff;
  TH1D *hNumEvtEff;
  TH1D *hRatioEvtEff;
  TH1D *hDenDEff;
  TH1D *hNumDEff;
  TH1D *hRatioDEff;

  // ================================================================ //
  // List of variables used for MC-data comparisons 

  // Event level variables:
  TH2D *hHFEmaxPlus_vs_EvtMult; // gap energy vs event multiplicity
  TH2D *hHFEmaxMinus_vs_EvtMult; // gap energy vs event multiplicity
  TH1D *hnVtx, *hVX, *hVY, *hVZ; // number of verticies and best vertex positions

  // Track level variables:
  THnSparseD *hSparseTrackInfo;
  TH2D *htrkPt_vs_trkEta;

  // D level variables: 

  TH1D *hDchi2cl, *hDalpha, *hDdtheta, *hDsvpvDistance, *hDsvpvDisErr, *hDsvpvSig;
  TH2D *hDtrk1Pt_vs_Dtrk2Pt, *hDpt_vs_Dy;

  // TODO: Add track and jet level variables when present in the skims
  // TODO: Add gen level plots

  // ================================================================ //

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

      // Check if the event is a signal MC event, a.k.a., having at least one gen-level D candidate that falls into the (pt,y) bin of interest
      bool isSigMCEvt = false;
      if (!par.IsData){
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
          isSigMCEvt = true;
        }
      }

      if (!par.IsData && isSigMCEvt) hDenEvtEff->Fill(1);

      // Check if the event passes the selection criteria
      if (eventSelection(MDzeroUPC, par)) {
        if (!par.IsData && isSigMCEvt) hNumEvtEff->Fill(1);

        bool passedDSel = false;

        for (unsigned long j = 0; j < MDzeroUPC->Dalpha->size(); j++) {
          if (MDzeroUPC->Dpt->at(j) < par.MinDzeroPT)
            continue;
          if (MDzeroUPC->Dpt->at(j) > par.MaxDzeroPT)
            continue;
          if (MDzeroUPC->Dy->at(j) < par.MinDzeroY)
            continue;
          if (MDzeroUPC->Dy->at(j) > par.MaxDzeroY)
            continue;
          // if (MDzeroUPC->DpassCut23PAS->at(j) == false)
          if (par.DoSystD==0 && MDzeroUPC->DpassCut23PAS->at(j) == false) continue;
          // if (par.DoSystD==0 && MDzeroUPC->DpassCut23LowPt->at(j) == false) continue;
          if (par.DoSystD==1 && MDzeroUPC->DpassCut23PASSystDsvpvSig->at(j) == false) continue;
          if (par.DoSystD==2 && MDzeroUPC->DpassCut23PASSystDtrkPt->at(j) == false) continue;
          if (par.DoSystD==3 && MDzeroUPC->DpassCut23PASSystDalpha->at(j) == false) continue;
          if (par.DoSystD==4 && MDzeroUPC->DpassCut23PASSystDchi2cl->at(j) == false) continue;

          passedDSel = true;
          hDmass->Fill((*MDzeroUPC->Dmass)[j]);
          fillRecoDLevelHistograms(j);
          if (!par.IsData) {
            nt->Fill((*MDzeroUPC->Dmass)[j], (*MDzeroUPC->Dgen)[j]);
            if (MDzeroUPC->Dgen->at(j) == 23333) {
              hNumDEff->Fill(1);
            }
          } else
            nt->Fill((*MDzeroUPC->Dmass)[j], 0);
        } // end of reco-level Dzero loop

        if (!par.IsData && isSigMCEvt) {
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
            hDenDEff->Fill(1);
          } // end of gen-level Dzero loop
        }   // end of gen-level Dzero loop

        if (passedDSel) {
          fillEventLevelHistograms();
        }
      }   // end of event selection
    }     // end of event loop
  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
    smartWrite(hDmass);
    hRatioEvtEff->Divide(hNumEvtEff, hDenEvtEff, 1, 1, "B");
    hRatioDEff->Divide(hNumDEff, hDenDEff, 1, 1, "B");
    hDenEvtEff->Write();
    hNumEvtEff->Write();
    hRatioEvtEff->Write();
    hDenDEff->Write();
    hNumDEff->Write();
    hRatioDEff->Write();
    smartWrite(nt);

    // Data-MC histograms
    smartWrite(hHFEmaxMinus_vs_EvtMult);
    smartWrite(hHFEmaxPlus_vs_EvtMult);
    smartWrite(hnVtx);
    smartWrite(hVX);
    smartWrite(hVY);
    smartWrite(hVZ);

    smartWrite(htrkPt_vs_trkEta);

    smartWrite(hDalpha);
    smartWrite(hDchi2cl);
    smartWrite(hDdtheta);
    smartWrite(hDsvpvDisErr);
    smartWrite(hDsvpvDistance);
    smartWrite(hDsvpvSig);
    smartWrite(hDtrk1Pt_vs_Dtrk2Pt);
    smartWrite(hDpt_vs_Dy);
  }

private:
  void deleteHistograms() {
    delete hDmass;
    delete hDenEvtEff;
    delete hNumEvtEff;
    delete hRatioEvtEff;
    delete hDenDEff;
    delete hNumDEff;
    delete hRatioDEff;

    // Data-MC histograms
    delete hHFEmaxMinus_vs_EvtMult;
    delete hHFEmaxPlus_vs_EvtMult;
    delete hnVtx;
    delete hVX;
    delete hVY;
    delete hVZ;

    if (htrkPt_vs_trkEta != nullptr) {
      delete htrkPt_vs_trkEta;
    }

    delete hDalpha;
    delete hDchi2cl;
    delete hDdtheta;
    delete hDsvpvDisErr;
    delete hDsvpvDistance;
    delete hDsvpvSig;
    delete hDtrk1Pt_vs_Dtrk2Pt;
    delete hDpt_vs_Dy;
  }

  void initiateHistograms(Parameters &par) {
    hDmass = new TH1D(Form("hDmass%s", title.c_str()), "", 60, 1.7, 2.0);
    hDenEvtEff = new TH1D(Form("hDenEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumEvtEff = new TH1D(Form("hNumEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioEvtEff = (TH1D*) hNumEvtEff->Clone(Form("hRatioEvtEff%s", title.c_str()));
    hDenDEff = new TH1D(Form("hDenDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumDEff = new TH1D(Form("hNumDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioDEff = (TH1D*) hNumDEff->Clone(Form("hRatioDEff%s",title.c_str()));

    // Data-MC histograms
    hHFEmaxMinus_vs_EvtMult = new TH2D(Form("hHFEmaxMinus_vs_EvtMult%s", title.c_str()), "", 80, 0, 20, 200, 0, 600);
    hHFEmaxPlus_vs_EvtMult = new TH2D(Form("hHFEmaxPlus_vs_EvtMult%s", title.c_str()), "", 80, 0, 20, 200, 0, 600);
    hnVtx = new TH1D(Form("hnVtx%s", title.c_str()), "", 6, -0.5, 5.5);
    hVX = new TH1D(Form("hVX%s", title.c_str()), "", 50, -5, 5);
    hVY = new TH1D(Form("hVY%s", title.c_str()), "", 50, -5, 5);
    hVZ = new TH1D(Form("hVZ%s", title.c_str()), "", 100, -15, 15);

    hSparseTrackInfo = (THnSparseD*) inf->Get("hSparseTrackInfo");
    htrkPt_vs_trkEta = nullptr;
    // If present in the file, project to trk variables
    if (hSparseTrackInfo != nullptr) {
      hSparseTrackInfo->GetAxis(2)->SetRange(hSparseTrackInfo->GetAxis(2)->FindBin(par.MinDzeroPT), hSparseTrackInfo->GetAxis(2)->FindBin(par.MaxDzeroPT));
      hSparseTrackInfo->GetAxis(3)->SetRange(hSparseTrackInfo->GetAxis(3)->FindBin(par.MinDzeroY), hSparseTrackInfo->GetAxis(3)->FindBin(par.MaxDzeroY));
      htrkPt_vs_trkEta = (TH2D*) hSparseTrackInfo->Projection(0, 1, "E");
      htrkPt_vs_trkEta->SetName("htrkPt_vs_trkEta");
    }

    hDalpha = new TH1D(Form("hDalpha%s",title.c_str()), "", 100, 0, 3.15);
    hDchi2cl = new TH1D(Form("hDchi2cl%s",title.c_str()), "", 100, 0, 1);
    hDdtheta = new TH1D(Form("hDdtheta%s",title.c_str()), "", 100, 0, 3.15);
    hDsvpvDisErr = new TH1D(Form("hDsvpvDisErr%s",title.c_str()), "", 100, 0, 5);
    hDsvpvDistance = new TH1D(Form("hDsvpvDistance%s",title.c_str()), "", 100, 0, 5);
    hDsvpvSig = new TH1D(Form("hDsvpvSig%s",title.c_str()), "", 150, 0, 15);
    hDtrk1Pt_vs_Dtrk2Pt = new TH2D(Form("hDtrk1Pt_vs_Dtrk2Pt%s",title.c_str()), "", 100, 0, 16, 100, 0, 16);
    hDpt_vs_Dy = new TH2D(Form("hDpt_vs_Dy%s",title.c_str()), "", 120, 0, 12, 40, -2, 2);


    hDmass->Sumw2();
    hDenEvtEff->Sumw2();
    hNumEvtEff->Sumw2();
    hRatioEvtEff->Sumw2();
    hDenDEff->Sumw2();
    hNumDEff->Sumw2();
    hRatioDEff->Sumw2();

    // Data-MC histograms
    hHFEmaxMinus_vs_EvtMult->Sumw2();
    hHFEmaxPlus_vs_EvtMult->Sumw2();
    hnVtx->Sumw2();
    hVX->Sumw2();
    hVY->Sumw2();
    hVZ->Sumw2();

    if (htrkPt_vs_trkEta != nullptr) {
      htrkPt_vs_trkEta->Sumw2();
    }

    hDalpha->Sumw2();
    hDchi2cl->Sumw2();
    hDdtheta->Sumw2();
    hDsvpvDisErr->Sumw2();
    hDsvpvDistance->Sumw2();
    hDsvpvSig->Sumw2();
    hDtrk1Pt_vs_Dtrk2Pt->Sumw2();
    hDpt_vs_Dy->Sumw2();
  }

  void fillEventLevelHistograms() {
    // TODO: Add weighted filling when wheigthing strategy is finalized
    hHFEmaxMinus_vs_EvtMult->Fill(MDzeroUPC->HFEMaxMinus, MDzeroUPC->nTrackInAcceptanceHP);
    hHFEmaxPlus_vs_EvtMult->Fill(MDzeroUPC->HFEMaxPlus, MDzeroUPC->nTrackInAcceptanceHP);
    hnVtx->Fill(MDzeroUPC->nVtx);
    hVX->Fill(MDzeroUPC->VX);
    hVY->Fill(MDzeroUPC->VY);
    hVZ->Fill(MDzeroUPC->VZ);
  }

  void fillRecoDLevelHistograms(int Dindex) {
    hDalpha->Fill(MDzeroUPC->Dalpha->at(Dindex));
    hDchi2cl->Fill(MDzeroUPC->Dchi2cl->at(Dindex));
    hDdtheta->Fill(MDzeroUPC->Ddtheta->at(Dindex));
    hDsvpvDisErr->Fill(MDzeroUPC->DsvpvDisErr->at(Dindex));
    hDsvpvDistance->Fill(MDzeroUPC->DsvpvDistance->at(Dindex));
    hDsvpvSig->Fill(MDzeroUPC->DsvpvDistance->at(Dindex) / MDzeroUPC->DsvpvDisErr->at(Dindex));
    hDtrk1Pt_vs_Dtrk2Pt->Fill(MDzeroUPC->Dtrk1Pt->at(Dindex), MDzeroUPC->Dtrk2Pt->at(Dindex));
    hDpt_vs_Dy->Fill(MDzeroUPC->Dpt->at(Dindex), MDzeroUPC->Dy->at(Dindex));
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
  int DoSystRapGap = CL.GetInt("DoSystRapGap", 0);   // Systematic study: apply the alternative event selections
                                                     // 0 = nominal, 1 = tight, -1: loose
                                                     // 9 < DoSystRapGap: use custom HF energy threshold, the threshold value will be DoSystRapGap/10.
  int DoSystD = CL.GetInt("DoSystD", 0);             // Systematic study: apply the alternative D selections
                                                     // 0 = nominal, 1 = Dsvpv variation, 2: DtrkPt variation
                                                     // 3 = Dalpha variation, 4: Dchi2cl variation

  bool IsData = CL.GetBool("IsData", 0);              // Data or MC
  Parameters par(MinDzeroPT, MaxDzeroPT, MinDzeroY, MaxDzeroY, IsGammaN, TriggerChoice, IsData, scaleFactor,
                 DoSystRapGap, DoSystD);
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
