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
  std::vector<TH2D *> hHFEmaxPlus_vs_EvtMult; // gap energy vs event multiplicity
  std::vector<TH2D *> hHFEmaxMinus_vs_EvtMult; // gap energy vs event multiplicity
  std::vector<TH1D *> hnVtx, hVX, hVY, hVZ; // number of verticies and best vertex positions

  // Track level variables:
  THnSparseD *hSparseTrackInfo;
  TH2D *htrkPt_vs_trkEta;

  // D level variables: 

  std::vector<TH1D *> hDchi2cl, hDalpha, hDdtheta, hDsvpvDistance, hDsvpvDisErr, hDsvpvSig, hDmass;
  std::vector<TH2D *> hDtrk1Pt_vs_Dtrk2Pt, hDpt_vs_Dy;

  // Gen D level variables

  TH2D *hGDpt_GDy; 

  // TODO: Add jet level variables when present in the skims

  const double sideBandMinusUpperEdge = 1.82;
  const double sideBandPlusLowerEdge = 1.9;

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

        bool hasSignalD = false;
        double maxDpt = -999.;
        double selectedDmass = -999.;

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

          fillRecoDLevelHistograms(j);
          if (!par.IsData) {
            nt->Fill((*MDzeroUPC->Dmass)[j], (*MDzeroUPC->Dgen)[j]);
            if (MDzeroUPC->Dgen->at(j) == 23333) {
              hNumDEff->Fill(1);
            }
          } else
            nt->Fill((*MDzeroUPC->Dmass)[j], 0);

          if (MDzeroUPC->Dmass->at(j) < sideBandMinusUpperEdge || MDzeroUPC->Dmass->at(j) > sideBandPlusLowerEdge) {
            // if the D meson corresponds to the sideband update max Dpt only if no signal D was found yet in the event
            if (!hasSignalD && MDzeroUPC->Dpt->at(j) > maxDpt) {
              maxDpt = MDzeroUPC->Dpt->at(j);
              selectedDmass = MDzeroUPC->Dmass->at(j);
            }
          } else {
            if (!hasSignalD) {
              // first signal D found always takes priority over sideband
              hasSignalD = true;
              maxDpt = MDzeroUPC->Dpt->at(j);
              selectedDmass = MDzeroUPC->Dmass->at(j);
            } else if (MDzeroUPC->Dpt->at(j) > maxDpt) {
              maxDpt = MDzeroUPC->Dpt->at(j);
              selectedDmass = MDzeroUPC->Dmass->at(j);
            }
          }
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
            fillGenDLevelHistograms(j);
          } // end of gen-level Dzero loop
        }   // end of gen-level Dzero loop

        if (maxDpt > 0) {
          fillEventLevelHistograms(selectedDmass);
        }
      }   // end of event selection
    }     // end of event loop
  }       // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
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
    smartWrite(htrkPt_vs_trkEta);
    smartWrite(hGDpt_GDy);

    for (int iHistName = 0; iHistName < Parameters::HISTO_NAMES.size(); iHistName++) {
      smartWrite(hHFEmaxMinus_vs_EvtMult[iHistName]);
      smartWrite(hHFEmaxPlus_vs_EvtMult[iHistName]);
      smartWrite(hnVtx[iHistName]);
      smartWrite(hVX[iHistName]);
      smartWrite(hVY[iHistName]);
      smartWrite(hVZ[iHistName]);


      smartWrite(hDmass[iHistName]);
      smartWrite(hDalpha[iHistName]);
      smartWrite(hDchi2cl[iHistName]);
      smartWrite(hDdtheta[iHistName]);
      smartWrite(hDsvpvDisErr[iHistName]);
      smartWrite(hDsvpvDistance[iHistName]);
      smartWrite(hDsvpvSig[iHistName]);
      smartWrite(hDtrk1Pt_vs_Dtrk2Pt[iHistName]);
      smartWrite(hDpt_vs_Dy[iHistName]);
    }
  }

private:
  void deleteHistograms() {
    delete hDenEvtEff;
    delete hNumEvtEff;
    delete hRatioEvtEff;
    delete hDenDEff;
    delete hNumDEff;
    delete hRatioDEff;

    // Data-MC histograms
    if (htrkPt_vs_trkEta != nullptr) {
      delete htrkPt_vs_trkEta;
    }
    delete hGDpt_GDy;

    for (int iHistName = 0; iHistName < Parameters::HISTO_NAMES.size(); iHistName++) {
      delete hHFEmaxMinus_vs_EvtMult[iHistName];
      delete hHFEmaxPlus_vs_EvtMult[iHistName];
      delete hnVtx[iHistName];
      delete hVX[iHistName];
      delete hVY[iHistName];
      delete hVZ[iHistName];

      delete hDmass[iHistName];
      delete hDalpha[iHistName];
      delete hDchi2cl[iHistName];
      delete hDdtheta[iHistName];
      delete hDsvpvDisErr[iHistName];
      delete hDsvpvDistance[iHistName];
      delete hDsvpvSig[iHistName];
      delete hDtrk1Pt_vs_Dtrk2Pt[iHistName];
      delete hDpt_vs_Dy[iHistName];
    }
  }

  void initiateHistograms(Parameters &par) {
    hDenEvtEff = new TH1D(Form("hDenEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumEvtEff = new TH1D(Form("hNumEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioEvtEff = (TH1D*) hNumEvtEff->Clone(Form("hRatioEvtEff%s", title.c_str()));
    hDenDEff = new TH1D(Form("hDenDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumDEff = new TH1D(Form("hNumDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioDEff = (TH1D*) hNumDEff->Clone(Form("hRatioDEff%s",title.c_str()));

    hGDpt_GDy = new TH2D(Form("hGDpt_GDy%s", title.c_str()), "", 120, 0, 12, 40, -2, 2);
    
    hSparseTrackInfo = (THnSparseD*) inf->Get("hSparseTrackInfo");
    htrkPt_vs_trkEta = nullptr;
    // If present in the file, project to trk variables
    if (hSparseTrackInfo != nullptr) {
      hSparseTrackInfo->GetAxis(2)->SetRange(hSparseTrackInfo->GetAxis(2)->FindBin(par.MinDzeroPT), hSparseTrackInfo->GetAxis(2)->FindBin(par.MaxDzeroPT));
      hSparseTrackInfo->GetAxis(3)->SetRange(hSparseTrackInfo->GetAxis(3)->FindBin(par.MinDzeroY), hSparseTrackInfo->GetAxis(3)->FindBin(par.MaxDzeroY));
      htrkPt_vs_trkEta = (TH2D*) hSparseTrackInfo->Projection(0, 1, "E");
      htrkPt_vs_trkEta->SetName("htrkPt_vs_trkEta");
    }

    hDenEvtEff->Sumw2();
    hNumEvtEff->Sumw2();
    hRatioEvtEff->Sumw2();
    hDenDEff->Sumw2();
    hNumDEff->Sumw2();
    hRatioDEff->Sumw2();

    if (htrkPt_vs_trkEta != nullptr) {
      htrkPt_vs_trkEta->Sumw2();
    }

    hGDpt_GDy->Sumw2();

    for (int iHistName = 0; iHistName < Parameters::HISTO_NAMES.size(); iHistName++) {
      // Data-MC histograms
      hHFEmaxMinus_vs_EvtMult.push_back(new TH2D(Form("hHFEmaxMinus_vs_EvtMult_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 80, 0, 20, 200, 0, 600));
      hHFEmaxPlus_vs_EvtMult.push_back(new TH2D(Form("hHFEmaxPlus_vs_EvtMult_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 80, 0, 20, 200, 0, 600));
      hnVtx.push_back(new TH1D(Form("hnVtx_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 6, -0.5, 5.5));
      hVX.push_back(new TH1D(Form("hVX_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 50, -1, 1));
      hVY.push_back(new TH1D(Form("hVY_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 50, -1, 1));
      hVZ.push_back(new TH1D(Form("hVZ_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, -15, 15));

      hDmass.push_back(new TH1D(Form("hDmass_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 60, 1.7, 2.0));
      hDalpha.push_back(new TH1D(Form("hDalpha_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, 0, 3.15));
      hDchi2cl.push_back(new TH1D(Form("hDchi2cl_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, 0, 1));
      hDdtheta.push_back(new TH1D(Form("hDdtheta_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, 0, 3.15));
      hDsvpvDisErr.push_back(new TH1D(Form("hDsvpvDisErr_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, 0, 5));
      hDsvpvDistance.push_back(new TH1D(Form("hDsvpvDistance_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, 0, 5));
      hDsvpvSig.push_back(new TH1D(Form("hDsvpvSig_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 150, 0, 15));
      hDtrk1Pt_vs_Dtrk2Pt.push_back(new TH2D(Form("hDtrk1Pt_vs_Dtrk2Pt_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 100, 0, 16, 100, 0, 16));
      hDpt_vs_Dy.push_back(new TH2D(Form("hDpt_vs_Dy_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), title.c_str()), "", 120, 0, 12, 40, -2, 2));

      // Data-MC histograms
      hHFEmaxMinus_vs_EvtMult[iHistName]->Sumw2();
      hHFEmaxPlus_vs_EvtMult[iHistName]->Sumw2();
      hnVtx[iHistName]->Sumw2();
      hVX[iHistName]->Sumw2();
      hVY[iHistName]->Sumw2();
      hVZ[iHistName]->Sumw2();

      hDmass[iHistName]->Sumw2();
      hDalpha[iHistName]->Sumw2();
      hDchi2cl[iHistName]->Sumw2();
      hDdtheta[iHistName]->Sumw2();
      hDsvpvDisErr[iHistName]->Sumw2();
      hDsvpvDistance[iHistName]->Sumw2();
      hDsvpvSig[iHistName]->Sumw2();
      hDtrk1Pt_vs_Dtrk2Pt[iHistName]->Sumw2();
      hDpt_vs_Dy[iHistName]->Sumw2();
    }
  }

  void fillEventLevelHistograms(double selectedDmass) {
    // select D mass region
    int sideBandIndex = -1;
    if (selectedDmass < sideBandMinusUpperEdge) {
      sideBandIndex = 0; // fill sideband minus histograms
    } else if (selectedDmass > sideBandPlusLowerEdge) {
      sideBandIndex = 2; // fill sideband plus histograms
    } else {
      sideBandIndex = 1; // fill signal histograms
    }

    // TODO: Add weighted filling when wheigthing strategy is finalized
    hHFEmaxMinus_vs_EvtMult[sideBandIndex]->Fill(MDzeroUPC->HFEMaxMinus, MDzeroUPC->nTrackInAcceptanceHP);
    hHFEmaxPlus_vs_EvtMult[sideBandIndex]->Fill(MDzeroUPC->HFEMaxPlus, MDzeroUPC->nTrackInAcceptanceHP);
    hnVtx[sideBandIndex]->Fill(MDzeroUPC->nVtx);
    hVX[sideBandIndex]->Fill(MDzeroUPC->VX);
    hVY[sideBandIndex]->Fill(MDzeroUPC->VY);
    hVZ[sideBandIndex]->Fill(MDzeroUPC->VZ);
  }

  void fillRecoDLevelHistograms(int Dindex) {
    // select D mass region
    int sideBandIndex = -1;
    if (MDzeroUPC->Dmass->at(Dindex) < sideBandMinusUpperEdge) {
      sideBandIndex = 0; // fill sideband minus histograms
    } else if (MDzeroUPC->Dmass->at(Dindex) > sideBandPlusLowerEdge) {
      sideBandIndex = 2; // fill sideband plus histograms
    } else {
      sideBandIndex = 1; // fill signal histograms
    }

    // Do the filling
    hDmass[sideBandIndex]->Fill(MDzeroUPC->Dmass->at(Dindex));
    hDalpha[sideBandIndex]->Fill(MDzeroUPC->Dalpha->at(Dindex));
    hDchi2cl[sideBandIndex]->Fill(MDzeroUPC->Dchi2cl->at(Dindex));
    hDdtheta[sideBandIndex]->Fill(MDzeroUPC->Ddtheta->at(Dindex));
    hDsvpvDisErr[sideBandIndex]->Fill(MDzeroUPC->DsvpvDisErr->at(Dindex));
    hDsvpvDistance[sideBandIndex]->Fill(MDzeroUPC->DsvpvDistance->at(Dindex));
    hDsvpvSig[sideBandIndex]->Fill(MDzeroUPC->DsvpvDistance->at(Dindex) / MDzeroUPC->DsvpvDisErr->at(Dindex));
    hDtrk1Pt_vs_Dtrk2Pt[sideBandIndex]->Fill(MDzeroUPC->Dtrk1Pt->at(Dindex), MDzeroUPC->Dtrk2Pt->at(Dindex));
    hDpt_vs_Dy[sideBandIndex]->Fill(MDzeroUPC->Dpt->at(Dindex), MDzeroUPC->Dy->at(Dindex));
  }

  void fillGenDLevelHistograms(int Dindex) {
    hGDpt_GDy->Fill(MDzeroUPC->Gpt->at(Dindex), MDzeroUPC->Gy->at(Dindex));
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
