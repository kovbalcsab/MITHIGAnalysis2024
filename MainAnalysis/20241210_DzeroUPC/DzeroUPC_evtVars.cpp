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
#include "histogramParameters.h"
#include "parameter.h" // The parameters used in the analysis
#include "utilities.h" // Yen-Jie's random utility functions

#include "WeightHandler1D.h"
#include "WeightHandler2D.h"

#define DMASS 1.86484

// Mass ranges and bins for Dmass monitoring histograms,
// intended for cross-checking datapoints from the massfit.C routine.
// Current analysis binning:
#define DMASS_HISTMIN 1.66
#define DMASS_HISTMAX 2.26
#define DMASS_HISTNBINS 48
// Bins matched to HIN-24-003 Dmass fit bins:
#define DMASS_HISTMIN_24003 1.68
#define DMASS_HISTMAX_24003 2.05
#define DMASS_HISTNBINS_24003 74
// Bins matched to HIN-25-002 Dmass fit bins:
#define DMASS_HISTMIN_25002 1.66
#define DMASS_HISTMAX_25002 2.16
#define DMASS_HISTNBINS_25002 48

//============================================================//
// Function to check for configuration errors
//============================================================//
bool checkError(const Parameters &par) { return false; }

//======= trackSelection =====================================//
// Check if the track pass selection criteria
//============================================================//
bool dzeroSelection(DzeroUPCMicroTreeMessenger *b, Parameters par, int j) { return true; }

//======= eventSelection =====================================//
// Check if the event pass eventSelection criteria
//============================================================//
bool eventSelection(DzeroUPCMicroTreeMessenger *b, const Parameters &par, TH1D *hEventsPassedSelection) {
  hEventsPassedSelection->Fill(0); // Total events before selection

  if (par.IsData) {
    if (par.TriggerChoice == 1 && b->isL1ZDCOr == false)
      return false;
    if (par.TriggerChoice == 2 && b->isL1ZDCXORJet8 == false)
      return false;
    if (par.TriggerChoice == 3 && b->isL1ZDCXORJet12 == false)
      return false;
    if (par.TriggerChoice == 4 && b->isL1ZDCXORJet16 == false)
      return false;
  }

  hEventsPassedSelection->Fill(1); // Events after trigger selection

  if (par.IsData) {
    if ((par.DoCCFsyst == 0 && b->selectedBkgFilter == false) ||
        (par.DoCCFsyst == 1 && !(b->cscTightHalo2015Filter == true && b->ClusterCompatibilityFilter == false)) ||
        (par.DoCCFsyst == 2 && (b->cscTightHalo2015Filter == false)))
      return false;
  }

  hEventsPassedSelection->Fill(2); // Events after CCF and CSC selection

  if (b->selectedVtxFilter == false)
    return false;

  hEventsPassedSelection->Fill(3); // Events after vertex selection

  // Set HF rapidity definition
  float HFEMaxPlus;
  float HFEMaxMinus;
  HFEMaxPlus = b->HFEMaxPlus_eta5;
  HFEMaxMinus = b->HFEMaxMinus_eta5;

  if (par.DoSystRapGap == -1) {
    // alternative (loose) rapidity gap selection
    if (par.IsGammaN && (b->ZDCgammaN && HFEMaxPlus < 15) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && HFEMaxMinus < 15) == false)
      return false;
  } else if (par.DoSystRapGap == 1) {
    // alternative (tight) rapidity gap selection
    if (par.IsGammaN && (b->ZDCgammaN && HFEMaxPlus < 5.5) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && HFEMaxMinus < 5.5) == false)
      return false;
  } else if (par.DoSystRapGap > 9) {
    // Custom rapidity gap threshold decision
    if (par.IsGammaN && (b->ZDCgammaN && HFEMaxPlus < (((float)par.DoSystRapGap) / 10.)) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && HFEMaxMinus < (((float)par.DoSystRapGap) / 10.)) == false)
      return false;
  } else if (par.DoSystRapGap < -9) {
    // Custom rapidity gap threshold in rejection mode
    if (par.IsGammaN && (b->ZDCgammaN && HFEMaxPlus > ((float)abs(par.DoSystRapGap)) / 10.) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && HFEMaxMinus > ((float)abs(par.DoSystRapGap)) / 10.) == false)
      return false;
  } else if (par.DoSystRapGap == -2) {
    if (par.IsGammaN && (b->ZDCgammaN && HFEMaxPlus >= 9.2) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && HFEMaxMinus >= 8.6) == false)
      return false;
  } else {
    // nominal rapidity gap selection
    if (par.IsGammaN && (b->ZDCgammaN && HFEMaxPlus < 9.2) == false)
      return false;
    if (!par.IsGammaN && (b->ZDCNgamma && HFEMaxMinus < 8.6) == false)
      return false;
  }

  hEventsPassedSelection->Fill(4); // Events after rapidity gap selection

  if (b->nVtx >= 3)
    return false;

  hEventsPassedSelection->Fill(5); // Events after vertex multiplicity selection

  return true;
}

class DataAnalyzer {
public:
  TFile *inf, *outf;
  TH1D *hDmass;
  TH1D *hDmass24003Bins;
  TH1D *hDmass25002Bins;
  DzeroUPCMicroTreeMessenger *MDzeroUPC;
  TNtuple *nt;
  string title;
  TH1D *hDenEvtEff;
  TH1D *hNumEvtEff;
  TH1D *hRatioEvtEff;
  TH1D *hDenDEff;
  TH1D *hNumDEff;
  TH1D *hRatioDEff;

  // Event counting histogram at each event selection step
  TH1D *hEventsPassedSelection;

  // ================================================================ //
  // List of variables used for MC-data comparisons

  // Event level variables:
  std::vector<TH1D *> hHFEmaxPlus;           // gap energy
  std::vector<TH1D *> hHFEmaxMinus;          // gap energy
  std::vector<TH1D *> hnTrackInAcceptanceHP; // event multiplicity
  std::vector<TH1D *> hnVtx, hVX, hVY, hVZ;  // number of verticies and best vertex positions
  std::vector<TH1D *> hDNumPassingDcuts;     // number of D candidates passing selection cuts

  // ================================================================ //

  DataAnalyzer(const char *filename, const char *outFilename, const char *mytitle = "")
      : inf(new TFile(filename)), MDzeroUPC(new DzeroUPCMicroTreeMessenger(*inf, string("Tree"))), title(mytitle),
        outf(new TFile(outFilename, "recreate")) {
    outf->cd();
    nt = new TNtuple("nt", "D0 mass tree", "Dmass:Dgen");
    TTree *tree = MDzeroUPC->Tree;
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

    initHistograms(par);

    WeightHandler2D GptGyWH;
    if (par.DoGptGyReweighting) {
      GptGyWH.LoadFromFile(par.GptGyWeightFileName);
    }

    WeightHandler1D MultWH;
    if (par.DoMultReweighting) {
      MultWH.LoadFromFile(par.MultWeightFileName);
    }

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

      // Check if the event is a signal MC event, a.k.a., having at least one gen-level D candidate that falls into the
      // (pt,y) bin of interest
      bool isSigMCEvt = false;
      double leadingGpt = -999.;
      double leadingGy = -999.;
      if (!par.IsData) {
        for (unsigned long j = 0; j < MDzeroUPC->Gpt->size(); j++) {
          if (MDzeroUPC->GisSignalCalc->at(j) == false)
            continue;
          if (MDzeroUPC->Gpt->at(j) > leadingGpt) {
            leadingGpt = MDzeroUPC->Gpt->at(j);
            leadingGy = MDzeroUPC->Gy->at(j);
          }

          if (MDzeroUPC->Gpt->at(j) < par.MinDzeroPT)
            continue;
          if (MDzeroUPC->Gpt->at(j) > par.MaxDzeroPT)
            continue;
          if (MDzeroUPC->Gy->at(j) < par.MinDzeroY)
            continue;
          if (MDzeroUPC->Gy->at(j) > par.MaxDzeroY)
            continue;
          isSigMCEvt = true;
        }
      }

      double GptGyWeight = 1.0;
      double MultWeight = 1.0;
      if (!par.IsData) {
        if (par.DoGptGyReweighting)
          GptGyWeight = GptGyWH.GetWeight(leadingGpt, leadingGy, (!par.IsGammaN));
        if (par.DoMultReweighting)
          MultWeight = MultWH.GetWeight(MDzeroUPC->nTrackInAcceptanceHP);
      }

      if (!par.IsData && isSigMCEvt)
        hDenEvtEff->Fill(1, GptGyWeight * MultWeight);

      // Check if the event passes the selection criteria
      if (eventSelection(MDzeroUPC, par, hEventsPassedSelection)) {
        bool recoDzeroFound = false;
        bool genDzeroFound = false;
        bool hasSignalD = false;
        double maxDpt = -999.;
        double selectedDmass = -999.;
        int numPassedD = 0;

        if (!par.IsData && isSigMCEvt)
          hNumEvtEff->Fill(1, GptGyWeight * MultWeight);
        bool doTrkFilter = false;
        if (MDzeroUPC->Dtrk1PtErr != nullptr && MDzeroUPC->Dtrk2PtErr != nullptr &&
            MDzeroUPC->Dtrk1PixelHit != nullptr && MDzeroUPC->Dtrk1StripHit != nullptr &&
            MDzeroUPC->Dtrk2PixelHit != nullptr && MDzeroUPC->Dtrk2StripHit != nullptr)
          doTrkFilter = true;

        for (unsigned long j = 0; j < MDzeroUPC->Dsize; j++) {
          if (MDzeroUPC->Dpt->at(j) < par.MinDzeroPT)
            continue;
          if (MDzeroUPC->Dpt->at(j) > par.MaxDzeroPT)
            continue;
          if (MDzeroUPC->Dy->at(j) < par.MinDzeroY)
            continue;
          if (MDzeroUPC->Dy->at(j) > par.MaxDzeroY)
            continue;
          if (par.DoSystD == 0 && MDzeroUPC->DpassCut23PAS->at(j) == false)
            continue;
          // if (par.DoSystD==0 && MDzeroUPC->DpassCut23LowPt->at(j) == false) continue;
          if (par.DoSystD == 1 && MDzeroUPC->DpassCut23PASSystDsvpvSig->at(j) == false)
            continue;
          if (par.DoSystD == 2 && MDzeroUPC->DpassCut23PASSystDtrkPt->at(j) == false)
            continue;
          if (par.DoSystD == 3 && MDzeroUPC->DpassCut23PASSystDalpha->at(j) == false)
            continue;
          if (par.DoSystD == 4 && MDzeroUPC->DpassCut23PASSystDchi2cl->at(j) == false)
            continue;
          //          if (par.DoSystD==0 && MDzeroUPC->DpassCutDefault->at(j) == false) continue;
          //          if (par.DoSystD==1 && MDzeroUPC->DpassCutSystDsvpvSig->at(j) == false) continue;
          //          if (par.DoSystD==2 && MDzeroUPC->DpassCutSystDtrkPt->at(j) == false) continue;
          //          if (par.DoSystD==3 && MDzeroUPC->DpassCutSystDalpha->at(j) == false) continue;
          //          if (par.DoSystD==4 && MDzeroUPC->DpassCutSystDchi2cl->at(j) == false) continue;
          if (doTrkFilter) {
            if ((MDzeroUPC->Dtrk1PtErr->at(j) / MDzeroUPC->Dtrk1Pt->at(j)) > 0.1 ||
                (MDzeroUPC->Dtrk2PtErr->at(j) / MDzeroUPC->Dtrk2Pt->at(j)) > 0.1)
              continue;
            if ((MDzeroUPC->Dtrk1PixelHit->at(j) + MDzeroUPC->Dtrk1StripHit->at(j)) < 11 ||
                (MDzeroUPC->Dtrk2PixelHit->at(j) + MDzeroUPC->Dtrk2StripHit->at(j)) < 11)
              continue;
          }
          recoDzeroFound = true;
          numPassedD++;

          hDmass->Fill((*MDzeroUPC->Dmass)[j]);
          hDmass24003Bins->Fill((*MDzeroUPC->Dmass)[j]);
          hDmass25002Bins->Fill((*MDzeroUPC->Dmass)[j]);
          if (!par.IsData) {
            nt->Fill((*MDzeroUPC->Dmass)[j], (*MDzeroUPC->Dgen)[j]);
            if (MDzeroUPC->Dgen->at(j) == 23333) {
              hNumDEff->Fill(1, GptGyWeight * MultWeight);
            }
          } else
            nt->Fill((*MDzeroUPC->Dmass)[j], 0);

          // Determine if D0 is signal or sideband
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
            hDenDEff->Fill(1, GptGyWeight * MultWeight);
            genDzeroFound = true;
          } // end of gen-level Dzero loop
        } // end of gen-level Dzero loop

        if (par.IsData) {
          fillEvtVarsData(par, recoDzeroFound, hasSignalD, numPassedD);
        } else {
          fillEvtVarsMC(par, recoDzeroFound, genDzeroFound, hasSignalD, numPassedD);
        }

      } // end of event selection
    } // end of event loop
  } // end of analyze

  void writeHistograms(TFile *outf) {
    outf->cd();
    smartWrite(hDmass);
    smartWrite(hDmass24003Bins);
    smartWrite(hDmass25002Bins);
    hRatioEvtEff->Divide(hNumEvtEff, hDenEvtEff, 1, 1, "B");
    hRatioDEff->Divide(hNumDEff, hDenDEff, 1, 1, "B");
    hDenEvtEff->Write();
    hNumEvtEff->Write();
    hRatioEvtEff->Write();
    hDenDEff->Write();
    hNumDEff->Write();
    hRatioDEff->Write();
    hEventsPassedSelection->Write();
    smartWrite(nt);
    for (int i = 0; i < hHFEmaxMinus.size(); i++) {
      smartWrite(hHFEmaxMinus[i]);
      smartWrite(hHFEmaxPlus[i]);
      smartWrite(hnTrackInAcceptanceHP[i]);
      smartWrite(hnVtx[i]);
      smartWrite(hVX[i]);
      smartWrite(hVY[i]);
      smartWrite(hVZ[i]);
      smartWrite(hDNumPassingDcuts[i]);
    }
  }

private:
  void deleteHistograms() {
    delete hDmass;
    delete hDmass24003Bins;
    delete hDmass25002Bins;
    delete hDenEvtEff;
    delete hNumEvtEff;
    delete hRatioEvtEff;
    delete hDenDEff;
    delete hNumDEff;
    delete hRatioDEff;
    delete hEventsPassedSelection;

    for (int i = 0; i < hHFEmaxMinus.size(); i++) {
      delete hHFEmaxMinus[i];
      delete hHFEmaxPlus[i];
      delete hnTrackInAcceptanceHP[i];
      delete hnVtx[i];
      delete hVX[i];
      delete hVY[i];
      delete hVZ[i];
      delete hDNumPassingDcuts[i];
    }
  }

  void initHistograms(Parameters &par) {
    // Initialize histograms here if needed
    hDmass = new TH1D(Form("hDmass%s", title.c_str()), "", DMASS_HISTNBINS, DMASS_HISTMIN, DMASS_HISTMAX);
    hDmass24003Bins = new TH1D(Form("hDmass%s_24003Bins", title.c_str()), "", DMASS_HISTNBINS_24003,
                               DMASS_HISTMIN_24003, DMASS_HISTMAX_24003);
    hDmass25002Bins = new TH1D(Form("hDmass%s_25002Bins", title.c_str()), "", DMASS_HISTNBINS_25002,
                               DMASS_HISTMIN_25002, DMASS_HISTMAX_25002);
    hDenEvtEff = new TH1D(Form("hDenEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumEvtEff = new TH1D(Form("hNumEvtEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioEvtEff = (TH1D *)hNumEvtEff->Clone("hRatioEvtEff");
    hDenDEff = new TH1D(Form("hDenDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hNumDEff = new TH1D(Form("hNumDEff%s", title.c_str()), "", 1, 0.5, 1.5);
    hRatioDEff = (TH1D *)hNumDEff->Clone("hRatioDEff");

    hEventsPassedSelection =
        new TH1D(Form("hEventsPassedSelection%s", title.c_str()), "Events Passed Selection", 6, -0.5, 5.5);
    hEventsPassedSelection->GetXaxis()->SetBinLabel(1, "Total events");
    hEventsPassedSelection->GetXaxis()->SetBinLabel(2, "After trigger");
    hEventsPassedSelection->GetXaxis()->SetBinLabel(3, "After CCF and CSC");
    hEventsPassedSelection->GetXaxis()->SetBinLabel(4, "After vertex selection");
    hEventsPassedSelection->GetXaxis()->SetBinLabel(5, "After rapidity gap cut");
    hEventsPassedSelection->GetXaxis()->SetBinLabel(6, "After vertex multiplicity cut");

    hDmass->Sumw2();
    hDmass25002Bins->Sumw2();
    hDmass24003Bins->Sumw2();
    hDenEvtEff->Sumw2();
    hNumEvtEff->Sumw2();
    hRatioEvtEff->Sumw2();
    hDenDEff->Sumw2();
    hNumDEff->Sumw2();
    hRatioDEff->Sumw2();
    hEventsPassedSelection->Sumw2();

    std::vector<std::string> histTypes = par.IsData ? evtHistTypesData : evtHistTypesMC;
    for (int i = 0; i < histTypes.size(); i++) {
      std::string type = histTypes[i];

      hHFEmaxMinus.push_back(new TH1D(Form("hHFEmaxMinus%s%s", title.c_str(), type.c_str()), "", 400, 0, 100));
      hHFEmaxPlus.push_back(new TH1D(Form("hHFEmaxPlus%s%s", title.c_str(), type.c_str()), "", 400, 0, 100));
      hnTrackInAcceptanceHP.push_back(
          new TH1D(Form("hnTrackInAcceptanceHP%s%s", title.c_str(), type.c_str()), "", 200, 0, 400));
      hnVtx.push_back(new TH1D(Form("hnVtx%s%s", title.c_str(), type.c_str()), "", 10, -0.5, 9.5));
      hVX.push_back(new TH1D(Form("hVX%s%s", title.c_str(), type.c_str()), "", 100, -0.5, 0.5));
      hVY.push_back(new TH1D(Form("hVY%s%s", title.c_str(), type.c_str()), "", 100, -0.5, 0.5));
      hVZ.push_back(new TH1D(Form("hVZ%s%s", title.c_str(), type.c_str()), "", 200, -20, 20));
      hDNumPassingDcuts.push_back(
          new TH1D(Form("hDNumPassingDcuts%s%s", title.c_str(), type.c_str()), "", 20, -0.5, 19.5));

      hHFEmaxMinus[i]->Sumw2();
      hHFEmaxPlus[i]->Sumw2();
      hnTrackInAcceptanceHP[i]->Sumw2();
      hnVtx[i]->Sumw2();
      hVX[i]->Sumw2();
      hVY[i]->Sumw2();
      hVZ[i]->Sumw2();
      hDNumPassingDcuts[i]->Sumw2();
    }
  }

  void fillEvtVarsData(Parameters &par, bool recoDzeroFound, bool hasSignalD, int numPassedD) {
    for (int i = 0; i < hHFEmaxMinus.size(); i++) {
      std::string histName = hHFEmaxMinus[i]->GetName();
      bool doFill = (histName.find("SignalRecoD") != std::string::npos && recoDzeroFound && hasSignalD) ||
                    (histName.find("BkgRecoD") != std::string::npos && recoDzeroFound && !hasSignalD);

      if (doFill) {
        hHFEmaxMinus[i]->Fill(MDzeroUPC->HFEMaxMinus_eta5);
        hHFEmaxPlus[i]->Fill(MDzeroUPC->HFEMaxPlus_eta5);
        hnTrackInAcceptanceHP[i]->Fill(MDzeroUPC->nTrackInAcceptanceHP);
        hnVtx[i]->Fill(MDzeroUPC->nVtx);
        hVX[i]->Fill(MDzeroUPC->VX);
        hVY[i]->Fill(MDzeroUPC->VY);
        hVZ[i]->Fill(MDzeroUPC->VZ);
        hDNumPassingDcuts[i]->Fill(numPassedD);
      }
    }
  }

  void fillEvtVarsMC(Parameters &par, bool recoDzeroFound, bool genDzeroFound, bool hasSignalD, int numPassedD) {
    for (int i = 0; i < hHFEmaxMinus.size(); i++) {
      std::string histName = hHFEmaxMinus[i]->GetName();
      bool doFill = (histName.find("SignalGenD") != std::string::npos && genDzeroFound) ||
                    (histName.find("BkgRecoD") != std::string::npos && recoDzeroFound && !hasSignalD) ||
                    (histName.find("SignalRecoD") != std::string::npos && recoDzeroFound && hasSignalD);

      if (doFill) {
        hHFEmaxMinus[i]->Fill(MDzeroUPC->HFEMaxMinus_eta5);
        hHFEmaxPlus[i]->Fill(MDzeroUPC->HFEMaxPlus_eta5);
        hnTrackInAcceptanceHP[i]->Fill(MDzeroUPC->nTrackInAcceptanceHP);
        hnVtx[i]->Fill(MDzeroUPC->nVtx);
        hVX[i]->Fill(MDzeroUPC->VX);
        hVY[i]->Fill(MDzeroUPC->VY);
        hVZ[i]->Fill(MDzeroUPC->VZ);
        hDNumPassingDcuts[i]->Fill(numPassedD);
      }
    }
  }
};

//============================================================//
// Main analysis
//============================================================//
int main(int argc, char *argv[]) {
  if (printHelpMessage(argc, argv))
    return 0;
  CommandLine CL(argc, argv);
  float MinDzeroPT = CL.GetDouble("MinDzeroPT", 2);  // Minimum Dzero transverse momentum threshold for Dzero selection.
  float MaxDzeroPT = CL.GetDouble("MaxDzeroPT", 5);  // Maximum Dzero transverse momentum threshold for Dzero selection.
  float MinDzeroY = CL.GetDouble("MinDzeroY", -2);   // Minimum Dzero rapidity threshold for Dzero selection.
  float MaxDzeroY = CL.GetDouble("MaxDzeroY", +2);   // Maximum Dzero rapidity threshold for Dzero selection.
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  int TriggerChoice = CL.GetInt("TriggerChoice", 2); // 0 = no trigger sel, 1 = isL1ZDCOr, 2 = isL1ZDCXORJet8
  float scaleFactor = CL.GetDouble("scaleFactor", 1); // Scale factor for the number of events to be processed.
  int DoSystRapGap = CL.GetInt(
      "DoSystRapGap", 0); // Systematic study: apply the alternative event selections
                          // 0 = nominal, 1 = tight, -1: loose
                          // 9 < DoSystRapGap: use custom HF energy threshold, the threshold value will be
                          // DoSystRapGap/10. -2 = nominal rejection mode -9 > DoSystRapGap: use custom HF energy
                          // threshold in rejection mode, the threshold value will be abs(DoSystRapGap)/10.
  int DoSystD = CL.GetInt("DoSystD", 0); // Systematic study: apply the alternative D selections
                                         // 0 = nominal, 1 = Dsvpv variation, 2: DtrkPt variation
                                         // 3 = Dalpha variation, 4: Dchi2cl variation
  int doCCFsyst =
      CL.GetInt("DoCCFsyst", 0); // Systematic study: CCFs variation
                                 // 0 = nominal, 1 = select events with failing ccf, 2 = no ccf condition applied
  bool DoGptGyReweighting = CL.GetBool("DoGptGyReweighting", false);
  string GptGyWeightFileName =
      CL.Get("GptGyWeightFileName", "../../WeightHandler/20250305_DzeroUPC_GptGyWeight/Weights/testWeight.root");
  bool DoMultReweighting = CL.GetBool("DoMultReweighting", false);
  string MultWeightFileName =
      CL.Get("MultWeightFileName", "../../WeightHandler/20250312_DzeroUPC_multiplicityWeight/Weights/testWeight.root");

  bool IsData = CL.GetBool("IsData", 0); // Data or MC
  Parameters par(MinDzeroPT, MaxDzeroPT, MinDzeroY, MaxDzeroY, IsGammaN, TriggerChoice, 1, IsData, 2025, scaleFactor,
                 DoSystRapGap, DoSystD, DoGptGyReweighting, GptGyWeightFileName, DoMultReweighting, MultWeightFileName,
                 doCCFsyst);
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

  return 0;
}
