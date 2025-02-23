#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TFile.h"
#include "TColor.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> // For std::max_element

using namespace std;
#include "CommandLine.h"
#include "parameter.h"

bool checkFileType(string fileType);
bool checkPlottingModes(std::vector<string> modes, bool& needMergedHistograms, std::vector<int> &modeIndexes);

void plotGraph1D(std::vector<TH1D*> inHist, std::vector<string> Modes, std::vector<int> modeIndexes, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy,
                 double xLowLim, const double xHighLim, double yLowLim, const double yHighLim,
                 double legendXLow=0.3, double legendXHigh=0.6, double legendYLow=0.7, double legendYHigh=0.9, double latexY=0.27);

void plotGraph2D(std::vector<TH2D*> inHist, std::vector<string> Modes, std::vector<int> modeIndexes, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy,
                 double xLowLim, const double xHighLim, double yLowLim, const double yHighLim, double latexY=0.27);

void plotOverlapGraph(TH1D* inHist, TH1D* inHist2, string fileType, string fileType2, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy, const double xLowLim, const double xHighLim, string title="",
                 double legendXLow=0.3, double legendXHigh=0.6, double legendYLow=0.7, double legendYHigh=0.9, double latexY=0.77);

void readInHistogramsFromFile(const char* histoSubname, bool needMergedHistograms, TFile* inFile, std::vector<TH2D*> &hHFEmaxPlus_vs_EvtMult, std::vector<TH2D*> &hHFEmaxMinus_vs_EvtMult, std::vector<TH1D*> &hHFEmaxMinus, std::vector<TH1D*> &hHFEmaxPlus, std::vector<TH1D*> &hEvtMult, std::vector<TH1D*> &hnVtx, std::vector<TH1D*> &hVX, std::vector<TH1D*> &hVY, std::vector<TH1D*> &hVZ, std::vector<TH1D*> &hDchi2cl, std::vector<TH1D*> &hDalpha, std::vector<TH1D*> &hDdtheta, std::vector<TH1D*> &hDsvpvDistance, std::vector<TH1D*> &hDsvpvDisErr, std::vector<TH1D*> &hDsvpvSig, std::vector<TH1D*> &hDmass, std::vector<TH2D*> &hDtrk1Pt_vs_Dtrk2Pt, std::vector<TH2D*> &hDpt_vs_Dy, std::vector<TH1D*> &hDtrk1Pt, std::vector<TH1D*> &hDtrk2Pt, std::vector<TH1D*> &hDpt, std::vector<TH1D*> &hDy, std::vector<TH2D*> &htrkPt_vs_trkEta, std::vector<TH1D*> &htrkPt, std::vector<TH1D*> &htrkEta, std::vector<TH2D*> &hGDpt_GDy, std::vector<TH1D*> &hGDpt, std::vector<TH1D*> &hGDy);

int main(int argc, char *argv[])
{
  CommandLine CL(argc, argv);
  string         PlotDir       = CL.Get    ("PlotDir", "");       // subdirectory under plot/
  string         inputDir       = CL.Get    ("InputDir", "");     // parent directory of the file you want to plot
  float MinDzeroPT = CL.GetDouble("MinDzeroPT", 2);  // Minimum Dzero transverse momentum threshold for Dzero selection.
  float MaxDzeroPT = CL.GetDouble("MaxDzeroPT", 5);  // Maximum Dzero transverse momentum threshold for Dzero selection.
  float MinDzeroY = CL.GetDouble("MinDzeroY", 0);  // Minimum Dzero rapidity threshold for Dzero selection.
  float MaxDzeroY = CL.GetDouble("MaxDzeroY", 1);  // Maximum Dzero rapidity threshold for Dzero selection.
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  std::vector<string> Modes = CL.GetStringVector("Modes", "");
  string fileType = CL.Get("fileType", ""); // determine wether the file is Data, MC or MC_inclusive
  string fileType2 = CL.Get("fileType2", ""); // give second fileType only if want to overlap two plots /w ratio determine wether the file is Data, MC or MC_inclusive

  if (!checkFileType(fileType) || (fileType2 != "" && !checkFileType(fileType2) )) { std::cout << "[Error] Invalid fileType entered in config!" << std::endl; return -1; }

  int HFEMax = CL.GetInt("HFEMax",-999); // Read in HFEMax values for plotting

  string inputFileName;
  // Select appropriate folder naming convention depending on if HFEMax was given or not
  if (HFEMax == -999) {
    inputFileName = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, fileType.c_str() );
  } else {
    inputFileName = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o_DoSystRapGap%d/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, (int) HFEMax, fileType.c_str() );
  }

  if (Modes.size() == 0 || Modes[0] == "") { std::cout << "[Error] Plotting mode (LowerSideBand, Signal, UpperSideBand, Merged) not selected! Please modify config file!" << std::endl; return -1; }
  else if (Modes.size() != 1 && fileType2 != "") { std::cout << "[Error] When making ratio plots please only select one plotting mode (LowerSideBand, Signal, UpperSideBand, Merged)!" << std::endl; return -1;}

  bool needMergedHistograms = false;
  std::vector<int> modeIndexes;
  if (!checkPlottingModes(Modes, needMergedHistograms, modeIndexes)) { std::cout << "[Error] Invalid plotting mode selected! Valid options are: LowerSideBand, Signal, UpperSideBand, Merged" << std::endl; return -1; }

  TFile *inFile = new TFile(inputFileName.c_str());

  // Event level variables:
  std::vector<TH2D *> hHFEmaxPlus_vs_EvtMult; // gap energy vs event multiplicity
  std::vector<TH2D *> hHFEmaxMinus_vs_EvtMult; // gap energy vs event multiplicity
  std::vector<TH1D *> hHFEmaxMinus, hHFEmaxPlus, hEvtMult; // projections of the above histograms
  std::vector<TH1D *> hnVtx, hVX, hVY, hVZ; // number of verticies and best vertex positions

  // Track level variables:
  std::vector<TH2D *> htrkPt_vs_trkEta;
  std::vector<TH1D *> htrkPt, htrkEta;

  // D level variables: 

  // reco D level
  std::vector<TH1D *> hDchi2cl, hDalpha, hDdtheta, hDsvpvDistance, hDsvpvDisErr, hDsvpvSig, hDmass;
  std::vector<TH2D *> hDtrk1Pt_vs_Dtrk2Pt, hDpt_vs_Dy;
  std::vector<TH1D *> hDtrk1Pt, hDtrk2Pt, hDpt, hDy;

  // gen D level
  std::vector<TH2D *> hGDpt_GDy;
  std::vector<TH1D *> hGDpt, hGDy;

  readInHistogramsFromFile("_1", needMergedHistograms, inFile, hHFEmaxPlus_vs_EvtMult, hHFEmaxMinus_vs_EvtMult, hHFEmaxMinus, hHFEmaxPlus, hEvtMult, hnVtx, hVX, hVY, hVZ, hDchi2cl, hDalpha, hDdtheta, hDsvpvDistance, hDsvpvDisErr, hDsvpvSig, hDmass, hDtrk1Pt_vs_Dtrk2Pt, hDpt_vs_Dy, hDtrk1Pt, hDtrk2Pt, hDpt, hDy, htrkPt_vs_trkEta, htrkPt, htrkEta, hGDpt_GDy, hGDpt, hGDy);

  /////////////////////////////////
  // 1. Plot the HF energy distributions
  /////////////////////////////////

  const string latexText = Form("%d < D_{p_{T}} < %d (GeV/#it{c})", (int) MinDzeroPT, (int) MaxDzeroPT);
  const string latexText2 = Form("%d < D_{y} < %d", (int) MinDzeroY, (int) MaxDzeroY);

  if (fileType2=="") {
    // Start plotting event level histograms
    for (int iMode = 0; iMode < Modes.size(); iMode++) {
      plotGraph2D(hHFEmaxPlus_vs_EvtMult, {Modes[iMode]}, {modeIndexes[iMode]}, "Charged hadron multiplicity", "HF E_{max} Plus [GeV]",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/HFEmaxPlus_vs_EvtMult/HFEmaxPlus_vs_EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str(), Modes[iMode].c_str()), false, true, 0, 20, 0, 20);
      
      plotGraph2D(hHFEmaxMinus_vs_EvtMult, {Modes[iMode]}, {modeIndexes[iMode]}, "Charged hadron multiplicity", "HF E_{max} Minus [GeV]",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/HFEmaxMinus_vs_EvtMult/HFEmaxMinus_vs_EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str(), Modes[iMode].c_str()), false, true, 0, 20, 0, 20);
    }

    plotGraph1D(hHFEmaxPlus, Modes, modeIndexes, "Normalized counts", "HF E_{max} Plus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxPlus/HFEmaxPlus_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 0, 20, 0, 0, 0.6, 0.8, 0.7, 0.9);

    plotGraph1D(hHFEmaxMinus, Modes, modeIndexes, "Normalized counts", "HF E_{max} Minus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxMinus/HFEmaxMinus_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 0, 20, 0, 0, 0.6, 0.8, 0.7, 0.9);
                  
    plotGraph1D(hEvtMult, Modes, modeIndexes, "Normalized counts", "Charged hadron multiplicity",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/EvtMult/EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 0, 80, 0, 0);

    plotGraph1D(hnVtx, Modes, modeIndexes, "Normalized counts", "Number of vertices",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/nVtx/nVtx_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, -1, 3, 0, 0, 0.15, 0.35);
    
    plotGraph1D(hVX, Modes, modeIndexes, "Normalized counts", "VX [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VX/VX_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, -1, 1, 0, 0, 0.15, 0.35);

    plotGraph1D(hVY, Modes, modeIndexes, "Normalized counts", "VY [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VY/VY_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, -1, 1, 0, 0, 0.15, 0.35);

    plotGraph1D(hVZ, Modes, modeIndexes, "Normalized counts", "VZ [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VZ/VZ_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, -15, 15, 0, 0, 0.15, 0.35);

    // Start plotting track level histograms

    if (htrkPt_vs_trkEta[0] != nullptr) {
      plotGraph2D(htrkPt_vs_trkEta, {"Merged"}, {0}, "track p_{T} [GeV]", "track #eta",
              latexText.c_str(), latexText2.c_str(),
                  Form("%s/trkPt_vs_trkEta/trkPt_vs_trkEta_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), false, false, 0, 15, -2.4, 2.4);

      plotGraph1D(htrkPt, {"Merged"}, {0}, "Normalized counts", "track p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
                  Form("%s/trkPt/trkPt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), false, true, 0, 15, 0, 0);

      plotGraph1D(htrkEta, {"Merged"}, {0}, "Normalized counts", "track #eta",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/trkEta/trkEta_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), false, false, -2.4, 2.4, 0, 0, 0.3, 0.6, 0.2, 0.4, 0.2);
    }

    // Start plotting D level histograms

    plotGraph1D(hDalpha, Modes, modeIndexes, "Normalized counts", "D^{0} #alpha",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dalpha/Dalpha_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 0, 0.5, 0, 0, 0.2, 0.4, 0.2, 0.4);
    
    plotGraph1D(hDchi2cl, Modes, modeIndexes, "Normalized counts", "D^{0} #chi^{2}",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dchi2cl/Dchi2cl_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 0, 1., 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph1D(hDdtheta, Modes, modeIndexes, "Normalized counts", "D^{0} #Delta#theta",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Ddtheta/Ddtheta_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 0, 0.6, 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph1D(hDsvpvDisErr, Modes, modeIndexes, "Normalized counts", "DsvpvDisErr",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDisErr/DsvpvDisErr_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 0, 1.2, 0, 0);

    plotGraph1D(hDsvpvDistance, Modes, modeIndexes, "Normalized counts", "DsvpvDistance",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDistance/DsvpvDistance_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 0, 5, 0, 0);

    plotGraph1D(hDsvpvSig, Modes, modeIndexes, "Normalized counts", "DsvpvSig",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/DsvpvSig/DsvpvSig_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), false, true, 0, 15, 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph1D(hDmass, Modes, modeIndexes, "Normalized counts", "D^{0} mass [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dmass/Dmass_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 1.7, 2., 0, 0);

    for (int iMode = 0; iMode < Modes.size(); iMode++) {
      plotGraph2D(hDtrk1Pt_vs_Dtrk2Pt, {Modes[iMode]}, {modeIndexes[iMode]}, "D^{0} trk2 p_{T} [GeV]", "D^{0} trk1 p_{T} [GeV]",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/Dtrk1Pt_vs_Dtrk2Pt/Dtrk1Pt_vs_Dtrk2Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str(), Modes[iMode].c_str()), false, false, 0, 5, 0, 5);
      
      plotGraph2D(hDpt_vs_Dy, {Modes[iMode]}, {modeIndexes[iMode]}, "D^{0} rapidity", "D^{0} p_{T} [GeV]",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/Dpt_vs_Dy/Dpt_vs_Dy_pt%d-%d_y%d-%d_IsGammaN%o_%s_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str(), Modes[iMode].c_str()), false, false, 200, -200, 200, -200);
    }

    plotGraph1D(hDtrk1Pt, Modes, modeIndexes, "Normalized counts", "D^{0} trk1 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk1Pt/Dtrk1Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 0, 5, 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph1D(hDtrk2Pt, Modes, modeIndexes, "Normalized counts", "D^{0} trk2 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk2Pt/Dtrk2Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 0, 5, 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph1D(hDpt, Modes, modeIndexes, "Normalized counts", "D^{0} p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dpt/Dpt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 200, -200, 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph1D(hDy, Modes, modeIndexes, "Normalized counts", "D^{0} y",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dy/Dy_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 200, -200, 0, 0, 0.2, 0.4, 0.2, 0.4);

    plotGraph2D(hGDpt_GDy, {"Merged"}, {0}, "Gen D^{0} p_{T} [GeV]", "Gen D^{0} rapidity",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/GDpt_GDy/GDpt_GDy_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), false, false, 200, -200, 200, -200);
    plotGraph1D(hGDpt, {"Merged"}, {0}, "Normalized counts", "Gen D^{0} pt [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/GDpt/GDpt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, true, 200, -200, 0, 0, 0.2, 0.4, 0.2, 0.4);
    plotGraph1D(hGDy, {"Merged"}, {0}, "Normalized counts", "Gen D^{0} y",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/GDy/GDy_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), false, false, 200, -200, 0, 0, 0.2, 0.4, 0.2, 0.4);

  } else {
    string inputFileName2;
    if (HFEMax == -999) {
      inputFileName2 = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, fileType2.c_str() );
    } else {
      inputFileName2 = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o_DoSystRapGap%d/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, (int) HFEMax, fileType2.c_str() );
    }
    TFile *inFile2 = new TFile(inputFileName2.c_str());
    
    // Event level variables:
    std::vector<TH2D *> hHFEmaxPlus_vs_EvtMult_2; // gap energy vs event multiplicity
    std::vector<TH2D *> hHFEmaxMinus_vs_EvtMult_2; // gap energy vs event multiplicity
    std::vector<TH1D *> hHFEmaxMinus_2, hHFEmaxPlus_2, hEvtMult_2; // projections of the above histograms
    std::vector<TH1D *> hnVtx_2, hVX_2, hVY_2, hVZ_2; // number of verticies and best vertex positions

    // Track level variables:
    std::vector<TH2D *> htrkPt_vs_trkEta_2;
    std::vector<TH1D *> htrkPt_2, htrkEta_2;

    // D level variables: 

    std::vector<TH1D *> hDchi2cl_2, hDalpha_2, hDdtheta_2, hDsvpvDistance_2, hDsvpvDisErr_2, hDsvpvSig_2, hDmass_2;
    std::vector<TH2D *> hDtrk1Pt_vs_Dtrk2Pt_2, hDpt_vs_Dy_2;
    std::vector<TH1D *> hDtrk1Pt_2, hDtrk2Pt_2, hDpt_2, hDy_2;

    // gen D level
    std::vector<TH2D *> hGDpt_GDy_2;
    std::vector<TH1D *> hGDpt_2, hGDy_2;

    readInHistogramsFromFile("_2", needMergedHistograms, inFile2, hHFEmaxPlus_vs_EvtMult_2, hHFEmaxMinus_vs_EvtMult_2, hHFEmaxMinus_2, hHFEmaxPlus_2, hEvtMult_2, hnVtx_2, hVX_2, hVY_2, hVZ_2, hDchi2cl_2, hDalpha_2, hDdtheta_2, hDsvpvDistance_2, hDsvpvDisErr_2, hDsvpvSig_2, hDmass_2, hDtrk1Pt_vs_Dtrk2Pt_2, hDpt_vs_Dy_2, hDtrk1Pt_2, hDtrk2Pt_2, hDpt_2, hDy_2, htrkPt_vs_trkEta_2, htrkPt_2, htrkEta_2, hGDpt_GDy_2, hGDpt_2, hGDy_2);

    // Start plotting event level histograms

    plotOverlapGraph(hHFEmaxPlus[modeIndexes[0]], hHFEmaxPlus_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "HF E_{max} Plus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxPlus/HFEmaxPlus_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 20, Modes[0]);

    plotOverlapGraph(hHFEmaxMinus[modeIndexes[0]], hHFEmaxMinus_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "HF E_{max} Minus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxMinus/HFEmaxMinus_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 20, Modes[0]);
                  
    plotOverlapGraph(hEvtMult[modeIndexes[0]], hEvtMult_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "Charged hadron multiplicity",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/EvtMult/EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 0, 80, Modes[0]);

    plotOverlapGraph(hnVtx[modeIndexes[0]], hnVtx_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "Number of vertices",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/nVtx/nVtx_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -1, 3, Modes[0], 0.15, 0.35);
    
    plotOverlapGraph(hVX[modeIndexes[0]], hVX_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "VX [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VX/VX_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -1, 1, Modes[0]);

    plotOverlapGraph(hVY[modeIndexes[0]], hVY_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "VY [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VY/VY_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -1, 1, Modes[0]);

    plotOverlapGraph(hVZ[modeIndexes[0]], hVZ_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "VZ [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VZ/VZ_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -15, 15, Modes[0], 0.15, 0.35);

    // Start plotting track level histograms

    if (htrkPt_vs_trkEta[0] != nullptr && htrkPt_vs_trkEta_2[0] != nullptr) {

      plotOverlapGraph(htrkPt[0], htrkPt_2[0], fileType, fileType2, "Normalized counts", "track p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/trkPt/trkPt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 15, "Merged");

      plotOverlapGraph(htrkEta[0], htrkEta_2[0], fileType, fileType2, "Normalized counts", "track #eta",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/trkEta/trkEta_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, -2.4, 2.4, "Merged", 0.3, 0.6, 0.2, 0.4, 0.2);
                    
    }

    // Start plotting D level histograms

    plotOverlapGraph(hDalpha[modeIndexes[0]], hDalpha_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} #alpha",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dalpha/Dalpha_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 0, 0.5, Modes[0], 0.6, 0.8, 0.2, 0.4);
    
    plotOverlapGraph(hDchi2cl[modeIndexes[0]], hDchi2cl_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} #chi^{2}",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dchi2cl/Dchi2cl_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 0, 1, Modes[0], 0.3, 0.6, 0.2, 0.4, 0.1);

    plotOverlapGraph(hDdtheta[modeIndexes[0]], hDdtheta_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} #Delta#theta",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Ddtheta/Ddtheta_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 0, 0.6, Modes[0], 0.2, 0.4, 0.2, 0.4);

    plotOverlapGraph(hDsvpvDisErr[modeIndexes[0]], hDsvpvDisErr_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "DsvpvDisErr",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDisErr/DsvpvDisErr_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 1.2, Modes[0]);

    plotOverlapGraph(hDsvpvDistance[modeIndexes[0]], hDsvpvDistance_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "DsvpvDistance",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDistance/DsvpvDistance_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 5, Modes[0]);

    plotOverlapGraph(hDsvpvSig[modeIndexes[0]], hDsvpvSig_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "DsvpvSig",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/DsvpvSig/DsvpvSig_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 15, Modes[0], 0.35, 0.6);

    plotOverlapGraph(hDmass[modeIndexes[0]], hDmass_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} mass [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dmass/Dmass_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 1.7, 2, Modes[0]);

    plotOverlapGraph(hDtrk1Pt[modeIndexes[0]], hDtrk1Pt_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} trk1 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk1Pt/Dtrk1Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 5, Modes[0], 0.3, 0.6, 0.2, 0.4);

    plotOverlapGraph(hDtrk2Pt[modeIndexes[0]], hDtrk2Pt_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} trk2 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk2Pt/Dtrk2Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 5, Modes[0], 0.3, 0.6, 0.2, 0.4);

    plotOverlapGraph(hDpt[modeIndexes[0]], hDpt_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dpt/Dpt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 200, -200, Modes[0]);

    plotOverlapGraph(hDy[modeIndexes[0]], hDy_2[modeIndexes[0]], fileType, fileType2, "Normalized counts", "D^{0} y [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dy/Dy_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 200, -200, Modes[0]);
  }

  return 0;
}

void plotGraph1D(std::vector<TH1D*> inHist, std::vector<string> Modes, std::vector<int> modeIndexes, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy, 
                 double xLowLim, const double xHighLim, double yLowLim, const double yHighLim, 
                 double legendXLow, double legendXHigh, double legendYLow, double legendYHigh, double latexY)
{
  std::vector<int> colors = {kBlue, kGreen, kRed, 1};

  // Create canvas
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.08);
  double minVal = 0;
  if (logx) { gPad->SetLogx(); if (!(xLowLim > 0)) { xLowLim = 1e-5; }}
  if (logy) { gPad->SetLogy(); minVal = 1e-5; }
  
  int overallFirstNonEmptyBin = -999;
  int overallLastNonEmptyBin = -999;
  double overallMax = -999;

  // Create and configure the histogram frame
  for (int iMode = 0; iMode < Modes.size(); iMode++) {
      int iHist = modeIndexes[iMode];
      inHist[iHist]->Scale(1./inHist[iHist]->Integral());
      inHist[iHist]->GetYaxis()->SetTitle(yAxisTitle);
      inHist[iHist]->GetXaxis()->SetTitle(xAxisTitle);

      int firstNonEmptyBin = inHist[iHist]->FindFirstBinAbove(0,1);
      int lastNonEmptyBin = inHist[iHist]->FindLastBinAbove(0,1);
      double maxVal = inHist[iHist]->GetMaximum();
      if (overallFirstNonEmptyBin < 0 || firstNonEmptyBin < overallFirstNonEmptyBin) { overallFirstNonEmptyBin = firstNonEmptyBin; }
      if (overallLastNonEmptyBin < 0 || lastNonEmptyBin > overallLastNonEmptyBin) { overallLastNonEmptyBin = lastNonEmptyBin; }
      if (maxVal > overallMax) { overallMax = maxVal; }

      inHist[iHist]->SetStats(0);
      inHist[iHist]->SetTitle("");
      inHist[iHist]->GetYaxis()->SetTitleOffset(1.5);
  }

  for (int iMode = 0; iMode < Modes.size(); iMode++) {
      int iHist = modeIndexes[iMode];
      inHist[iHist]->GetYaxis()->SetRangeUser(minVal,overallMax*1.2);
      inHist[iHist]->GetXaxis()->SetRangeUser(TMath::Min(xLowLim, inHist[iHist]->GetXaxis()->GetBinLowEdge(overallFirstNonEmptyBin)), TMath::Max(xHighLim, inHist[iHist]->GetXaxis()->GetBinUpEdge(overallLastNonEmptyBin)));
      
      inHist[iHist]->SetLineColor(colors[iMode]);
      if (iMode == 0) { inHist[iHist]->Draw(); }
      else { inHist[iHist]->Draw("same"); }
  }

  // Create legend
  TLegend* legend = new TLegend(legendXLow, legendYLow, legendXHigh, legendYHigh);
  for (int iMode = 0; iMode < Modes.size(); iMode++) { legend->AddEntry(inHist[modeIndexes[iMode]], Modes[iMode].c_str(), "l"); }
  legend->SetFillStyle(0);  // Transparent fill
  legend->SetLineStyle(0);  // Transparent border
  legend->SetBorderSize(0); // No border
  legend->SetTextSize(0.035);
  legend->Draw();

  // Add TLatex for additional text
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.SetTextFont(42);
  latex.DrawLatex(0.6, latexY+0.05, latexText);
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.035);
  latex2.SetTextFont(42);
  latex2.DrawLatex(0.6, latexY, latexText2);

  // Update and save the canvas
  canvas->Update();
  canvas->SaveAs(plotname);

  // Clean up
  delete legend;
  delete canvas;
}

void plotGraph2D(std::vector<TH2D*> inHist, std::vector<string> Modes, std::vector<int> modeIndexes, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy, 
                 double xLowLim, const double xHighLim, double yLowLim, const double yHighLim, double latexY)
{
  if (Modes.size() > 1) { std::cout << "[Warning] Multiple 2D histograms plotted on top of eachother!" << std::endl; }

  // Create canvas
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.12);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.08);
  if (logx) { gPad->SetLogx(); if (!(xLowLim > 0)) { xLowLim = 1e-5; }}
  if (logy) { gPad->SetLogy(); if (!(yLowLim > 0)) { yLowLim = 1e-5; }}
  gPad->SetLogz();

  int overallFirstNonEmptyBin = -999;
  int overallLastNonEmptyBin = -999;

  // Create and configure the histogram frame
  for (int iMode = 0; iMode < Modes.size(); iMode++) {
      int iHist = modeIndexes[iMode];
      inHist[iHist]->Scale(1./inHist[iHist]->Integral());
      inHist[iHist]->GetYaxis()->SetTitle(yAxisTitle);
      inHist[iHist]->GetXaxis()->SetTitle(xAxisTitle);

      int firstNonEmptyBin = TMath::Max(inHist[iHist]->GetXaxis()->GetFirst(),(inHist[iHist]->FindFirstBinAbove(0,1)));
      int lastNonEmptyBin = TMath::Min(inHist[iHist]->GetXaxis()->GetLast(),inHist[iHist]->FindLastBinAbove(0,1));
      double maxVal = inHist[iHist]->GetMaximum();
      if (overallFirstNonEmptyBin < 0 || firstNonEmptyBin < overallFirstNonEmptyBin) { overallFirstNonEmptyBin = firstNonEmptyBin; }
      if (overallLastNonEmptyBin < 0 || lastNonEmptyBin > overallLastNonEmptyBin) { overallLastNonEmptyBin = lastNonEmptyBin; }

      int firstNonEmptyBiny = TMath::Max(inHist[iHist]->GetYaxis()->GetFirst(),inHist[iHist]->FindFirstBinAbove(0,2));
      int lastNonEmptyBiny = TMath::Min(inHist[iHist]->GetYaxis()->GetLast(),inHist[iHist]->FindLastBinAbove(0,2));
      inHist[iHist]->GetYaxis()->SetRangeUser(TMath::Min(yLowLim, inHist[iHist]->GetYaxis()->GetBinLowEdge(firstNonEmptyBiny)), TMath::Max(yHighLim, inHist[iHist]->GetYaxis()->GetBinUpEdge(lastNonEmptyBiny))); 

      inHist[iHist]->SetStats(0);
      inHist[iHist]->SetTitle("");
      inHist[iHist]->GetYaxis()->SetTitleOffset(1.5);
  }

  for (int iMode = 0; iMode < Modes.size(); iMode++) {
      int iHist = modeIndexes[iMode];
      inHist[iHist]->GetXaxis()->SetRangeUser(TMath::Min(xLowLim, inHist[iHist]->GetXaxis()->GetBinLowEdge(overallFirstNonEmptyBin)), TMath::Max(xHighLim, inHist[iHist]->GetXaxis()->GetBinUpEdge(overallLastNonEmptyBin)));
      inHist[iHist]->Draw("colz");
  }

  // Add TLatex for additional text
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.SetTextFont(42);
  latex.DrawLatex(0.6, latexY+0.05, latexText);
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.035);
  latex2.SetTextFont(42);
  latex2.DrawLatex(0.6, latexY, latexText2);

  // Update and save the canvas
  // FIXME: this Update() call sometimes freezes for almost empty histograms
  canvas->Update();
  canvas->SaveAs(plotname);

  // Clean up
  delete canvas;
}

void plotOverlapGraph(TH1D* inHist, TH1D* inHist2, string fileType, string fileType2, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy, const double xLowLim, const double xHighLim, string title,
                 double legendXLow, double legendXHigh, double legendYLow, double legendYHigh, double latexY)
{
  // Create canvas
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.08);

  TPad* padUp = new TPad("padUp", "padUp", 0, 0.3, 1, 1);
  padUp->SetBottomMargin(0);
  padUp->Draw();
  padUp->cd();
  if (logx) { gPad->SetLogx(); }
  if (logy) { gPad->SetLogy(); }
  

  auto styleUpHist = [&title](TH1 *hist) {
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->ChangeLabel(1,0,0);
    hist->SetTitle(title.c_str());
  };

  auto styleDownHist = [](TH1 *hist) {
    hist->GetXaxis()->SetTitleSize(0.05 / 0.4);
    hist->GetXaxis()->SetLabelSize(0.045 / 0.4);
    hist->GetYaxis()->SetTitleSize(0.05 / 0.4);
    hist->GetYaxis()->SetLabelSize(0.045 / 0.4);
    hist->GetXaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleOffset(.3);
    hist->GetYaxis()->ChangeLabel(-1,0,0);
    hist->SetTitle("");
    hist->GetYaxis()->SetNdivisions(505, true);
  };

  // Configure the histogram frame
  inHist->Scale(1./inHist->Integral());
  inHist2->Scale(1./inHist2->Integral());
  int firstNonEmptyBin = TMath::Min(inHist->FindFirstBinAbove(0,1),inHist2->FindFirstBinAbove(0,1));
  int lastNonEmptyBin = TMath::Max(inHist->FindLastBinAbove(0,1),inHist2->FindLastBinAbove(0,1));
  inHist->GetXaxis()->SetRangeUser(TMath::Min(xLowLim, inHist->GetXaxis()->GetBinLowEdge(firstNonEmptyBin)), TMath::Max(xHighLim, inHist->GetXaxis()->GetBinUpEdge(lastNonEmptyBin)));

  inHist->GetYaxis()->SetTitle(yAxisTitle);
  inHist->GetXaxis()->SetTitle(xAxisTitle);
  inHist->SetStats(0);
  inHist->SetLineColor(kBlue);
  inHist->GetYaxis()->SetTitleOffset(1.);
  inHist->GetYaxis()->SetRangeUser(1e-4,1.2*TMath::Max(inHist->GetMaximum(),inHist2->GetMaximum()));
  styleUpHist(inHist);
  inHist->Draw("hist f e");

  inHist2->GetYaxis()->SetTitle(yAxisTitle);
  inHist2->GetXaxis()->SetTitle(xAxisTitle);
  inHist2->SetStats(0);
  inHist2->SetLineColor(kRed);
  inHist2->GetYaxis()->SetTitleOffset(1.);
  inHist2->Draw("same");

  TLegend* legend = new TLegend(legendXLow, legendYLow, legendXHigh, legendYHigh);
  legend->AddEntry(inHist, fileType.c_str(), "l");
  legend->AddEntry(inHist2, fileType2.c_str(), "l");
  legend->SetFillStyle(0);  // Transparent fill
  legend->SetLineStyle(0);  // Transparent border
  legend->SetBorderSize(0); // No border
  legend->SetTextSize(0.045);
  legend->Draw();

  // Add TLatex for additional text
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.045);
  latex.SetTextFont(42);
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.045);
  latex2.SetTextFont(42);
  latex.DrawLatex(0.63, latexY, latexText);
  latex2.DrawLatex(0.63, latexY+0.05, latexText2);

  canvas->cd();
  TPad* padDown = new TPad("padDown", "padDown", 0, 0.05, 1, .3);
  padDown->SetTopMargin(0);
  padDown->SetBottomMargin(.4);
  padDown->Draw();
  padDown->cd();

  TH1D* ratioHist = (TH1D*) inHist->Clone("ratioHist");
  ratioHist->SetStats(0);
  //ratioHist->Sumw2();
  ratioHist->GetYaxis()->SetTitle(Form("%s/%s",fileType.c_str(),fileType2.c_str()));
  ratioHist->Divide(inHist2);
  ratioHist->GetYaxis()->SetRangeUser(0, 2);
  styleDownHist(ratioHist);
  ratioHist->Draw("ep");

  // Update and save the canvas
  canvas->Update();
  canvas->SaveAs(plotname);

  // Clean up
  delete padUp;
  delete padDown;
  delete legend;
  delete canvas;
}

bool checkFileType(string fileType) {
  if (fileType == "Data" || fileType == "MC" || fileType == "MC_inclusive") {
    return true;
  }
  return false;
}

void readInHistogramsFromFile(const char* histoSubname, bool needMergedHistograms, TFile* inFile, std::vector<TH2D*> &hHFEmaxPlus_vs_EvtMult, std::vector<TH2D*> &hHFEmaxMinus_vs_EvtMult, std::vector<TH1D*> &hHFEmaxMinus, std::vector<TH1D*> &hHFEmaxPlus, std::vector<TH1D*> &hEvtMult, std::vector<TH1D*> &hnVtx, std::vector<TH1D*> &hVX, std::vector<TH1D*> &hVY, std::vector<TH1D*> &hVZ, std::vector<TH1D*> &hDchi2cl, std::vector<TH1D*> &hDalpha, std::vector<TH1D*> &hDdtheta, std::vector<TH1D*> &hDsvpvDistance, std::vector<TH1D*> &hDsvpvDisErr, std::vector<TH1D*> &hDsvpvSig, std::vector<TH1D*> &hDmass, std::vector<TH2D*> &hDtrk1Pt_vs_Dtrk2Pt, std::vector<TH2D*> &hDpt_vs_Dy, std::vector<TH1D*> &hDtrk1Pt, std::vector<TH1D*> &hDtrk2Pt, std::vector<TH1D*> &hDpt, std::vector<TH1D*> &hDy, std::vector<TH2D*> &htrkPt_vs_trkEta, std::vector<TH1D*> &htrkPt, std::vector<TH1D*> &htrkEta, std::vector<TH2D*> &hGDpt_GDy, std::vector<TH1D*> &hGDpt, std::vector<TH1D*> &hGDy) {
  htrkPt_vs_trkEta.push_back( (TH2D*) inFile->Get("htrkPt_vs_trkEta") );
  htrkPt.push_back(nullptr);
  htrkEta.push_back(nullptr);
  if (htrkPt_vs_trkEta[0] != nullptr) {
      htrkPt_vs_trkEta[0]->SetName(Form("htrkPt_vs_trkEta%s",histoSubname));
      htrkPt[0] = (TH1D*) htrkPt_vs_trkEta[0]->ProjectionX(Form("htrkPt%s",histoSubname));
      htrkEta[0] = (TH1D*) htrkPt_vs_trkEta[0]->ProjectionY(Form("htrkEta%s",histoSubname));
  }

  hGDpt_GDy.push_back((TH2D*) inFile->Get("hGDpt_GDy"));
  hGDpt_GDy[0]->SetName(Form("hGDpt_GDy%s", histoSubname));
  hGDpt.push_back((TH1D*) hGDpt_GDy[0]->ProjectionX(Form("hGDpt%s",histoSubname)));
  hGDy.push_back((TH1D*) hGDpt_GDy[0]->ProjectionY(Form("hGDy%s",histoSubname)));

  for (int iHistName = 0; iHistName < Parameters::HISTO_NAMES.size(); iHistName++) {
      hHFEmaxPlus_vs_EvtMult.push_back( (TH2D*) inFile->Get(Form("hHFEmaxPlus_vs_EvtMult_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hHFEmaxPlus_vs_EvtMult[iHistName]->SetName(Form("hHFEmaxPlus_vs_EvtMult_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hHFEmaxPlus.push_back( (TH1D*) hHFEmaxPlus_vs_EvtMult[iHistName]->ProjectionX(Form("hHFEmaxPlus_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );
      hHFEmaxPlus[iHistName]->Rebin(2);

      hHFEmaxMinus_vs_EvtMult.push_back( (TH2D*) inFile->Get(Form("hHFEmaxMinus_vs_EvtMult_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hHFEmaxMinus_vs_EvtMult[iHistName]->SetName(Form("hHFEmaxMinus_vs_EvtMult_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hHFEmaxMinus.push_back( (TH1D*) hHFEmaxMinus_vs_EvtMult[iHistName]->ProjectionX(Form("hHFEmaxMinus_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );
      hHFEmaxMinus[iHistName]->Rebin(2);

      hEvtMult.push_back( (TH1D*) hHFEmaxPlus_vs_EvtMult[iHistName]->ProjectionY(Form("hEvtMult_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );

      hnVtx.push_back( (TH1D*) inFile->Get(Form("hnVtx_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hnVtx[iHistName]->SetName(Form("hnVtx_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hVX.push_back( (TH1D*) inFile->Get(Form("hVX_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hVX[iHistName]->SetName(Form("hVX_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hVY.push_back( (TH1D*) inFile->Get(Form("hVY_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hVY[iHistName]->SetName(Form("hVY_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hVZ.push_back( (TH1D*) inFile->Get(Form("hVZ_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hVZ[iHistName]->SetName(Form("hVZ_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hVZ[iHistName]->Rebin(2);

      hDchi2cl.push_back( (TH1D*) inFile->Get(Form("hDchi2cl_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDchi2cl[iHistName]->SetName(Form("hDchi2cl_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDchi2cl[iHistName]->Rebin(2);
      hDalpha.push_back( (TH1D*) inFile->Get(Form("hDalpha_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDalpha[iHistName]->SetName(Form("hDalpha_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDdtheta.push_back( (TH1D*) inFile->Get(Form("hDdtheta_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDdtheta[iHistName]->SetName(Form("hDdtheta_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDmass.push_back( (TH1D*) inFile->Get(Form("hDmass_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDmass[iHistName]->SetName(Form("hDmass_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDmass[iHistName]->Rebin(2);
      hDsvpvSig.push_back( (TH1D*) inFile->Get(Form("hDsvpvSig_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDsvpvSig[iHistName]->SetName(Form("hDsvpvSig_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDsvpvSig[iHistName]->Rebin(2);
      hDsvpvDistance.push_back( (TH1D*) inFile->Get(Form("hDsvpvDistance_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDsvpvDistance[iHistName]->SetName(Form("hDsvpvDistance_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDsvpvDistance[iHistName]->Rebin(2);
      hDsvpvDisErr.push_back( (TH1D*) inFile->Get(Form("hDsvpvDisErr_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDsvpvDisErr[iHistName]->SetName(Form("hDsvpvDisErr_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDsvpvDisErr[iHistName]->Rebin(2);

      hDtrk1Pt_vs_Dtrk2Pt.push_back( (TH2D*) inFile->Get(Form("hDtrk1Pt_vs_Dtrk2Pt_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDtrk1Pt_vs_Dtrk2Pt[iHistName]->SetName(Form("hDtrk1Pt_vs_Dtrk2Pt_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));
      hDpt_vs_Dy.push_back( (TH2D*) inFile->Get(Form("hDpt_vs_Dy_%s", Parameters::HISTO_NAMES[iHistName].c_str())) );
      hDpt_vs_Dy[iHistName]->SetName(Form("hDpt_vs_Dy_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname));

      hDtrk1Pt.push_back( (TH1D*) hDtrk1Pt_vs_Dtrk2Pt[iHistName]->ProjectionX(Form("hDtrk1Pt_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );
      hDtrk2Pt.push_back( (TH1D*) hDtrk1Pt_vs_Dtrk2Pt[iHistName]->ProjectionY(Form("hDtrk2Pt_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );
      hDpt.push_back( (TH1D*) hDpt_vs_Dy[iHistName]->ProjectionX(Form("hDpt_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );
      hDpt[iHistName]->Rebin(5);
      hDy.push_back( (TH1D*) hDpt_vs_Dy[iHistName]->ProjectionY(Form("hDy_%s%s", Parameters::HISTO_NAMES[iHistName].c_str(), histoSubname)) );
      hDy[iHistName]->Rebin(2);
  }

  // compute merged histograms if required
  if (needMergedHistograms) {
      hHFEmaxPlus_vs_EvtMult.push_back((TH2D*) hHFEmaxPlus_vs_EvtMult[0]->Clone(Form("hHFEmaxPlus_vs_EvtMult_Merged%s", histoSubname)));
      hHFEmaxPlus.push_back((TH1D*) hHFEmaxPlus[0]->Clone(Form("hHFEmaxPlus_Merged%s", histoSubname)));

      hHFEmaxMinus_vs_EvtMult.push_back((TH2D*) hHFEmaxMinus_vs_EvtMult[0]->Clone(Form("hHFEmaxMinus_vs_EvtMult_Merged%s", histoSubname)));
      hHFEmaxMinus.push_back((TH1D*) hHFEmaxMinus[0]->Clone(Form("hHFEmaxMinus_Merged%s", histoSubname)));

      hEvtMult.push_back((TH1D*) hEvtMult[0]->Clone(Form("hEvtMult_Merged%s", histoSubname)));

      hnVtx.push_back((TH1D*) hnVtx[0]->Clone(Form("hnVtx_Merged%s", histoSubname)));
      hVX.push_back((TH1D*) hVX[0]->Clone(Form("hVX_Merged%s", histoSubname)));
      hVY.push_back((TH1D*) hVY[0]->Clone(Form("hVY_Merged%s", histoSubname)));
      hVZ.push_back((TH1D*) hVZ[0]->Clone(Form("hVZ_Merged%s", histoSubname)));

      hDchi2cl.push_back((TH1D*) hDchi2cl[0]->Clone(Form("hDchi2cl_Merged%s", histoSubname)));
      hDalpha.push_back((TH1D*) hDalpha[0]->Clone(Form("hDalpha_Merged%s", histoSubname)));
      hDdtheta.push_back((TH1D*) hDdtheta[0]->Clone(Form("hDdtheta_Merged%s", histoSubname)));
      hDmass.push_back((TH1D*) hDmass[0]->Clone(Form("hDmass_Merged%s", histoSubname)));
      hDsvpvSig.push_back((TH1D*) hDsvpvSig[0]->Clone(Form("hDsvpvSig_Merged%s", histoSubname)));
      hDsvpvDistance.push_back((TH1D*) hDsvpvDistance[0]->Clone(Form("hDsvpvDistance_Merged%s", histoSubname)));
      hDsvpvDisErr.push_back((TH1D*) hDsvpvDisErr[0]->Clone(Form("hDsvpvDisErr_Merged%s", histoSubname)));

      hDtrk1Pt_vs_Dtrk2Pt.push_back((TH2D*) hDtrk1Pt_vs_Dtrk2Pt[0]->Clone(Form("hDtrk1Pt_vs_Dtrk2Pt_Merged%s", histoSubname)));
      hDpt_vs_Dy.push_back((TH2D*) hDpt_vs_Dy[0]->Clone(Form("hDpt_vs_Dy_Merged%s", histoSubname)));

      hDtrk1Pt.push_back((TH1D*) hDtrk1Pt[0]->Clone(Form("hDtrk1Pt_Merged%s", histoSubname)));
      hDtrk2Pt.push_back((TH1D*) hDtrk2Pt[0]->Clone(Form("hDtrk2Pt_Merged%s", histoSubname)));
      hDpt.push_back((TH1D*) hDpt[0]->Clone(Form("hDpt_Merged%s", histoSubname)));
      hDy.push_back((TH1D*) hDy[0]->Clone(Form("hDy_Merged%s", histoSubname)));

      for (int iHistName = 1; iHistName < Parameters::HISTO_NAMES.size(); iHistName++) {
            hHFEmaxPlus_vs_EvtMult[Parameters::HISTO_NAMES.size()]->Add(hHFEmaxPlus_vs_EvtMult[iHistName]);
            hHFEmaxPlus[Parameters::HISTO_NAMES.size()]->Add(hHFEmaxPlus[iHistName]);

            hHFEmaxMinus_vs_EvtMult[Parameters::HISTO_NAMES.size()]->Add(hHFEmaxMinus_vs_EvtMult[iHistName]);
            hHFEmaxMinus[Parameters::HISTO_NAMES.size()]->Add(hHFEmaxMinus[iHistName]);

            hEvtMult[Parameters::HISTO_NAMES.size()]->Add(hEvtMult[iHistName]);

            hnVtx[Parameters::HISTO_NAMES.size()]->Add(hnVtx[iHistName]);
            hVX[Parameters::HISTO_NAMES.size()]->Add(hVX[iHistName]);
            hVY[Parameters::HISTO_NAMES.size()]->Add(hVY[iHistName]);
            hVZ[Parameters::HISTO_NAMES.size()]->Add(hVZ[iHistName]);

            hDchi2cl[Parameters::HISTO_NAMES.size()]->Add(hDchi2cl[iHistName]);
            hDalpha[Parameters::HISTO_NAMES.size()]->Add(hDalpha[iHistName]);
            hDdtheta[Parameters::HISTO_NAMES.size()]->Add(hDdtheta[iHistName]);
            hDmass[Parameters::HISTO_NAMES.size()]->Add(hDmass[iHistName]);
            hDsvpvSig[Parameters::HISTO_NAMES.size()]->Add(hDsvpvSig[iHistName]);
            hDsvpvDistance[Parameters::HISTO_NAMES.size()]->Add(hDsvpvDistance[iHistName]);
            hDsvpvDisErr[Parameters::HISTO_NAMES.size()]->Add(hDsvpvDisErr[iHistName]);

            hDtrk1Pt_vs_Dtrk2Pt[Parameters::HISTO_NAMES.size()]->Add(hDtrk1Pt_vs_Dtrk2Pt[iHistName]);
            hDpt_vs_Dy[Parameters::HISTO_NAMES.size()]->Add(hDpt_vs_Dy[iHistName]);

            hDtrk1Pt[Parameters::HISTO_NAMES.size()]->Add(hDtrk1Pt[iHistName]);
            hDtrk2Pt[Parameters::HISTO_NAMES.size()]->Add(hDtrk2Pt[iHistName]);
            hDpt[Parameters::HISTO_NAMES.size()]->Add(hDpt[iHistName]);
            hDy[Parameters::HISTO_NAMES.size()]->Add(hDy[iHistName]);
      }      
  }
}

bool checkPlottingModes(std::vector<string> modes, bool& needMergedHistograms, std::vector<int> &modeIndexes) {
      for (int i = 0; i < modes.size(); i++) {
            if (std::find(Parameters::HISTO_NAMES.begin(), Parameters::HISTO_NAMES.end(), modes[i]) == Parameters::HISTO_NAMES.end()) {
                  if (modes[i] == "Merged") {
                        needMergedHistograms = true;
                        modeIndexes.push_back(Parameters::HISTO_NAMES.size());
                  } else {
                        return false;
                  }
            } else {
                  modeIndexes.push_back((int)(std::find(Parameters::HISTO_NAMES.begin(), Parameters::HISTO_NAMES.end(), modes[i])- Parameters::HISTO_NAMES.begin()));
            }
      }
      return true;
}
