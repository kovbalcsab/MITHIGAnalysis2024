#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TFile.h"

#include "CommandLine.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> // For std::max_element

using namespace std;

bool checkFileType(string fileType);

void plotGraph(TH1* inHist, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const int nDim, const bool logx, const bool logy,
                 const double xLowLim, const double xHighLim, const double yLowLim, const double yHighLim, double latexY=0.27);

void plotOverlapGraph(TH1D* inHist, TH1D* inHist2, string fileType, string fileType2, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy, const double xLowLim, const double xHighLim,
                 double legendXLow=0.3, double legendXHigh=0.6, double legendYLow=0.7, double legendYHigh=0.9, double latexY=0.77);

void readInHistogramsFromFile(const char* histoSubname, TFile* inFile, TH2D* &hHFEmaxPlus_vs_EvtMult, TH2D* &hHFEmaxMinus_vs_EvtMult, TH1D* &hHFEmaxMinus, TH1D* &hHFEmaxPlus, TH1D* &hEvtMult, TH1D* &hnVtx, TH1D* &hVX, TH1D* &hVY, TH1D* &hVZ, TH1D* &hDchi2cl, TH1D* &hDalpha, TH1D* &hDdtheta, TH1D* &hDsvpvDistance, TH1D* &hDsvpvDisErr, TH1D* &hDsvpvSig, TH1D* &hDmass, TH2D* &hDtrk1Pt_vs_Dtrk2Pt, TH2D* &hDpt_vs_Dy, TH1D* &hDtrk1Pt, TH1D* &hDtrk2Pt, TH1D* &hDpt, TH1D* &hDy, TH2D* &htrkPt_vs_trkEta, TH1D* &htrkPt, TH1D* &htrkEta);

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
  string fileType = CL.Get("fileType", ""); // determine wether the file is Data, MC or MC_inclusive
  string fileType2 = CL.Get("fileType2", ""); // give second fileType only if want to overlap two plots /w ratio determine wether the file is Data, MC or MC_inclusive

  if (!checkFileType(fileType) || (fileType2 != "" && !checkFileType(fileType2) )) { std::cout << "[Error] Invalid fileType entered in config!" << std::endl; return -1; }

  int HFEMax = CL.GetInt("HFEMax",-999); // Read in HFEMax values for plotting

  string inputFileName;
  if (HFEMax == -999) {
    inputFileName = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, fileType.c_str() );
  } else {
    inputFileName = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o_DoSystRapGap%d/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, (int) HFEMax, fileType.c_str() );
  }


  TFile *inFile = new TFile(inputFileName.c_str());

  // Event level variables:
  TH2D *hHFEmaxPlus_vs_EvtMult; // gap energy vs event multiplicity
  TH2D *hHFEmaxMinus_vs_EvtMult; // gap energy vs event multiplicity
  TH1D *hHFEmaxMinus, *hHFEmaxPlus, *hEvtMult; // projections of the above histograms
  TH1D *hnVtx, *hVX, *hVY, *hVZ; // number of verticies and best vertex positions

  // Track level variables:
  TH2D *htrkPt_vs_trkEta;
  TH1D *htrkPt, *htrkEta;

  // D level variables: 

  TH1D *hDchi2cl, *hDalpha, *hDdtheta, *hDsvpvDistance, *hDsvpvDisErr, *hDsvpvSig, *hDmass;
  TH2D *hDtrk1Pt_vs_Dtrk2Pt, *hDpt_vs_Dy;
  TH1D *hDtrk1Pt, *hDtrk2Pt, *hDpt, *hDy;

  readInHistogramsFromFile("_1", inFile, hHFEmaxPlus_vs_EvtMult, hHFEmaxMinus_vs_EvtMult, hHFEmaxMinus, hHFEmaxPlus, hEvtMult, hnVtx, hVX, hVY, hVZ, hDchi2cl, hDalpha, hDdtheta, hDsvpvDistance, hDsvpvDisErr, hDsvpvSig, hDmass, hDtrk1Pt_vs_Dtrk2Pt, hDpt_vs_Dy, hDtrk1Pt, hDtrk2Pt, hDpt, hDy, htrkPt_vs_trkEta, htrkPt, htrkEta);

  /////////////////////////////////
  // 1. Plot the HF energy distributions
  /////////////////////////////////

  const string latexText = Form("%d < D_{p_{T}} < %d (GeV/#it{c})", (int) MinDzeroPT, (int) MaxDzeroPT);
  const string latexText2 = Form("%d < D_{y} < %d", (int) MinDzeroY, (int) MaxDzeroY);

  if (fileType2=="") {
    // Start plotting event level histograms
    plotGraph(hHFEmaxPlus_vs_EvtMult, "Charged hadron multiplicity", "HF E_{max} Plus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxPlus_vs_EvtMult/HFEmaxPlus_vs_EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 2, false, true, 0, 20, 0, 20);
    
    plotGraph(hHFEmaxMinus_vs_EvtMult, "Charged hadron multiplicity", "HF E_{max} Minus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxMinus_vs_EvtMult/HFEmaxMinus_vs_EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 2, false, true, 0, 20, 0, 20);

    plotGraph(hHFEmaxPlus, "Normalized counts", "HF E_{max} Plus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxPlus/HFEmaxPlus_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 20, 0, 0);

    plotGraph(hHFEmaxMinus, "Normalized counts", "HF E_{max} Minus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxMinus/HFEmaxMinus_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 20, 0, 0);
                  
    plotGraph(hEvtMult, "Normalized counts", "Charged hadron multiplicity",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/EvtMult/EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 80, 0, 0);

    plotGraph(hnVtx, "Normalized counts", "Number of vertices",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/nVtx/nVtx_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, -1, 3, 0, 0);
    
    plotGraph(hVX, "Normalized counts", "VX [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VX/VX_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, -1, 1, 0, 0);

    plotGraph(hVY, "Normalized counts", "VY [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VY/VY_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, -1, 1, 0, 0);

    plotGraph(hVZ, "Normalized counts", "VZ [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VZ/VZ_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, -15, 15, 0, 0);

    // Start plotting track level histograms

    if (htrkPt_vs_trkEta != nullptr) {
      plotGraph(htrkPt_vs_trkEta, "track p_{T} [GeV]", "track #eta",
              latexText.c_str(), latexText2.c_str(),
                  Form("%s/trkPt_vs_trkEta/trkPt_vs_trkEta_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), 2, false, false, 0, 15, -2.4, 2.4);

      plotGraph(htrkPt, "Normalized counts", "track p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
                  Form("%s/trkPt/trkPt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), 1, false, true, 0, 15, 0, 0);

      plotGraph(htrkEta, "Normalized counts", "track #eta",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/trkEta/trkEta_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), 1, false, true, -2.4, 2.4, 0, 0);
    }

    // Start plotting D level histograms

    plotGraph(hDalpha, "Normalized counts", "D^{0} #alpha",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dalpha/Dalpha_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 0.5, 0, 0);
    
    plotGraph(hDchi2cl, "Normalized counts", "D^{0} #chi^{2}",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dchi2cl/Dchi2cl_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 1., 0, 0);

    plotGraph(hDdtheta, "Normalized counts", "D^{0} #Delta#theta",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Ddtheta/Ddtheta_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 0.6, 0, 0);

    plotGraph(hDsvpvDisErr, "Normalized counts", "DsvpvDisErr",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDisErr/DsvpvDisErr_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 1.2, 0, 0);

    plotGraph(hDsvpvDistance, "Normalized counts", "DsvpvDistance",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDistance/DsvpvDistance_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 5, 0, 0);

    plotGraph(hDsvpvSig, "Normalized counts", "DsvpvSig",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/DsvpvSig/DsvpvSig_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str()), 1, false, true, 0, 15, 0, 0);

    plotGraph(hDmass, "Normalized counts", "D^{0} mass [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dmass/Dmass_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, false, 1.7, 2., 0, 0);

    plotGraph(hDtrk1Pt_vs_Dtrk2Pt, "D^{0} trk1 p_{T} [GeV]", "D^{0} trk2 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk1Pt_vs_Dtrk2Pt/Dtrk1Pt_vs_Dtrk2Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 2, false, false, 0, 5, 0, 5);
    
    plotGraph(hDpt_vs_Dy, "D^{0} p_{T} [GeV]", "D^{0} rapidity",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dpt_vs_Dy/Dpt_vs_Dy_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 2, false, false, 200, -200, 200, -200);

    plotGraph(hDtrk1Pt, "Normalized counts", "D^{0} trk1 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk1Pt/Dtrk1Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 5, 0, 0);

    plotGraph(hDtrk2Pt, "Normalized counts", "D^{0} trk2 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk2Pt/Dtrk2Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, true, 0, 5, 0, 0);

    plotGraph(hDpt, "Normalized counts", "D^{0} p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dpt/Dpt_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, false, 200, -200, 0, 0);

    plotGraph(hDy, "Normalized counts", "D^{0} y [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dy/Dy_pt%d-%d_y%d-%d_IsGammaN%o_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str()), 1, false, false, 200, -200, 0, 0);
    
  } else {
    string inputFileName2;
    if (HFEMax == -999) {
      inputFileName2 = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, fileType2.c_str() );
    } else {
      inputFileName2 = Form("%s/pt%d-%d_y%d-%d_IsGammaN%o_DoSystRapGap%d/%s.root", inputDir.c_str(), (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, (int) HFEMax, fileType2.c_str() );
    }
    TFile *inFile2 = new TFile(inputFileName2.c_str());
    
    // Event level variables:
    TH2D *hHFEmaxPlus_vs_EvtMult_2; // gap energy vs event multiplicity
    TH2D *hHFEmaxMinus_vs_EvtMult_2; // gap energy vs event multiplicity
    TH1D *hHFEmaxMinus_2, *hHFEmaxPlus_2, *hEvtMult_2; // projections of the above histograms
    TH1D *hnVtx_2, *hVX_2, *hVY_2, *hVZ_2; // number of verticies and best vertex positions

    // Track level variables:
    TH2D *htrkPt_vs_trkEta_2;
    TH1D *htrkPt_2, *htrkEta_2;

    // D level variables: 

    TH1D *hDchi2cl_2, *hDalpha_2, *hDdtheta_2, *hDsvpvDistance_2, *hDsvpvDisErr_2, *hDsvpvSig_2, *hDmass_2;
    TH2D *hDtrk1Pt_vs_Dtrk2Pt_2, *hDpt_vs_Dy_2;
    TH1D *hDtrk1Pt_2, *hDtrk2Pt_2, *hDpt_2, *hDy_2;

    readInHistogramsFromFile("_2", inFile2, hHFEmaxPlus_vs_EvtMult_2, hHFEmaxMinus_vs_EvtMult_2, hHFEmaxMinus_2, hHFEmaxPlus_2, hEvtMult_2, hnVtx_2, hVX_2, hVY_2, hVZ_2, hDchi2cl_2, hDalpha_2, hDdtheta_2, hDsvpvDistance_2, hDsvpvDisErr_2, hDsvpvSig_2, hDmass_2, hDtrk1Pt_vs_Dtrk2Pt_2, hDpt_vs_Dy_2, hDtrk1Pt_2, hDtrk2Pt_2, hDpt_2, hDy_2, htrkPt_vs_trkEta_2, htrkPt_2, htrkEta_2);

    // Start plotting event level histograms

    plotOverlapGraph(hHFEmaxPlus, hHFEmaxPlus_2, fileType, fileType2, "Normalized counts", "HF E_{max} Plus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxPlus/HFEmaxPlus_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 20);

    plotOverlapGraph(hHFEmaxMinus, hHFEmaxMinus_2, fileType, fileType2, "Normalized counts", "HF E_{max} Minus [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/HFEmaxMinus/HFEmaxMinus_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 20);
                  
    plotOverlapGraph(hEvtMult, hEvtMult_2, fileType, fileType2, "Normalized counts", "Charged hadron multiplicity",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/EvtMult/EvtMult_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 80);

    plotOverlapGraph(hnVtx, hnVtx_2, fileType, fileType2, "Normalized counts", "Number of vertices",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/nVtx/nVtx_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -1, 3, 0.15, 0.35);
    
    plotOverlapGraph(hVX, hVX_2, fileType, fileType2, "Normalized counts", "VX [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VX/VX_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -1, 1);

    plotOverlapGraph(hVY, hVY_2, fileType, fileType2, "Normalized counts", "VY [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VY/VY_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -1, 1);

    plotOverlapGraph(hVZ, hVZ_2, fileType, fileType2, "Normalized counts", "VZ [cm]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/VZ/VZ_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -15, 15, 0.15, 0.35);

    // Start plotting track level histograms

    if (htrkPt_vs_trkEta != nullptr && htrkPt_vs_trkEta_2 != nullptr) {

      plotOverlapGraph(htrkPt, htrkPt_2, fileType, fileType2, "Normalized counts", "track p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/trkPt/trkPt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 15);

      plotOverlapGraph(htrkEta, htrkEta_2, fileType, fileType2, "Normalized counts", "track #eta",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/trkEta/trkEta_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, -2.4, 2.4, 0.3, 0.6, 0.2, 0.4, 0.2);
                    
    }


    // Start plotting D level histograms

    plotOverlapGraph(hDalpha, hDalpha_2, fileType, fileType2, "Normalized counts", "D^{0} #alpha",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dalpha/Dalpha_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 0.5, 0.6, 0.8, 0.2, 0.4);
    
    plotOverlapGraph(hDchi2cl, hDchi2cl_2, fileType, fileType2, "Normalized counts", "D^{0} #chi^{2}",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dchi2cl/Dchi2cl_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 1, 0.3, 0.6, 0.2, 0.4, 0.1);

    plotOverlapGraph(hDdtheta, hDdtheta_2, fileType, fileType2, "Normalized counts", "D^{0} #Delta#theta",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Ddtheta/Ddtheta_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 0.6, 0.2, 0.4, 0.2, 0.4);

    plotOverlapGraph(hDsvpvDisErr, hDsvpvDisErr_2, fileType, fileType2, "Normalized counts", "DsvpvDisErr",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDisErr/DsvpvDisErr_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 1.2);

    plotOverlapGraph(hDsvpvDistance, hDsvpvDistance_2, fileType, fileType2, "Normalized counts", "DsvpvDistance",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/DsvpvDistance/DsvpvDistance_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 5);

    plotOverlapGraph(hDsvpvSig, hDsvpvSig_2, fileType, fileType2, "Normalized counts", "DsvpvSig",
                  latexText.c_str(), latexText2.c_str(),
                  Form("%s/DsvpvSig/DsvpvSig_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                        PlotDir.c_str(),
                        (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                        IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 15, 0.35, 0.6);

    plotOverlapGraph(hDmass, hDmass_2, fileType, fileType2, "Normalized counts", "D^{0} mass [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dmass/Dmass_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 1.7, 2);

    plotOverlapGraph(hDtrk1Pt, hDtrk1Pt_2, fileType, fileType2, "Normalized counts", "D^{0} trk1 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk1Pt/Dtrk1Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 5, 0.3, 0.6, 0.2, 0.4);

    plotOverlapGraph(hDtrk2Pt, hDtrk2Pt_2, fileType, fileType2, "Normalized counts", "D^{0} trk2 p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dtrk2Pt/Dtrk2Pt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, true, 0, 5, 0.3, 0.6, 0.2, 0.4);

    plotOverlapGraph(hDpt, hDpt_2, fileType, fileType2, "Normalized counts", "D^{0} p_{T} [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dpt/Dpt_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 200, -200);

    plotOverlapGraph(hDy, hDy_2, fileType, fileType2, "Normalized counts", "D^{0} y [GeV]",
              latexText.c_str(), latexText2.c_str(),
              Form("%s/Dy/Dy_pt%d-%d_y%d-%d_IsGammaN%o_%s_over_%s.pdf",
                    PlotDir.c_str(),
                    (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                    IsGammaN, fileType.c_str(), fileType2.c_str()), false, false, 200, -200);
  }

  return 0;
}


void plotGraph(TH1* inHist, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const int nDim, const bool logx, const bool logy, 
                 const double xLowLim, const double xHighLim, const double yLowLim, const double yHighLim, double latexY)
{
  // Create canvas
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLeftMargin(0.13);
  if (nDim==1) {
    canvas->SetRightMargin(0.04);
  } else if (nDim==2) {
    canvas->SetRightMargin(0.12);
  }
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.08);
  if (logx) { gPad->SetLogx(); }
  if (logy) { gPad->SetLogy(); }
  if (nDim==2) {
    gPad->SetLogz();
  }

  // Create and configure the histogram frame
  inHist->GetYaxis()->SetTitle(yAxisTitle);
  inHist->GetXaxis()->SetTitle(xAxisTitle);

  int firstNonEmptyBin = inHist->FindFirstBinAbove(0,1);
  int lastNonEmptyBin = inHist->FindLastBinAbove(0,1);
  inHist->GetXaxis()->SetRangeUser(TMath::Min(xLowLim, inHist->GetXaxis()->GetBinLowEdge(firstNonEmptyBin)), TMath::Max(xHighLim, inHist->GetXaxis()->GetBinUpEdge(lastNonEmptyBin)));

  if (nDim==2) {
    int firstNonEmptyBiny = inHist->FindFirstBinAbove(0,2);
    int lastNonEmptyBiny = inHist->FindLastBinAbove(0,2);
    inHist->GetYaxis()->SetRangeUser(TMath::Min(yLowLim, inHist->GetYaxis()->GetBinLowEdge(firstNonEmptyBiny)), TMath::Max(yHighLim, inHist->GetYaxis()->GetBinUpEdge(lastNonEmptyBiny))); 
  }

  inHist->SetStats(0);
  inHist->SetTitle("");
  inHist->GetYaxis()->SetTitleOffset(1.5);
  if (nDim==1) {
    inHist->Draw();
  } else if (nDim==2) {
    inHist->Draw("colz");
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
  canvas->Update();
  canvas->SaveAs(plotname);

  // Clean up
  delete canvas;
}

void plotOverlapGraph(TH1D* inHist, TH1D* inHist2, string fileType, string fileType2, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname, const bool logx, const bool logy, const double xLowLim, const double xHighLim,
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
  

  auto styleUpHist = [](TH1 *hist) {
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->ChangeLabel(1,0,0);
    hist->SetTitle("");
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

void readInHistogramsFromFile(const char* histoSubname, TFile* inFile, TH2D* &hHFEmaxPlus_vs_EvtMult, TH2D* &hHFEmaxMinus_vs_EvtMult, TH1D* &hHFEmaxMinus, TH1D* &hHFEmaxPlus, TH1D* &hEvtMult, TH1D* &hnVtx, TH1D* &hVX, TH1D* &hVY, TH1D* &hVZ, TH1D* &hDchi2cl, TH1D* &hDalpha, TH1D* &hDdtheta, TH1D* &hDsvpvDistance, TH1D* &hDsvpvDisErr, TH1D* &hDsvpvSig, TH1D* &hDmass, TH2D* &hDtrk1Pt_vs_Dtrk2Pt, TH2D* &hDpt_vs_Dy, TH1D* &hDtrk1Pt, TH1D* &hDtrk2Pt, TH1D* &hDpt, TH1D* &hDy, TH2D* &htrkPt_vs_trkEta, TH1D* &htrkPt, TH1D* &htrkEta) {
  hHFEmaxPlus_vs_EvtMult = (TH2D*) inFile->Get("hHFEmaxPlus_vs_EvtMult");
  hHFEmaxPlus_vs_EvtMult->SetName(Form("hHFEmaxPlus_vs_EvtMult%s",histoSubname));
  hHFEmaxPlus_vs_EvtMult->Scale(1./hHFEmaxPlus_vs_EvtMult->Integral());
  hHFEmaxPlus = (TH1D*) hHFEmaxPlus_vs_EvtMult->ProjectionX(Form("hHFEmaxPlus%s",histoSubname));
  hHFEmaxPlus->Rebin(2);

  hHFEmaxMinus_vs_EvtMult = (TH2D*) inFile->Get("hHFEmaxMinus_vs_EvtMult");
  hHFEmaxMinus_vs_EvtMult->SetName(Form("hHFEmaxMinus_vs_EvtMult%s",histoSubname));
  hHFEmaxMinus_vs_EvtMult->Scale(1./hHFEmaxMinus_vs_EvtMult->Integral());
  hHFEmaxMinus = (TH1D*) hHFEmaxMinus_vs_EvtMult->ProjectionX(Form("hHFEmaxMinus%s",histoSubname));
  hHFEmaxMinus->Rebin(2);

  hEvtMult = (TH1D*) hHFEmaxPlus_vs_EvtMult->ProjectionY(Form("hEvtMult%s",histoSubname));

  hnVtx = (TH1D*) inFile->Get("hnVtx");
  hnVtx->SetName(Form("hnVtx%s",histoSubname));
  hnVtx->Scale(1./hnVtx->Integral());
  hVX = (TH1D*) inFile->Get("hVX");
  hVX->SetName(Form("hVX%s",histoSubname));
  hVX->Scale(1./hVX->Integral());
  hVY = (TH1D*) inFile->Get("hVY");
  hVY->SetName(Form("hVY%s",histoSubname));
  hVY->Scale(1./hVY->Integral());
  hVZ = (TH1D*) inFile->Get("hVZ");
  hVZ->SetName(Form("hVZ%s",histoSubname));
  hVZ->Rebin(2);
  hVZ->Scale(1./hVZ->Integral());

  htrkPt_vs_trkEta = (TH2D*) inFile->Get("htrkPt_vs_trkEta");
  htrkPt = nullptr;
  htrkEta = nullptr;
  if (htrkPt_vs_trkEta != nullptr) {
      htrkPt_vs_trkEta->SetName(Form("htrkPt_vs_trkEta%s",histoSubname));
      htrkPt_vs_trkEta->Scale(1./htrkPt_vs_trkEta->Integral());
      htrkPt = (TH1D*) htrkPt_vs_trkEta->ProjectionX(Form("htrkPt%s",histoSubname));
      htrkEta = (TH1D*) htrkPt_vs_trkEta->ProjectionY(Form("htrkEta%s",histoSubname));
  }

  hDchi2cl = (TH1D*) inFile->Get("hDchi2cl");
  hDchi2cl->SetName(Form("hDchi2cl%s",histoSubname));
  hDchi2cl->Rebin(2);
  hDchi2cl->Scale(1./hDchi2cl->Integral());
  hDalpha = (TH1D*) inFile->Get("hDalpha");
  hDalpha->SetName(Form("hDalpha%s",histoSubname));
  hDalpha->Scale(1./hDalpha->Integral());
  hDdtheta = (TH1D*) inFile->Get("hDdtheta");
  hDdtheta->SetName(Form("hDdtheta%s",histoSubname));
  hDdtheta->Scale(1./hDdtheta->Integral());
  hDmass = (TH1D*) inFile->Get("hDmass");
  hDmass->SetName(Form("hDmass%s",histoSubname));
  hDmass->Rebin(2);
  hDmass->Scale(1./hDmass->Integral());
  hDsvpvSig = (TH1D*) inFile->Get("hDsvpvSig");
  hDsvpvSig->SetName(Form("hDsvpvSig%s",histoSubname));
  hDsvpvSig->Rebin(2);
  hDsvpvSig->Scale(1./hDsvpvSig->Integral());
  hDsvpvDistance = (TH1D*) inFile->Get("hDsvpvDistance");
  hDsvpvDistance->SetName(Form("hDsvpvDistance%s",histoSubname));
  hDsvpvDistance->Rebin(2);
  hDsvpvDistance->Scale(1./hDsvpvDistance->Integral());
  hDsvpvDisErr = (TH1D*) inFile->Get("hDsvpvDisErr");
  hDsvpvDisErr->SetName(Form("hDsvpvDisErr%s",histoSubname));
  hDsvpvDisErr->Rebin(2);
  hDsvpvDisErr->Scale(1./hDsvpvDisErr->Integral());

  hDtrk1Pt_vs_Dtrk2Pt = (TH2D*) inFile->Get("hDtrk1Pt_vs_Dtrk2Pt");
  hDtrk1Pt_vs_Dtrk2Pt->SetName(Form("hDtrk1Pt_vs_Dtrk2Pt%s",histoSubname));
  hDtrk1Pt_vs_Dtrk2Pt->Scale(1./hDtrk1Pt_vs_Dtrk2Pt->Integral());
  hDpt_vs_Dy = (TH2D*) inFile->Get("hDpt_vs_Dy");
  hDpt_vs_Dy->SetName(Form("hDpt_vs_Dy%s",histoSubname));
  hDpt_vs_Dy->Scale(1./hDpt_vs_Dy->Integral());

  hDtrk1Pt = (TH1D*) hDtrk1Pt_vs_Dtrk2Pt->ProjectionX(Form("hDtrk1Pt%s",histoSubname));
  hDtrk2Pt = (TH1D*) hDtrk1Pt_vs_Dtrk2Pt->ProjectionY(Form("hDtrk2Pt%s",histoSubname));
  hDpt = (TH1D*) hDpt_vs_Dy->ProjectionX(Form("hDpt%s",histoSubname));
  hDpt->Rebin(4);
  hDy = (TH1D*) hDpt_vs_Dy->ProjectionY(Form("hDy%s",histoSubname));
  hDy->Rebin(2);
}
