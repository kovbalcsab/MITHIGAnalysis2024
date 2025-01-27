#include "TCanvas.h"
#include "TH1D.h"
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

int main(int argc, char *argv[])
{
  CommandLine CL(argc, argv);
  string         PlotDir       = CL.Get    ("PlotDir", "");       // subdirectory under plot/
  float MinDzeroPT = CL.GetDouble("MinDzeroPT", 2);  // Minimum Dzero transverse momentum threshold for Dzero selection.
  float MaxDzeroPT = CL.GetDouble("MaxDzeroPT", 5);  // Maximum Dzero transverse momentum threshold for Dzero selection.
  float MinDzeroY = CL.GetDouble("MinDzeroY", 0);  // Minimum Dzero rapidity threshold for Dzero selection.
  float MaxDzeroY = CL.GetDouble("MaxDzeroY", 1);  // Maximum Dzero rapidity threshold for Dzero selection.
  bool IsGammaN = CL.GetBool("IsGammaN", true);      // GammaN analysis (or NGamma)
  string fileType = CL.Get("fileType", ""); // determine wether the file is Data, MC or MC_inclusive

  int HFEMax = CL.GetInt("HFEMax",-999); // Read in HFEMax values for plotting

  string inputFileName = Form("rapGapScan/pt%d-%d_y%d-%d_IsGammaN%o_DoSystRapGap%d/%s.root", (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN, (int) HFEMax, fileType.c_str() );

  TFile *inFile = new TFile(inputFileName.c_str());
  TH1D *hHFEmaxPlus = (TH1D*) inFile->Get("hHFEmaxPlus");
  hHFEmaxPlus->Scale(1./hHFEmaxPlus->Integral());
  TH1D *hHFEmaxMinus = (TH1D*) inFile->Get("hHFEmaxMinus");
  hHFEmaxMinus->Scale(1./hHFEmaxMinus->Integral());

  /////////////////////////////////
  // 1. Plot the HF energy distributions
  /////////////////////////////////

  auto plotGraph = [](TH1D* inHist, const char* yAxisTitle, const char* xAxisTitle,
                 const char* latexText, const char* latexText2, const char* plotname)
  {
    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.04);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.08);

    // Create and configure the histogram frame
    inHist->GetYaxis()->SetTitle(yAxisTitle);
    inHist->GetXaxis()->SetTitle(xAxisTitle);
    inHist->SetStats(0);
    inHist->GetYaxis()->SetTitleOffset(1.5);
    inHist->Draw();

    // Add TLatex for additional text
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.DrawLatex(0.6, 0.32, latexText);
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.035);
    latex2.SetTextFont(42);
    latex2.DrawLatex(0.6, 0.27, latexText2);

    // Update and save the canvas
    canvas->Update();
    canvas->SaveAs(plotname);

    // Clean up
    delete canvas;
  };

  const char* latexText = Form("%d < D_{p_{T}} < %d (GeV/#it{c})", (int) MinDzeroPT, (int) MaxDzeroPT);
  const char* latexText2 = Form("%d < D_{y} < %d", (int) MinDzeroY, (int) MaxDzeroY);

  plotGraph(hHFEmaxPlus, "Normalized counts", "HF E_{max} Plus [GeV]",
            latexText, latexText2,
            Form("%s/HFEmaxPlus_pt%d-%d_y%d-%d_IsGammaN%o_HFETh%d_%s.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN, (int) HFEMax, fileType.c_str()));

  plotGraph(hHFEmaxMinus, "Normalized counts", "HF E_{max} Minus [GeV]",
            latexText, latexText2,
            Form("%s/HFEmaxMinus_pt%d-%d_y%d-%d_IsGammaN%o_HFETh%d_%s.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN, (int) HFEMax, fileType.c_str()));

  return 0;
}
