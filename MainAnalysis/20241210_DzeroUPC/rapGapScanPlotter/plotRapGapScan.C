#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"

#include "CommandLine.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> // For std::max_element

#include "plotRapGapScan.h"

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

  vector<int> HFEMax = CL.GetIntVector("HFEMax",""); // Read in HFEMax values for plotting

  vector<string> inputPoints; // Input corrected yields md files
  vector<string> inputPointPaths; // Input fit olders
  vector<double> yminVec;
  vector<double> ymaxVec;

  for (int iFile=0; iFile < HFEMax.size(); iFile++) {
    inputPoints.push_back( (string) Form("rapGapScan_threshold_%d/pt%d-%d_y%d-%d_IsGammaN%o/MassFit/correctedYields.md", (int) HFEMax[iFile], (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN ));
    inputPointPaths.push_back( (string) Form("rapGapScan_threshold_%d/pt%d-%d_y%d-%d_IsGammaN%o/MassFit/", (int) HFEMax[iFile], (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY, IsGammaN ));
    yminVec.push_back(MinDzeroY);
    ymaxVec.push_back(MaxDzeroY);
  }

  /////////////////////////////////
  // 0. Extract the points from the vector of .md
  /////////////////////////////////

  // nominal central values
  const int nPoints = inputPoints.size();
  std::vector<Point> PointsArr = getPointArr(MinDzeroPT, MaxDzeroPT, IsGammaN, HFEMax, inputPoints);
  std::vector<fitPoint> FitPointsArr = getFitPointArr(MinDzeroPT, MaxDzeroPT, IsGammaN, HFEMax, yminVec, ymaxVec, inputPointPaths);

  vector<double> HFEMaxValues = getDoubleArr(PointsArr, 
                           [](Point& p) -> double { return p.HFEMax;} );
  vector<double> HFEMaxErrors = getDoubleArr(PointsArr, 
                           [](Point& p) -> double { return 0.5;} );
  vector<double> correctedYieldValues = getDoubleArr(PointsArr, 
                           [](Point& p) -> double { return p.correctedYield;} );
  vector<double> correctedYieldErrors = getDoubleArr(PointsArr, 
                           [](Point& p) -> double { return p.correctedYieldError;} );
  vector<double> rawYieldValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.rawYield;} );
  vector<double> rawYieldErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.rawYieldError;} );
  vector<double> effEvtValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.effEvt;} );
  vector<double> effEvtErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.effEvtError;} );
  vector<double> numEvtValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.numEvt;} );
  vector<double> numEvtErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return TMath::Sqrt(p.numEvt);} );
  vector<double> denEvtValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.denEvt;} );
  vector<double> denEvtErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return TMath::Sqrt(p.denEvt);} );
  vector<double> effDValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.effD;} );
  vector<double> effDErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.effDError;} );
  vector<double> numDValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.numD;} );
  vector<double> numDErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return TMath::Sqrt(p.numD);} );
  vector<double> denDValues = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return p.denD;} );
  vector<double> denDErrors = getDoubleArr(PointsArr,
                           [](Point& p) -> double { return TMath::Sqrt(p.denD);} );

  vector<double> fitLambdaValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.lambda;} );
  vector<double> fitLambdaErrors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.lambdaError;} );
  vector<double> fitSignalAlphaValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalAlpha;} );
  vector<double> fitSignalAlphaErrors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalAlphaError;} );
  vector<double> fitSignalFractionValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalFraction;} );
  vector<double> fitSignalFractionErrors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalFractionError;} );
  vector<double> fitSignalMeanValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalMean;} );
  vector<double> fitSignalMeanErrors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalMeanError;} );
  vector<double> fitSignalSigma1Values = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalSigma1;} );
  vector<double> fitSignalSigma1Errors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalSigma1Error;} );
  vector<double> fitSignalSigma2Values = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalSigma2;} );
  vector<double> fitSignalSigma2Errors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.signalSigma2Error;} );
  vector<double> fitSwapMeanValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.swapMean;} );
  vector<double> fitSwapMeanErrors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.swapMeanError;} );
  vector<double> fitSwapSigmaValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.swapSigma;} );
  vector<double> fitSwapSigmaErrors = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.swapSigmaError;} );
  vector<double> fitSwpFractionValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.swpFraction;} );
  vector<double> fitPkppFractionValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.pkppFraction;} );
  vector<double> fitPkkkFractionValues = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return p.pkkkFraction;} );
  vector<double> nullVec = getDoubleArr(FitPointsArr,
                           [](fitPoint& p) -> double { return 0;} );

  printArr(correctedYieldValues, ", ", "correctedYieldValues: ");

  /////////////////////////////////
  // 1. Plot the cross section
  /////////////////////////////////
  // Create a canvas
  TCanvas* c1 = new TCanvas("c1", "D0 Cross Section", 800, 800);
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.04);
  c1->SetBottomMargin(0.12);
  c1->SetTopMargin(0.08);

  TH1F* hFrame = new TH1F("hFrame", " ", 100, 0, 32);
  hFrame->GetYaxis()->SetTitle("d^{2}#sigma/dydp_{T} (mb/GeV)");
  hFrame->GetXaxis()->SetTitle("HF Energy threshold [GeV]");
  hFrame->SetStats(0);
  hFrame->GetYaxis()->SetTitleOffset(1.5);
  hFrame->GetYaxis()->SetRangeUser(0, (*std::max_element(correctedYieldValues.begin(), correctedYieldValues.end()))*1.2);
  hFrame->Draw();

  TGraphErrors* gr = new TGraphErrors(nPoints, HFEMaxValues.data(), correctedYieldValues.data(), HFEMaxErrors.data(), correctedYieldErrors.data());
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);
  gr->SetLineColor(kRed);
  gr->SetMarkerColor(kRed);
  gr->SetLineWidth(2);

  gr->Draw("P E1 SAME");

  /*TLegend* leg = new TLegend(0.1, 0.15, 0.48, 0.28);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(gr, "New framework", "P");
  leg->Draw();*/

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.SetTextFont(42);
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.035);
  latex2.SetTextFont(42);
  // latex.DrawLatex(0.15, 0.92, "CMS #it{Preliminary} 1.38 nb^{-1} (5.36 TeV PbPb)");
  // latex.DrawLatex(0.15, 0.86, "UPCs, ZDC Xn0n w/ gap");
  // latex.DrawLatex(0.15, 0.82, "Global uncert. #pm 5.05%");
  latex.DrawLatex(0.6, 0.32, Form("%d < D_{p_{T}} < %d (GeV/#it{c})", (int) MinDzeroPT, (int) MaxDzeroPT));
  latex2.DrawLatex(0.6, 0.27, Form("%d < D_{y} < %d", (int) MinDzeroY, (int) MaxDzeroY));

  c1->Update();
  c1->SaveAs(Form("%s/rapGapScanCorrectedYields_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));
  delete gr;
  delete hFrame;

  /////////////////////////////////
  // 4. Plot other relevant plots
  /////////////////////////////////
  auto plotGraph = [](const char* yAxisTitle, const char* xAxisTitle,
                 double Emin, double Emax,
                 const std::vector<double>& xValues, const std::vector<double>& yValues,
                 const std::vector<double>& xErrors, const std::vector<double>& yErrors,
                 const char* latexText, const char* latexText2, const char* plotname,
                 int nBinsX=100, double xMin=0, double xMax=32)
  {
    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.04);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.08);

    // Create and configure the histogram frame
    TH1F* hFrame = new TH1F("hFrame", "", nBinsX, xMin, xMax);
    hFrame->GetYaxis()->SetTitle(yAxisTitle);
    hFrame->GetXaxis()->SetTitle(xAxisTitle);
    hFrame->SetStats(0);
    hFrame->GetYaxis()->SetTitleOffset(1.5);
    hFrame->GetYaxis()->SetRangeUser(Emin, Emax);
    hFrame->Draw();

    // Create TGraphErrors for data points
    TGraphErrors* graph = new TGraphErrors(xValues.size(), xValues.data(), yValues.data(),
                                           xErrors.data(), yErrors.data());
    graph->Draw("P E1 SAME");

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
    delete graph;
    delete hFrame;
    delete canvas;
  };

  const char* latexText = Form("%d < D p_{T} < %d (GeV/#it{c})", (int) MinDzeroPT, (int) MaxDzeroPT);
  const char* latexText2 = Form("%d < D y < %d", (int) MinDzeroY, (int) MaxDzeroY);

  plotGraph("#varepsilon_{event}", "HF Energy threshold [GeV]",
            0.8, 1.05,
            HFEMaxValues, effEvtValues, HFEMaxErrors, effEvtErrors,
            latexText, latexText2,
            Form("%s/evtEff_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Numerator N_{event}", "HF Energy threshold [GeV]",
            0, (*std::max_element(numEvtValues.begin(), numEvtValues.end()))*1.3,
            HFEMaxValues, numEvtValues, HFEMaxErrors, numEvtErrors,
            latexText, latexText2,
            Form("%s/evtNum_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Denominator N_{event}", "HF Energy threshold [GeV]",
            0, (*std::max_element(denEvtValues.begin(), denEvtValues.end()))*1.3,
            HFEMaxValues, denEvtValues, HFEMaxErrors, denEvtErrors,
            latexText, latexText2,
            Form("%s/evtDen_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("#varepsilon_{D}", "HF Energy threshold [GeV]",
            0, 1.05,
            HFEMaxValues, effDValues, HFEMaxErrors, effDErrors,
            latexText, latexText2,
            Form("%s/DEff_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("#varepsilon_{D}", "HF Energy threshold [GeV]",
            0, 0.2, //(*std::max_element(effDValues.begin(), effDValues.end()))*1.3,
            HFEMaxValues, effDValues, HFEMaxErrors, effDErrors,
            latexText, latexText2,
            Form("%s/DEff_zoom_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Numerator N_{D}", "HF Energy threshold [GeV]",
            0, (*std::max_element(numDValues.begin(), numDValues.end()))*1.3,
            HFEMaxValues, numDValues, HFEMaxErrors, numDErrors,
            latexText, latexText2,
            Form("%s/DNum_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Denominator N_{D}", "HF Energy threshold [GeV]",
            0, (*std::max_element(denDValues.begin(), denDValues.end()))*1.3,
            HFEMaxValues, denDValues, HFEMaxErrors, denDErrors,
            latexText, latexText2,
            Form("%s/DDen_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Raw yield", "HF Energy threshold [GeV]",
            0, (*std::max_element(rawYieldValues.begin(), rawYieldValues.end()))*1.3,
            HFEMaxValues, rawYieldValues, HFEMaxErrors, rawYieldErrors,
            latexText, latexText2,
            Form("%s/RawYield_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Signal alpha", "HF Energy threshold [GeV]",
            0, (*std::max_element(fitSignalAlphaValues.begin(), fitSignalAlphaValues.end()))*1.3,
            HFEMaxValues, fitSignalAlphaValues, HFEMaxErrors, fitSignalAlphaErrors,
            latexText, latexText2,
            Form("%s/SignalAlpha_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Lambda", "HF Energy threshold [GeV]",
            std::min((std::min_element(fitLambdaValues.begin(), fitLambdaValues.end())[0])*1.3,0.), std::max((std::max_element(fitLambdaValues.begin(), fitLambdaValues.end())[0])*1.3,0.),
            HFEMaxValues, fitLambdaValues, HFEMaxErrors, fitLambdaErrors,
            latexText, latexText2,
            Form("%s/Lambda_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Signal fraction", "HF Energy threshold [GeV]",
            0, 1.05,
            HFEMaxValues, fitSignalFractionValues, HFEMaxErrors, fitSignalFractionErrors,
            latexText, latexText2,
            Form("%s/SignalFraction_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Signal mean", "HF Energy threshold [GeV]",
            1.7, 2.0,
            HFEMaxValues, fitSignalMeanValues, HFEMaxErrors, fitSignalMeanErrors,
            latexText, latexText2,
            Form("%s/SignalMean_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  const char* latexText3 = Form("%d < D p_{T} < %d (GeV/#it{c})", (int) MinDzeroPT, (int) MaxDzeroPT);
  const char* latexText4 = Form("%d < D y < %d", (int) MinDzeroY, (int) MaxDzeroY);

  plotGraph("Signal sigma1", "HF Energy threshold [GeV]",
            0, (*std::max_element(fitSignalSigma1Values.begin(), fitSignalSigma1Values.end()))*1.3,
            HFEMaxValues, fitSignalSigma1Values, HFEMaxErrors, fitSignalSigma1Errors,
            latexText3, latexText4,
            Form("%s/SignalSigma1_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Signal sigma2", "HF Energy threshold [GeV]",
            0, std::max_element(fitSignalSigma2Values.begin(), fitSignalSigma2Values.end())[0]*1.3,
            HFEMaxValues, fitSignalSigma2Values, HFEMaxErrors, fitSignalSigma2Errors,
            latexText3, latexText4,
            Form("%s/SignalSigma2_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));
  plotGraph("Swap mean", "HF Energy threshold [GeV]",
            1.7, 2.0,
            HFEMaxValues, fitSwapMeanValues, HFEMaxErrors, fitSwapMeanErrors,
            latexText3, latexText4,
            Form("%s/SwapMean_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Swap Sigma", "HF Energy threshold [GeV]",
            0, (*std::max_element(fitSwapSigmaValues.begin(), fitSwapSigmaValues.end()))*1.3,
            HFEMaxValues, fitSwapSigmaValues, HFEMaxErrors, fitSwapSigmaErrors,
            latexText3, latexText4,
            Form("%s/SwapSigma_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("Swapped Fraction", "HF Energy threshold [GeV]",
            0, std::max_element(fitSwpFractionValues.begin(), fitSwpFractionValues.end())[0]*1.3,
            HFEMaxValues, fitSwpFractionValues, HFEMaxErrors, nullVec,
            latexText3, latexText4,
            Form("%s/SwappedFraction_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("PK PP Fraction", "HF Energy threshold [GeV]",
            0, std::max_element(fitPkppFractionValues.begin(), fitPkppFractionValues.end())[0]*1.3,
            HFEMaxValues, fitPkppFractionValues, HFEMaxErrors, nullVec,
            latexText3, latexText4,
            Form("%s/PKPPFraction_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  plotGraph("PK KK Fraction", "HF Energy threshold [GeV]",
            0, std::max_element(fitPkkkFractionValues.begin(), fitPkkkFractionValues.end())[0]*1.3,
            HFEMaxValues, fitPkkkFractionValues, HFEMaxErrors, nullVec,
            latexText3, latexText4,
            Form("%s/PKKKFraction_pt%d-%d_y%d-%d_IsGammaN%o.pdf",
                  PlotDir.c_str(),
                  (int) MinDzeroPT, (int) MaxDzeroPT, (int) MinDzeroY, (int) MaxDzeroY,
                  IsGammaN));

  return 0;
}
