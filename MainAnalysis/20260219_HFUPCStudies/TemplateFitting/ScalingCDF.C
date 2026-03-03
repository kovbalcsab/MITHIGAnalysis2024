#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooCustomizer.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF1.h"
#include <map>
#include <iostream>

#include "CommandLine.h" // Yi's Commandline bundle

using namespace RooFit;

int main(int argc, char** argv) {
    // Set up command line parser
    CommandLine CL(argc, argv);
    std::string targetFileName = CL.Get("TargetFile", "target.root");
    std::string fileToFitName = CL.Get("FileToFit", "fit.root");
    std::string outputFileName = CL.Get("OutputFileName", "fit_result.root");
    std::string varFitName = CL.Get("VarFitName", "HFEMaxPlus_forest");
    std::string varTargetName = CL.Get("VarTargetName", "HFEMaxPlus_forest");
    int NbinsTarget = CL.GetDouble("NBinsTarget", 1000);
    int NbinsFit = CL.GetDouble("NBinsFit", 1000);
    double EMaxTarget = CL.GetDouble("EMaxTarget", 250.0);
    double EMaxFit = CL.GetDouble("EMaxFit", 250.0);
    double plotXMax = CL.GetDouble("PlotXMax", 25.0);
    double plotYMax = CL.GetDouble("PlotYMax", 25.);


    // ---------------------------------------------------------
    // 1. Read in the TFiles and TTrees
    // ---------------------------------------------------------
    TFile* targetFile = new TFile(targetFileName.c_str(), "READ");
    if (!targetFile || targetFile->IsZombie()) {
        std::cerr << "Error opening target file: " << targetFileName << std::endl;
        return -1;
    }
    TTree* targetTree = dynamic_cast<TTree*>(targetFile->Get("OutputTree"));
    TH1D* targetHist = new TH1D("targetHist", "targetHist", NbinsTarget, 0, EMaxTarget);
    targetTree->Draw((varTargetName + ">>targetHist").c_str(), "", "goff");

    TFile* fileToFit = new TFile(fileToFitName.c_str(), "READ");
    if (!fileToFit || fileToFit->IsZombie()) {
        std::cerr << "Error opening file to fit: " << fileToFitName << std::endl;
        return -1;
    }
    TTree* fitTree = dynamic_cast<TTree*>(fileToFit->Get("OutputTree"));
    TH1D* fitHist = new TH1D("fitHist", "fitHist", NbinsFit, 0, EMaxFit);
    fitTree->Draw((varFitName + ">>fitHist").c_str(), "", "goff");

    // ---------------------------------------------------------
    // 2. Create CDF of histograms
    // ---------------------------------------------------------
    targetHist->Scale(1.0 / targetHist->Integral(""));
    fitHist->Scale(1.0 / fitHist->Integral(""));

    TH1D* CDF_MC= (TH1D*)fitHist->GetCumulative();

    // ---------------------------------------------------------
    // 3. Create E_MC points to evaluate the CDF
    // ---------------------------------------------------------
    TGraph* cdfGraph = new TGraph();
    double p[1];
    double xp[1];
    for (int i = 1; i <= CDF_MC->GetNbinsX(); ++i) {
        int n = cdfGraph->GetN();
        double E_MC = CDF_MC->GetBinCenter(i);
        double p_MC = CDF_MC->GetBinContent(i);
        p[0] = p_MC;
        targetHist->GetQuantiles(1, xp, p);
        cdfGraph->SetPoint(n, E_MC, xp[0]);
    }

    // ---------------------------------------------------------
    // 4. make polinom fit to the CDF points
    // ---------------------------------------------------------

    TF1* fitFunc = new TF1("fitFunc", "[2]*x*x*x + [1]*x*x + [0]*x + [3]", 1, 6);
    cdfGraph->Fit(fitFunc, "R");

     // Print fit parameters
    double a = fitFunc->GetParameter(0);
    double b = fitFunc->GetParameter(1);
    double c = fitFunc->GetParameter(2);
    double d = fitFunc->GetParameter(3);
    std::cout << "Fit parameters: a = " << a << ", b = " << b << ", c = " << c << ", d = " << d << std::endl;

    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        return -1;
    }
    outputFile->cd();
    cdfGraph->Write("cdfGraph");
    fitFunc->Write("fitFunc");

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();
    cdfGraph->SetMarkerStyle(20);
    cdfGraph->SetMarkerColor(kBlue);
    cdfGraph->SetLineColor(kBlue);
    cdfGraph->SetLineWidth(2);
    cdfGraph->GetXaxis()->SetRangeUser(0, plotXMax);
    cdfGraph->GetYaxis()->SetRangeUser(0, plotYMax);
    cdfGraph->SetTitle("CDF of MC;E_{MC} (GeV);E_{Data} (GeV)");
    cdfGraph->Draw("ALP");
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("Same");

    c1->SaveAs("CDF_MC_vs_Data.pdf");

    
    outputFile->Close();
    

    return 0;
}