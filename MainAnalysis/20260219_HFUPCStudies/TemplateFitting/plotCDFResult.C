#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include <map>
#include <iostream>

#include "CommandLine.h" // Yi's Commandline bundle

int main(int argc, char** argv) {
    // Set up command line parser
    CommandLine CL(argc, argv);
    std::string targetFileName = CL.Get("TargetFile", "target.root");
    std::string fileToFitName = CL.Get("FileToFit", "fit.root");
    std::string fitResultFileName = CL.Get("FitResultFileName", "fit_result.root");
    std::string outputFileName = CL.Get("OutputFileName", "fit_result_CDF.pdf");
    std::string varFitName = CL.Get("VarFitName", "HFEMaxPlus_forest");
    std::string varTargetName = CL.Get("VarTargetName", "HFEMaxPlus_forest");
    int NbinsTarget = CL.GetInt("NBinsTarget", 1000);
    int NbinsFit = CL.GetInt("NBinsFit", 1000);
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
    float E_fit;
    fitTree->SetBranchAddress(varFitName.c_str(), &E_fit);
    TH1D* fitHist = new TH1D("fitHist", "fitHist", NbinsFit, 0, EMaxFit);
    fitTree->Draw((varFitName + ">>fitHist").c_str(), "", "goff");
    
    TFile* fitResultFile = new TFile(fitResultFileName.c_str(), "READ");
    if (!fitResultFile || fitResultFile->IsZombie()) {
        std::cerr << "Error opening fit result file: " << fitResultFileName << std::endl;
        return -1;
    }
    TF1* fitFunc = dynamic_cast<TF1*>(fitResultFile->Get("fitFunc"));
    if (!fitFunc) {
        std::cerr << "Error retrieving fit function from file: " << fitResultFileName << std::endl;
        return -1;
    }

    TH1D* fitHistRescaled = new TH1D("fitHistRescaled", "fitHistRescaled", NbinsFit, 0, EMaxFit);
    int nFitEnt = fitTree->GetEntries();
    for (int i = 0; i < nFitEnt; ++i) {
        fitTree->GetEntry(i);
        if (E_fit < 1.0) continue; // Skip entries below 1 GeV to avoid issues with the fit function
        if (E_fit > 6.0) continue; // Skip entries above 6 GeV to avoid issues with the fit function
        fitHistRescaled->Fill(fitFunc->Eval(E_fit));
    }

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
    c2->cd();
    gStyle->SetOptStat(0);
    c2->SetLogy();
    fitHistRescaled->Scale(1.0 / fitHistRescaled->Integral("width"));
    targetHist->Scale(1.0 / targetHist->Integral("width"));
    fitHist->Scale(1.0 / fitHist->Integral("width"));
    fitHistRescaled->SetLineColor(kRed);
    fitHistRescaled->GetXaxis()->SetRangeUser(0, plotXMax);
    fitHistRescaled->GetYaxis()->SetRangeUser(1e-5, plotYMax);
    fitHistRescaled->SetLineWidth(2);
    fitHistRescaled->Draw("HIST");
    targetHist->SetLineColor(kBlue);
    targetHist->SetLineWidth(2);
    targetHist->Draw("HIST SAME");
    fitHist->SetLineColor(kGreen);
    fitHist->SetLineWidth(2);
    fitHist->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(targetHist, "Target", "l");
    legend->AddEntry(fitHist, "Fit", "l");
    legend->AddEntry(fitHistRescaled, "Fit Rescaled", "l");
    legend->Draw();

    c2->SaveAs(outputFileName.c_str());
}