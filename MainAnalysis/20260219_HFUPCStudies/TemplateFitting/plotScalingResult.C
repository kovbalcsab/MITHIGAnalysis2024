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
#include "TLegend.h"
#include <map>
#include <iostream>

#include "CommandLine.h" // Yi's Commandline bundle

using namespace RooFit;

int main(int argc, char** argv) {
    // Set up command line parser
    CommandLine CL(argc, argv);
    std::string fitFileName = CL.Get("FitFile", "fit_result.root");
    std::string plotFileName = CL.Get("PlotFileName", "scaling_fit_result.pdf");
    TFile* fitFile = new TFile(fitFileName.c_str(), "READ");
    if (!fitFile || fitFile->IsZombie()) {
        std::cerr << "Error opening fit file: " << fitFileName << std::endl;
        return -1;
    }
    fitFile->cd();

    RooFitResult* fitResult = dynamic_cast<RooFitResult*>(fitFile->Get("fitResult"));
    RooAbsPdf* fittedModel = dynamic_cast<RooAbsPdf*>(fitFile->Get("fittedModel"));
    RooDataHist* targetHist = dynamic_cast<RooDataHist*>(fitFile->Get("targetHist"));
    RooDataHist* fitHist = dynamic_cast<RooDataHist*>(fitFile->Get("fitHist"));
    if (!fitResult || !fittedModel) {
        std::cerr << "Error retrieving fit result or fitted model from file: " << fitFileName << std::endl;
        return -1;
    }

    // Print fit results
    fitResult->Print("v");

    // Create a plot of the fitted model
    TCanvas *c1 = new TCanvas("c1", "Fitted Model", 800, 800);
    RooRealVar* E_fit = dynamic_cast<RooRealVar*>(fittedModel->getVariables()->find("E_fit"));

    RooHistPdf targetHist_mapped_pdf("targetHist_mapped_pdf", "Mapped Template PDF", 
                                  RooArgList(*E_fit), 
                                  *targetHist, 2);

    RooHistPdf fitHist_pdf("fitHist_pdf", "Fit Data PDF", RooArgList(*E_fit), *fitHist, 2);
    if (E_fit) {
        RooPlot* frame = E_fit->frame();
        c1->SetLogy(1);
        fittedModel->plotOn(frame, LineColor(kBlack), Name("FittedModel"));
        targetHist_mapped_pdf.plotOn(frame, MarkerStyle(20), MarkerColor(kRed), LineColor(kRed), Name("Target"));
        fitHist_pdf.plotOn(frame, MarkerStyle(20), MarkerColor(kBlue), LineColor(kBlue), Name("FitData"));
        frame->Draw();

        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(frame->findObject("FittedModel"), "Fitted Model", "l");
        legend->AddEntry(frame->findObject("Target"), "Target Data", "p");
        legend->AddEntry(frame->findObject("FitData"), "Fit Data", "p");
        legend->Draw();
        c1->SaveAs(plotFileName.c_str());
    }

    fitFile->Close();
    delete fitFile;
    delete c1;
}