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
    double aBase = CL.GetDouble("aBase", 1.0);
    double bBase = CL.GetDouble("bBase", 0.0);
    int BinsPerGeV = CL.GetDouble("BinsPerGeV", 2);

    // ---------------------------------------------------------
    // 1. Read in the TFiles and TTrees
    // ---------------------------------------------------------
    TFile* targetFile = new TFile(targetFileName.c_str(), "READ");
    if (!targetFile || targetFile->IsZombie()) {
        std::cerr << "Error opening target file: " << targetFileName << std::endl;
        return -1;
    }
    TTree* targetTree = dynamic_cast<TTree*>(targetFile->Get("OutputTree"));
    TBranch *brTarget = targetTree->GetBranch(varTargetName.c_str());
    if (brTarget) {
        brTarget->SetName("E_template");
    } else {
        std::cerr << "Could not find branch: " << varTargetName << std::endl;
    }

    TFile* fileToFit = new TFile(fileToFitName.c_str(), "READ");
    if (!fileToFit || fileToFit->IsZombie()) {
        std::cerr << "Error opening file to fit: " << fileToFitName << std::endl;
        return -1;
    }
    TTree* fitTree = dynamic_cast<TTree*>(fileToFit->Get("OutputTree"));
    TBranch *brFit = fitTree->GetBranch(varFitName.c_str());
    if (brFit) {
        brFit->SetName("E_fit");
    } else {
        std::cerr << "Could not find branch: " << varFitName << std::endl;
    }

    // ---------------------------------------------------------
    // 2. Set up RooFit Variables and Datasets
    // ---------------------------------------------------------
    // The variable for the fit dataset MUST match what the final PDF uses
    RooRealVar E_fit("E_fit", "E_fit", 0, 100); 
    E_fit.setBins(100 * BinsPerGeV);
    RooDataSet* fitData = new RooDataSet("fitData", "fitData", RooArgSet(E_fit), Import(*fitTree));

    // A temporary variable representing the static template shape from the target tree
    RooRealVar E_template_var("E_template", "E_template", 0, 150);
    E_template_var.setBins(150 * BinsPerGeV);
    RooDataSet* targetData = new RooDataSet("targetData", "targetData", RooArgSet(E_template_var), Import(*targetTree));
    TH1D* targetHistTemp = new TH1D("targetHistTemp", "targetHistTemp", 800, 0, 150);
    targetData->fillHistogram(targetHistTemp, RooArgList(E_template_var));
    for (int i = 1; i <= targetHistTemp->GetNbinsX(); ++i) {
        double content = targetHistTemp->GetBinContent(i);
        if (content == 0) {
            targetHistTemp->SetBinContent(i, 1); // Avoid zero bins
        }
    }

    // ---------------------------------------------------------
    // 3. Create the Base Template PDF
    // ---------------------------------------------------------
    RooDataHist targetHist("targetHist", "targetHist", RooArgSet(E_template_var), targetHistTemp);
    RooHistPdf targetPdf("targetPdf", "targetPdf", RooArgSet(E_template_var), targetHist, 1);
    RooDataHist fitHist("fitHist", "fitHist", RooArgSet(E_fit), *fitData);
    E_fit.setRange("fitRange", 1, 25);

    // ---------------------------------------------------------
    // 4. Define the Transformation and Apply it via Customizer
    // ---------------------------------------------------------
    RooRealVar a("a", "a", aBase, 0., 5); // Scale
    RooRealVar b("b", "b", bBase, -10, 10); // Non-linear offset
    
    // scaledE: what value of the 'target' corresponds to the current 'fit' value
    RooFormulaVar scaledE("scaledE", "@1 * @0 + @2 * @0 * @0", RooArgSet(E_fit, a, b));

    // Swap the static template variable for the functional formula
    RooCustomizer cust(targetPdf, "scaled");
    cust.replaceArg(E_template_var, scaledE);
    RooAbsPdf* Basemodel = dynamic_cast<RooAbsPdf*>(cust.build());

    // ---------------------------------------------------------
    // 5. Apply the Jacobian and Build Final Model
    // ---------------------------------------------------------
    // Jacobian: d(scaledE)/d(E_fit) = a + 2*b*E_fit
    RooFormulaVar jac("jac", "abs(@1 + 2*@2*@0)", RooArgList(E_fit, a, b));
    
    // Combine base model with Jacobian
    RooEffProd model("model", "model", *Basemodel, jac);

    // ---------------------------------------------------------
    // 6. Execute the Fit
    // ---------------------------------------------------------
    RooFitResult* fitResult = model.fitTo(*fitData, 
                                          Save(true), 
                                          Range("fitRange"),
                                          Minimizer("Minuit2"), 
                                          Strategy(2), 
                                          PrintLevel(-1));

    if (fitResult) fitResult->Print("v");

    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        return -1;
    }
    outputFile->cd();
    fitResult->Write("fitResult");
    model.Write("fittedModel");
    targetHist.Write("targetHist");
    fitHist.Write("fitHist");

    // ---------------------------------------------------------
    // 7. Visual Check & Output
    // ---------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "Fit", 800, 600);
    RooPlot* frame = E_fit.frame();
    fitData->plotOn(frame);
    model.plotOn(frame);
    frame->Draw();
    c1->Write("fitPlot");

    outputFile->Close();

    // Clean up (optional but good practice)
    delete targetData;
    delete fitData;
    delete targetFile;
    delete fileToFit;
    delete outputFile;

    return 0;
}