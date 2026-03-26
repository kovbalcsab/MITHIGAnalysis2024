#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "RooCustomizer.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
#include "RooKeysPdf.h"
#include "RooNumIntConfig.h"
#include "RooRealVar.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#include "CommandLine.h" // Yi's Commandline bundle
#include "InfoManager.h"
#include "RootIOUtils.h"

using namespace RooFit;

static std::string StripExtension(const std::string &path) {
    const size_t slashPos = path.find_last_of('/');
    const size_t dotPos = path.find_last_of('.');
    if (dotPos == std::string::npos)
        return path;
    if (slashPos != std::string::npos && dotPos < slashPos)
        return path;
    return path.substr(0, dotPos);
}

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
    int BinsPerGeV = CL.GetInt("BinsPerGeV", 2);
    if (BinsPerGeV <= 0) {
        std::cerr << "Invalid BinsPerGeV: must be > 0." << std::endl;
        return -1;
    }

    int fitNumCPU = CL.GetInt("FitNumCPU", 4);
    fitNumCPU = std::max(1, fitNumCPU);
    int fitStrategy = CL.GetInt("FitStrategy", 1);
    int fitPrintLevel = CL.GetInt("FitPrintLevel", -1);
    bool fitOffset = CL.GetBool("FitOffset", true);
    bool clipScaledE = CL.GetBool("ClipScaledE", true);
    bool useAbsJacobian = CL.GetBool("UseAbsJacobian", true);
    double aMin = CL.GetDouble("AMin", 0.0);
    double aMax = CL.GetDouble("AMax", 5.0);
    double bMin = CL.GetDouble("BMin", -10.0);
    double bMax = CL.GetDouble("BMax", 10.0);
    if (!(aMin < aMax)) {
        std::cerr << "Invalid a-range: AMin must be < AMax." << std::endl;
        return -1;
    }
    if (!(bMin < bMax)) {
        std::cerr << "Invalid b-range: BMin must be < BMax." << std::endl;
        return -1;
    }
    if (aBase < aMin || aBase > aMax) {
        std::cerr << "Warning: aBase is outside [AMin, AMax], clamping to valid range." << std::endl;
        aBase = std::clamp(aBase, aMin, aMax);
    }
    if (bBase < bMin || bBase > bMax) {
        std::cerr << "Warning: bBase is outside [BMin, BMax], clamping to valid range." << std::endl;
        bBase = std::clamp(bBase, bMin, bMax);
    }

    // ---------------------------------------------------------
    // 1. Read in the TFiles and TTrees
    // ---------------------------------------------------------
    TFile* targetFile = RootIOUtils::OpenFileOrNull(targetFileName, "READ", "target file");
    if (targetFile == nullptr) {
        return -1;
    }
    TTree* targetTree = RootIOUtils::GetTreeOrNull(targetFile, "OutputTree", "target file: " + targetFileName);
    if (targetTree == nullptr) {
        RootIOUtils::CloseAndDeleteFile(targetFile);
        return -1;
    }
    TBranch *brTarget = RootIOUtils::RequireBranchOrNull(targetTree, varTargetName, "target tree");
    if (brTarget != nullptr) {
        brTarget->SetName("E_template");
    } else {
        RootIOUtils::CloseAndDeleteFile(targetFile);
        return -1;
    }

    TFile* fileToFit = RootIOUtils::OpenFileOrNull(fileToFitName, "READ", "fit file");
    if (fileToFit == nullptr) {
        RootIOUtils::CloseAndDeleteFile(targetFile);
        return -1;
    }
    TTree* fitTree = RootIOUtils::GetTreeOrNull(fileToFit, "OutputTree", "fit file: " + fileToFitName);
    if (fitTree == nullptr) {
        RootIOUtils::CloseAndDeleteFile(targetFile);
        RootIOUtils::CloseAndDeleteFile(fileToFit);
        return -1;
    }
    TBranch *brFit = RootIOUtils::RequireBranchOrNull(fitTree, varFitName, "fit tree");
    if (brFit != nullptr) {
        brFit->SetName("E_fit");
    } else {
        RootIOUtils::CloseAndDeleteFile(targetFile);
        RootIOUtils::CloseAndDeleteFile(fileToFit);
        return -1;
    }

    // ---------------------------------------------------------
    // 2. Set up RooFit Variables and Datasets
    // ---------------------------------------------------------
    // The variable for the fit dataset MUST match what the final PDF uses
    RooRealVar E_fit("E_fit", "E_fit", 0, 100); 
    E_fit.setBins(100 * BinsPerGeV);
    RooDataSet* fitData = new RooDataSet("fitData", "fitData", RooArgSet(E_fit), Import(*fitTree));
    if (fitData == nullptr || fitData->numEntries() <= 0) {
        std::cerr << "Fit dataset is empty or invalid for tree branch '" << varFitName << "'." << std::endl;
        return -1;
    }

    // A temporary variable representing the static template shape from the target tree
    RooRealVar E_template_var("E_template", "E_template", 0, 150);
    E_template_var.setBins(150 * BinsPerGeV);
    RooDataSet* targetData = new RooDataSet("targetData", "targetData", RooArgSet(E_template_var), Import(*targetTree));
    if (targetData == nullptr || targetData->numEntries() <= 0) {
        std::cerr << "Target dataset is empty or invalid for tree branch '" << varTargetName << "'." << std::endl;
        return -1;
    }
    RooDataSet* targetDataReduced = nullptr;
    int maxKDEEntries = CL.GetInt("MaxKDEEntries", 100000); // Limit for KDE to speed up fit
    maxKDEEntries = std::min(maxKDEEntries, (int)targetData->numEntries());
    maxKDEEntries = std::max(maxKDEEntries, 1); // Ensure it's at least 1 to avoid issues

    targetDataReduced = dynamic_cast<RooDataSet*>(
        targetData->reduce(SelectVars(RooArgSet(E_template_var)),
                        EventRange(0, maxKDEEntries))
    );
    if (targetDataReduced == nullptr || targetDataReduced->numEntries() <= 0) {
        std::cerr << "Reduced target dataset is empty or invalid." << std::endl;
        return -1;
    }
    TH1D* targetHistTemp = new TH1D("targetHistTemp", "targetHistTemp", 150 * BinsPerGeV, 0, 150);
    targetData->fillHistogram(targetHistTemp, RooArgList(E_template_var));

    // ---------------------------------------------------------
    // 3. Create the Base Template PDF
    // ---------------------------------------------------------
    std::cout << "Creating base template PDF..." << std::endl;
    RooDataHist targetHist("targetHist", "targetHist", RooArgSet(E_template_var), targetHistTemp);
    std::unique_ptr<RooAbsPdf> targetPdf;
    std::cout << "Using RooKeysPdf template (high quality, slower)." << std::endl;
    targetPdf = std::make_unique<RooKeysPdf>("targetPdf", "targetPdf",
                    E_template_var,
                    *targetDataReduced,
                    RooKeysPdf::MirrorBoth,
                    2.0);
  
    RooDataHist fitHist("fitHist", "fitHist", RooArgSet(E_fit), *fitData);
    E_fit.setRange("fitRange", 0, 25);

    // ---------------------------------------------------------
    // 4. Define the Transformation and Apply it via Customizer
    // ---------------------------------------------------------
    RooRealVar a("a", "a", aBase, aMin, aMax); // Scale
    RooRealVar b("b", "b", bBase, bMin, bMax); // Non-linear offset
    
    // scaledE: what value of the 'target' corresponds to the current 'fit' value
    std::string scaledEExpr = "@1*@0 + @2*@0*@0";
    if (clipScaledE)
        scaledEExpr = "TMath::Max(0.0, TMath::Min(150.0, @1*@0 + @2*@0*@0))";
    RooFormulaVar scaledE("scaledE", scaledEExpr.c_str(), RooArgList(E_fit, a, b));

    // Swap the static template variable for the functional formula
    RooCustomizer cust(*targetPdf, "scaled");
    cust.replaceArg(E_template_var, scaledE);
    RooAbsPdf* Basemodel = dynamic_cast<RooAbsPdf*>(cust.build());
    if (Basemodel == nullptr) {
        std::cerr << "Failed to build customized base model." << std::endl;
        return -1;
    }

    // ---------------------------------------------------------
    // 5. Apply the Jacobian and Build Final Model
    // ---------------------------------------------------------
    // Jacobian: d(scaledE)/d(E_fit) = a + 2*b*E_fit
    std::string jacExpr = "(@1 + 2*@2*@0)";
    if (useAbsJacobian)
        jacExpr = "abs(@1 + 2*@2*@0)";
    RooFormulaVar jac("jac", jacExpr.c_str(), RooArgList(E_fit, a, b));
    
    // Combine base model with Jacobian
    RooEffProd model("model", "model", *Basemodel, jac);

    // ---------------------------------------------------------
    // 6. Execute the Fit
    // ---------------------------------------------------------
    std::cout << "Starting fit..." << std::endl;
    RooFitResult* fitResult = model.fitTo(*fitData, 
                                          Save(true), 
                                          Range("fitRange"),
                                          NormRange("fitRange"),
                                          Minimizer("Minuit2"), 
                                          Strategy(fitStrategy),
                                          Offset(fitOffset),
                                          NumCPU(fitNumCPU, RooFit::Hybrid),
                                          Optimize(true),
                                          PrintLevel(fitPrintLevel));

    if (fitResult == nullptr) {
        std::cerr << "Fit failed and returned a null RooFitResult." << std::endl;
        return -1;
    }
    fitResult->Print("v");

    TFile* outputFile = RootIOUtils::OpenFileOrNull(outputFileName, "RECREATE", "output file");
    if (outputFile == nullptr) {
        RootIOUtils::CloseAndDeleteFile(targetFile);
        RootIOUtils::CloseAndDeleteFile(fileToFit);
        return -1;
    }
    TTimeStamp *currentTime = new TTimeStamp();
    GeneralInfoManager man(outputFile, "InfoDir", false);
    man.AddSourceFile(targetFileName, currentTime);
    man.AddSourceFile(fileToFitName, currentTime);
    man.AddCutParameter("BinsPerGeV", BinsPerGeV, currentTime);
    man.AddCutParameter("FitNumCPU", fitNumCPU, currentTime);
    man.AddCutParameter("FitStrategy", fitStrategy, currentTime);
    man.AddCutParameter("ClipScaledE", clipScaledE ? 1.0 : 0.0, currentTime);
    man.AddCutParameter("UseAbsJacobian", useAbsJacobian ? 1.0 : 0.0, currentTime);
    man.AddCutParameter("aBase", aBase, currentTime);
    man.AddCutParameter("bBase", bBase, currentTime);
    outputFile->cd();
    fitResult->Write("fitResult");
    model.Write("fittedModel");
    targetHist.Write("targetHist");
    fitHist.Write("fitHist");
    // Print the fit parameters to a text file alongside the ROOT output
    const std::string paramFileName = StripExtension(outputFileName) + "_fit_parameters.txt";
    std::ofstream paramFile(paramFileName);
    if (paramFile.is_open()) {
        paramFile << "Fit Parameters:\n";
        paramFile << "a: " << a.getVal() << " ± " << a.getError() << "\n";
        paramFile << "b: " << b.getVal() << " ± " << b.getError() << "\n";
        paramFile.close();
        std::cout << "Saved fit parameters to: " << paramFileName << std::endl;
    } else {
        std::cerr << "Unable to open file to write fit parameters: " << paramFileName << std::endl;
    }

    // Clean up (optional but good practice)
    delete targetData;
    delete targetDataReduced;
    delete targetHistTemp;
    delete fitData;
    RootIOUtils::CloseAndDeleteFile(targetFile);
    RootIOUtils::CloseAndDeleteFile(fileToFit);
    RootIOUtils::CloseAndDeleteFile(outputFile);
    delete currentTime;

    return 0;
}
