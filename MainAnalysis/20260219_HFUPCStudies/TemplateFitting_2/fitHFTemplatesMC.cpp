#include "CommandLine.h"
#include "InfoManager.h"
#include "RooGlobalFunc.h"
#include "RootIOUtils.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"

#include <vector>

const std::vector<double>
    defaultBinEdges = {0,   0.25, 0.5, 0.75, 1,   1.25, 1.5, 1.75, 2,  2.25, 2.5, 2.75, 3, 3.25,
                       3.5, 3.75, 4,   4.25, 4.5, 4.75, 5,   5.5,  6,  6.5,  7,   7.5,  8, 9,
                       10,  11,   12,  13,   14,  16,   18,  20,   25, 30,   40}; // 60,   80,   100,  120,  140,  160,
                                                                                  // 180,  200};

int main(int argc, char *argv[]) {
  CommandLine cmdLine(argc, argv);
  std::string inputFile = cmdLine.Get("inputFile", "input.root");
  bool isPlus = cmdLine.GetBool("isPlus", false);
  std::string outputFile = cmdLine.Get("outputFile", "output.root");
  std::vector<std::string> templatePaths = cmdLine.GetStringVector("templatePaths", std::vector<std::string>());
  std::vector<std::string> templateNames = cmdLine.GetStringVector("templateNames", std::vector<std::string>());

  if (templatePaths.size() != templateNames.size()) {
    std::cerr << "Error: The number of template paths must match the number of template names." << std::endl;
    return 1;
  }

  TFile *input = RootIOUtils::OpenFileOrNull(inputFile, "READ", "input file");

  TFile *output = TFile::Open(outputFile.c_str(), "RECREATE");

  GeneralInfoManager infoManager(output, "InfoDir", false);
  TTimeStamp *currentTime = new TTimeStamp();
  infoManager.AddSourceFile(inputFile, currentTime);
  infoManager.AddCutParameter("isPlus", isPlus, currentTime);

  // Load dataset
  std::string variableName = isPlus ? "HFEMaxPlus_forest" : "HFEMaxMinus_forest";
  RooRealVar variable(variableName.c_str(), variableName.c_str(), 0, defaultBinEdges.back());
  TTree *inTree = RootIOUtils::GetTreeOrNull(input, "OutputTree", "input file");
  RooDataSet dataset("dataset", "dataset", inTree, RooArgSet(variable));

  // Load input templates
  int numTemplates = templatePaths.size();
  std::vector<TH1D *> template_hists;
  std::vector<RooDataHist *> template_datahists;
  std::vector<RooHistPdf *> template_pdfs;
  RooArgList templatePdfList;
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    TFile *templateFile = RootIOUtils::OpenFileOrNull(templatePaths[iTemplate], "READ", "template file");
    TTree *templateTree = RootIOUtils::GetTreeOrNull(templateFile, "OutputTree", "template file");
    TH1D *hist = new TH1D(Form("template_%i", iTemplate), Form("template_%s", templateNames[iTemplate].c_str()),
                          defaultBinEdges.size() - 1, defaultBinEdges.data());
    hist->Sumw2();
    templateTree->Draw(Form("%s>>template_%i", variableName.c_str(), iTemplate), "", "goff");
    hist->SetDirectory(0); // Detach histogram from file to prevent it from being deleted when file is closed
    if (hist) {
      template_hists.push_back(hist);
      RooDataHist *dataHist =
          new RooDataHist(Form("templateDist_%i", iTemplate), Form("templateDist_%s", templateNames[iTemplate].c_str()),
                          RooArgSet(variable), hist);
      template_datahists.push_back(dataHist);
      RooHistPdf *histPdf =
          new RooHistPdf(Form("templatePdf_%i", iTemplate), Form("templatePdf_%s", templateNames[iTemplate].c_str()),
                         RooArgSet(variable), *dataHist);
      template_pdfs.push_back(histPdf);
      templatePdfList.add(*histPdf);
    }
    templateFile->Close();
  }

  std::vector<RooRealVar *> template_fractions;
  RooArgList fractionList;
  for (int iTemplate = 0; iTemplate < numTemplates - 1; ++iTemplate) {
    // TODO: Add default fraction as input for fitting
    RooRealVar *fraction = new RooRealVar(Form("fraction_%i", iTemplate),
                                          Form("fraction_%s", templateNames[iTemplate].c_str()), 0.5, 0, 1);
    template_fractions.push_back(fraction);
    fractionList.add(*fraction);
  }

  RooAddPdf model("model", "model", templatePdfList, fractionList);
  RooFitResult *fitResult = model.fitTo(dataset, RooFit::Save(true));
  fitResult->Print("v");
  double fracLast = 1.0 - std::accumulate((template_fractions).begin(), (template_fractions).end(), 0.0,
                                          [](double sum, RooRealVar *frac) { return sum + frac->getVal(); });
  for (size_t iTemplate = 0; iTemplate < template_fractions.size(); iTemplate++) {
    std::cout << "Fraction of " << templateNames[iTemplate] << ": " << template_fractions[iTemplate]->getVal()
              << std::endl;
  }
  std::cout << "Fraction of " << templateNames.back() << ": " << fracLast << std::endl;

  // Plotting
  RooPlot *frame = variable.frame();
  dataset.plotOn(frame);
  model.plotOn(frame);
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    model.plotOn(frame, RooFit::Components(template_pdfs[iTemplate]->GetName()), RooFit::LineStyle(kDashed),
                 RooFit::Name(templateNames[iTemplate].c_str()));
  }

  TCanvas *canvas = new TCanvas("canvas", "Template Fit", 800, 800);
  frame->Draw();
  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(frame->getObject(0), "Data", "ep");
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    legend->AddEntry(frame->getObject(iTemplate + 1), templateNames[iTemplate].c_str(), "l");
  }
  legend->Draw();
  std::string strippedPath = outputFile.substr(0, outputFile.find_last_of("."));
  // canvas->SaveAs((strippedPath + "_fit.png").c_str());

  output->cd();
  frame->Write("fitFrame");
  model.Write("fitModel");
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    template_hists[iTemplate]->Write();
    template_pdfs[iTemplate]->Write();
  }
  infoManager.SaveToFile();
  output->Close();
  input->Close();

  return 0;
}
