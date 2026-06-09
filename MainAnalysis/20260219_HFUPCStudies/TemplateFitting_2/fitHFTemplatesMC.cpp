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
  RooRealVar variable(variableName.c_str(), variableName.c_str(), 0, 200);
  TTree *inTree = RootIOUtils::GetTreeOrNull(input, "OutputTree", "input file");
  RooDataSet dataset("dataset", "dataset", inTree, RooArgSet(variable));

  // Load input templates
  int numTemplates = templatePaths.size();
  std::vector<TH1D *> template_hists;
  std::vector<RooDataHist> template_datahists;
  std::vector<RooHistPdf> template_pdfs;
  RooArgList templatePdfList;
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    TFile *templateFile = RootIOUtils::OpenFileOrNull(templatePaths[iTemplate], "READ", "template file");
    TTree *templateTree = RootIOUtils::GetTreeOrNull(templateFile, "OutputTree", "template file");
    // TODO: Add custom binning option for templates
    TH1D *hist =
        new TH1D(Form("template_%i", iTemplate), Form("template_%s", templateNames[iTemplate].c_str()), 100, 0, 200);
    templateTree->Draw(Form("%s>>template_%i", variableName.c_str(), iTemplate), "", "goff");
    if (hist) {
      template_hists.push_back(hist);
      RooDataHist dataHist(Form("templateDist_%i", iTemplate),
                           Form("templateDist_%s", templateNames[iTemplate].c_str()), RooArgSet(variable), hist);
      template_datahists.push_back(dataHist);
      RooHistPdf histPdf(Form("templatePdf_%i", iTemplate), Form("templatePdf_%s", templateNames[iTemplate].c_str()),
                         RooArgSet(variable), dataHist);
      template_pdfs.push_back(histPdf);
      templatePdfList.add(histPdf);
    }
    templateFile->Close();
  }

  std::vector<RooRealVar> template_fractions;
  RooArgList fractionList;
  for (int iTemplate = 0; iTemplate < numTemplates - 1; ++iTemplate) {
    // TODO: Add default fraction as input for fitting
    RooRealVar fraction(Form("fraction_%i", iTemplate), Form("fraction_%s", templateNames[iTemplate].c_str()), 0.5, 0,
                        1);
    template_fractions.push_back(fraction);
    fractionList.add(fraction);
  }

  RooAddPdf model("model", "model", templatePdfList, fractionList);
  RooFitResult *fitResult = model.fitTo(dataset, RooFit::Save(true));
  fitResult->Print("v");
  double fracLast = 1.0 - std::accumulate(template_fractions.begin(), template_fractions.end(), 0.0,
                                          [](double sum, const RooRealVar &frac) { return sum + frac.getVal(); });
  for (size_t iTemplate = 0; iTemplate < template_fractions.size() - 1; ++iTemplate) {
    std::cout << "Fraction of " << templateNames[iTemplate] << ": " << template_fractions[iTemplate].getVal()
              << std::endl;
  }
  std::cout << "Fraction of " << templateNames.back() << ": " << fracLast << std::endl;

  // Plotting
  RooPlot *frame = variable.frame();
  dataset.plotOn(frame);
  model.plotOn(frame);
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    model.plotOn(frame, RooFit::Components(Form("templatePdf_%i", iTemplate)), RooFit::LineStyle(kDashed),
                 RooFit::Name(templateNames[iTemplate].c_str()));
  }

  TCanvas *canvas = new TCanvas("canvas", "Template Fit", 800, 800);
  frame->Draw();
  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(frame->getObject(0), "Data", "ep");
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    legend->AddEntry(frame->getObject(iTemplate + 1), templateNames[iTemplate].c_str(), "l");
  }
  std::string strippedPath = outputFile.substr(0, outputFile.find_last_of("."));
  canvas->SaveAs((strippedPath + "_fit.png").c_str());

  output->cd();
  infoManager.SaveToFile();
  frame->Write("fitFrame");
  model.Write("fitModel");
  for (int iTemplate = 0; iTemplate < numTemplates; ++iTemplate) {
    template_hists[iTemplate]->Write(Form("templateHist_%s", templateNames[iTemplate].c_str()));
    template_pdfs[iTemplate].Write(Form("templatePdf_%s", templateNames[iTemplate].c_str()));
  }
  input->Close();

  return 0;
}
