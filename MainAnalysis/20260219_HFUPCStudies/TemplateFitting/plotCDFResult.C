#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TString.h"

#include "CommandLine.h"
#include "RootIOUtils.h"

static double ComputeChi2(const TH1D *target, const TH1D *mapped, int &ndf)
{
  ndf = 0;
  if(target == nullptr || mapped == nullptr)
    return 0.0;

  const int nBins = std::min(target->GetNbinsX(), mapped->GetNbinsX());
  double chi2 = 0.0;
  for(int i = 1; i <= nBins; ++i)
  {
    const double t = target->GetBinContent(i);
    const double m = mapped->GetBinContent(i);
    const double et = target->GetBinError(i);
    const double em = mapped->GetBinError(i);
    const double variance = et * et + em * em;
    if(variance <= 0.0)
      continue;

    const double diff = t - m;
    chi2 += diff * diff / variance;
    ndf++;
  }
  return chi2;
}

int main(int argc, char *argv[])
{
  CommandLine CL(argc, argv);

  const std::string targetFileName = CL.Get("TargetFileName", CL.Get("TargetFile", ""));
  const std::string fitFileName = CL.Get("FitFileName", CL.Get("FileToFit", ""));
  const std::string fitResultFileName = CL.Get("FitResultFileName", CL.Get("FitResultFile", "output_CDF.root"));
  const std::string outputName = CL.Get("OutputFileName", "output_CDF_plot.root");

  const std::string treeName = CL.Get("TreeName", "OutputTree");
  const std::string fitVarName = CL.Get("VarFitName", CL.Get("VariableName", "HFEMaxMinus"));
  const std::string targetVarName = CL.Get("VarTargetName", fitVarName);
  const int nBins = CL.GetInt("NBins", CL.GetInt("NBinsTarget", CL.GetInt("NBinsFit", 1000)));

  const double eMinTarget = CL.GetDouble("EMinTarget", CL.GetDouble("EMin", 0.0));
  const double eMinFit = CL.GetDouble("EMinFit", CL.GetDouble("EMin", 0.0));
  const double eMaxTarget = CL.GetDouble("EMaxTarget", CL.GetDouble("EMax", 100.0));
  const double eMaxFit = CL.GetDouble("EMaxFit", CL.GetDouble("EMax", 100.0));
  const double eMin = std::min(eMinTarget, eMinFit);
  const double eMax = std::max(eMaxTarget, eMaxFit);
  const double plotYMax = CL.GetDouble("PlotYMax", 0.5);

  const double roiMin = CL.GetDouble("ROIMin", eMinFit);
  const double roiMaxInput = CL.GetDouble("ROIMax", eMaxFit);
  const double roiMax = std::max(roiMin, roiMaxInput);
  const bool useROICut = CL.GetBool("UseROICut", true);

  if(targetFileName.empty() || fitFileName.empty())
  {
    std::cerr << "Error: both TargetFile(TargetFileName) and FileToFit(FitFileName) must be provided." << std::endl;
    return 1;
  }

  gStyle->SetOptStat(0);

  std::string outputPDFName;
  std::string outputRootName;
  if(outputName.size() >= 4 && outputName.substr(outputName.size() - 4) == ".pdf")
  {
    outputPDFName = outputName;
    outputRootName = outputName.substr(0, outputName.size() - 4) + ".root";
  }
  else if(outputName.size() >= 4 && outputName.substr(outputName.size() - 4) == ".png")
  {
    outputPDFName = outputName.substr(0, outputName.size() - 4) + ".pdf";
    outputRootName = outputName.substr(0, outputName.size() - 4) + ".root";
  }
  else if(outputName.size() >= 5 && outputName.substr(outputName.size() - 5) == ".root")
  {
    outputRootName = outputName;
    outputPDFName = outputName.substr(0, outputName.size() - 5) + ".pdf";
  }
  else
  {
    outputRootName = outputName + ".root";
    outputPDFName = outputName + ".pdf";
  }

  const std::size_t slashPos = outputPDFName.find_last_of('/');
  const std::string outputDir = (slashPos == std::string::npos) ? "." : outputPDFName.substr(0, slashPos);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile *targetFile = RootIOUtils::OpenFileOrNull(targetFileName, "READ", "target input file");
  TFile *fitFile = RootIOUtils::OpenFileOrNull(fitFileName, "READ", "fit input file");
  if(targetFile == nullptr || fitFile == nullptr)
  {
    if(targetFile != nullptr)
    {
      targetFile->Close();
      delete targetFile;
    }
    if(fitFile != nullptr)
    {
      fitFile->Close();
      delete fitFile;
    }
    return 1;
  }

  TTree *targetTree = RootIOUtils::GetTreeOrNull(targetFile, treeName, "target input file: " + targetFileName);
  TTree *fitTree = RootIOUtils::GetTreeOrNull(fitFile, treeName, "fit input file: " + fitFileName);
  if(targetTree == nullptr || fitTree == nullptr)
  {
    targetFile->Close();
    fitFile->Close();
    delete targetFile;
    delete fitFile;
    return 1;
  }
  if(RootIOUtils::RequireBranchOrNull(targetTree, targetVarName, "target tree") == nullptr ||
     RootIOUtils::RequireBranchOrNull(fitTree, fitVarName, "fit tree") == nullptr)
  {
    targetFile->Close();
    fitFile->Close();
    delete targetFile;
    delete fitFile;
    return 1;
  }

  TFile *fitOutput = RootIOUtils::OpenFileOrNull(fitResultFileName, "READ", "fit result file");
  if(fitOutput == nullptr)
  {
    targetFile->Close();
    fitFile->Close();
    delete targetFile;
    delete fitFile;
    return 1;
  }

  TF1 *fitFunc = dynamic_cast<TF1 *>(fitOutput->Get("fitFunc"));
  TGraph *cdfMap = dynamic_cast<TGraph *>(fitOutput->Get("cdfMapGraph"));
  if(fitFunc == nullptr || cdfMap == nullptr)
  {
    std::cerr << "Error: fitFunc and/or cdfMapGraph missing in output_CDF.root." << std::endl;
    fitOutput->Close();
    targetFile->Close();
    fitFile->Close();
    delete fitOutput;
    delete targetFile;
    delete fitFile;
    return 1;
  }

  std::cout << "Loaded fit function: " << fitFunc->GetExpFormula("P") << std::endl;
  std::cout << "Map points: " << cdfMap->GetN() << std::endl;

  TH1D hTarget("hTarget", "Energy Distribution;E;Normalized entries", nBins, eMin, eMax);
  TH1D hFit("hFit", "Energy Distribution;E;Normalized entries", nBins, eMin, eMax);
  TH1D hMapped("hMapped", "Energy Distribution;E;Normalized entries", nBins, eMin, eMax);
  hTarget.Sumw2();
  hFit.Sumw2();
  hMapped.Sumw2();

  const std::string drawTarget = targetVarName + ">>hTarget";
  const std::string drawFit = fitVarName + ">>hFit";
  targetTree->Draw(drawTarget.c_str(), "", "goff");
  fitTree->Draw(drawFit.c_str(), "", "goff");

  float eFit = 0;
  if(fitTree->SetBranchAddress(fitVarName.c_str(), &eFit) < 0)
  {
    std::cerr << "Error: could not set branch address for " << fitVarName << std::endl;
    return 1;
  }

  const Long64_t nEvents = fitTree->GetEntries();
  Long64_t usedEvents = 0;
  for(Long64_t i = 0; i < nEvents; ++i)
  {
    fitTree->GetEntry(i);
    if(useROICut && (eFit < roiMin || eFit >= roiMax))
      continue;

    const double mappedE = fitFunc->Eval(eFit);
    hMapped.Fill(mappedE);
    ++usedEvents;
  }

  auto normalize = [](TH1D &h)
  {
    const double integral = h.Integral(1, h.GetNbinsX(), "width");
    if(integral > 0)
      h.Scale(1.0 / integral);
  };
  normalize(hTarget);
  normalize(hFit);
  normalize(hMapped);

  int chi2NDF = 0;
  const double chi2 = ComputeChi2(&hTarget, &hMapped, chi2NDF);
  const double chi2OverNDF = (chi2NDF > 0) ? (chi2 / chi2NDF) : 0.0;

  hTarget.SetLineColor(kBlue + 1);
  hTarget.SetLineWidth(3);
  hFit.SetLineColor(kRed + 1);
  hFit.SetLineWidth(2);
  hFit.SetLineStyle(2);
  hMapped.SetLineColor(kGreen + 2);
  hMapped.SetLineWidth(3);

  TCanvas c1("c1", "Energy comparison", 900, 700);
  c1.SetLogy(1);
  hTarget.SetMaximum(plotYMax);
  hTarget.Draw("hist");
  hFit.Draw("hist same");
  hMapped.Draw("hist same");

  TLegend legend(0.58, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&hTarget, "Target", "l");
  legend.AddEntry(&hFit, "FileToFit (raw)", "l");
  legend.AddEntry(&hMapped, "FileToFit (mapped)", "l");
  legend.Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.032);
  label.DrawLatex(0.13, 0.86, Form("ROI cut on fit sample: %s, [%.3g, %.3g)", useROICut ? "on" : "off", roiMin, roiMax));
  label.DrawLatex(0.13, 0.82, Form("Mapped events used: %lld / %lld", usedEvents, nEvents));
  label.DrawLatex(0.13, 0.78,
                  Form("#chi^{2}/NDF (Target vs Mapped) = %.4g / %d = %.4g", chi2, chi2NDF, chi2OverNDF));

  std::cout << "Target vs mapped chi2 = " << chi2 << ", NDF = " << chi2NDF
            << ", chi2/NDF = " << chi2OverNDF << std::endl;

  c1.SaveAs(outputPDFName.c_str());
  std::string outputPNGName = outputPDFName;
  if(outputPNGName.size() >= 4 && outputPNGName.substr(outputPNGName.size() - 4) == ".pdf")
    outputPNGName = outputPNGName.substr(0, outputPNGName.size() - 4) + ".png";
  c1.SaveAs(outputPNGName.c_str());

  TCanvas c2("c2", "CDF map and fit", 900, 700);
  cdfMap->SetMarkerStyle(20);
  cdfMap->SetMarkerSize(0.8);
  cdfMap->SetTitle("CDF mapping graph;E_{fit};E_{target}");
  cdfMap->Draw("AP");
  fitFunc->SetLineColor(kRed + 1);
  fitFunc->SetLineWidth(2);
  fitFunc->Draw("same");

  TLegend legend2(0.15, 0.74, 0.45, 0.88);
  legend2.SetBorderSize(0);
  legend2.SetFillStyle(0);
  legend2.AddEntry(cdfMap, "CDF map points", "p");
  legend2.AddEntry(fitFunc, "Polynomial fit", "l");
  legend2.Draw();

  std::string mapPDFName = outputPDFName;
  if(mapPDFName.size() >= 4 && mapPDFName.substr(mapPDFName.size() - 4) == ".pdf")
    mapPDFName = mapPDFName.substr(0, mapPDFName.size() - 4) + "_map.pdf";
  std::string mapPNGName = mapPDFName;
  if(mapPNGName.size() >= 4 && mapPNGName.substr(mapPNGName.size() - 4) == ".pdf")
    mapPNGName = mapPNGName.substr(0, mapPNGName.size() - 4) + ".png";
  c2.SaveAs(mapPDFName.c_str());
  c2.SaveAs(mapPNGName.c_str());

  TFile out(outputRootName.c_str(), "RECREATE");
  hTarget.Write();
  hFit.Write();
  hMapped.Write();
  cdfMap->Write("cdfMapGraph");
  fitFunc->Write("fitFunc");
  c1.Write("cEnergy");
  c2.Write("cMap");
  TObjString runInfo(Form("UseROICut=%d, ROIMin=%.8g, ROIMax=%.8g, FitResultFile=%s",
                          useROICut ? 1 : 0, roiMin, roiMax, fitResultFileName.c_str()));
  runInfo.Write("RunInfo");
  TObjString chi2Info(Form("Chi2=%.17g,NDF=%d,Chi2NDF=%.17g", chi2, chi2NDF, chi2OverNDF));
  chi2Info.Write("Chi2Info");
  out.Close();

  fitOutput->Close();
  targetFile->Close();
  fitFile->Close();
  delete fitOutput;
  delete targetFile;
  delete fitFile;
  return 0;
}
