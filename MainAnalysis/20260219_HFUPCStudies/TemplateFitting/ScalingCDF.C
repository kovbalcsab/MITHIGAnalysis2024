#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "InfoManager.h"
#include "RootIOUtils.h"

namespace
{
struct TFileCloser
{
  void operator()(TFile *file) const
  {
    if(file == nullptr)
      return;
    if(file->IsOpen())
      file->Close();
    delete file;
  }
};

class TH1AddDirectoryGuard
{
public:
  explicit TH1AddDirectoryGuard(bool addDirectory)
    : PreviousState(TH1::AddDirectoryStatus())
  {
    TH1::AddDirectory(addDirectory);
  }

  ~TH1AddDirectoryGuard()
  {
    TH1::AddDirectory(PreviousState);
  }

private:
  bool PreviousState;
};

std::string BuildPolynomialFormula(int order)
{
  std::string formula = "[0]";
  for(int i = 1; i <= order; ++i)
    formula += " + [" + std::to_string(i) + "]*pow(x," + std::to_string(i) + ")";
  return formula;
}

std::string StripExtension(const std::string &path)
{
  const std::size_t slashPos = path.find_last_of("/\\");
  const std::size_t dotPos = path.find_last_of('.');
  if(dotPos == std::string::npos || (slashPos != std::string::npos && dotPos < slashPos))
    return path;
  return path.substr(0, dotPos);
}

std::string DirectoryFromPath(const std::string &path)
{
  const std::size_t slashPos = path.find_last_of("/\\");
  if(slashPos == std::string::npos)
    return ".";
  if(slashPos == 0)
    return "/";
  return path.substr(0, slashPos);
}

TH1D *BuildROIHistogram(TTree *tree, const std::string &varName, const std::string &histName,
                        int nBins, double roiMin, double roiMax)
{
  if(tree == nullptr)
    return nullptr;

  TH1D *h = new TH1D(histName.c_str(), (histName + ";" + varName + ";Entries").c_str(), nBins, roiMin, roiMax);
  h->Sumw2();

  std::ostringstream cutExpr;
  cutExpr << varName << ">=" << roiMin << " && " << varName << "<" << roiMax;
  tree->Draw((varName + ">>" + histName).c_str(), cutExpr.str().c_str(), "goff");

  // Explicitly suppress underflow/overflow content for ROI-only workflow.
  h->SetBinContent(0, 0.0);
  h->SetBinContent(h->GetNbinsX() + 1, 0.0);
  return h;
}
} // namespace

int main(int argc, char **argv)
{
  CommandLine CL(argc, argv);

  const std::string targetFileName = CL.Get("TargetFile", "target.root");
  const std::string fileToFitName = CL.Get("FileToFit", "fit.root");
  const std::string outputFileName = CL.Get("OutputFileName", "ScalingCDFResult.root");
  const std::string varFitName = CL.Get("VarFitName", "HFEMaxPlus_forest");
  const std::string varTargetName = CL.Get("VarTargetName", "HFEMaxPlus_forest");
  const std::string treeName = CL.Get("TreeName", "OutputTree");

  const int nBinsTarget = CL.GetInt("NBinsTarget", 1000);
  const int nBinsFit = CL.GetInt("NBinsFit", 1000);
  const double eMinTarget = CL.GetDouble("EMinTarget", 0.0);
  const double eMaxTargetROI = CL.GetDouble("EMaxTargetROI", 250.0);
  const double eMinFit = CL.GetDouble("EMinFit", 0.0);
  const double eMaxFitROI = CL.GetDouble("EMaxFitROI", 250.0);

  const int polyOrder = CL.GetInt("PolyOrder", 3);
  const bool forceZeroConstant = CL.GetBool("ForceZeroConstant", false);
  const double mapFitMin = CL.GetDouble("MapFitMin", std::max(eMinFit, eMinTarget));
  const double mapFitMax = CL.GetDouble("MapFitMax", std::min(eMaxFitROI, eMaxTargetROI));
  const double plotXMax = CL.GetDouble("PlotXMax", 25.0);
  const double plotYMax = CL.GetDouble("PlotYMax", 25.0);

  if(nBinsTarget <= 0 || nBinsFit <= 0)
  {
    std::cerr << "NBinsTarget and NBinsFit must be > 0." << std::endl;
    return -1;
  }
  if(!(eMaxTargetROI > eMinTarget) || !(eMaxFitROI > eMinFit))
  {
    std::cerr << "Invalid ROI limits: require EMaxTargetROI>EMinTarget and EMaxFitROI>EMinFit." << std::endl;
    return -1;
  }
  if(polyOrder < 1)
  {
    std::cerr << "PolyOrder must be >= 1." << std::endl;
    return -1;
  }
  if(!(mapFitMax > mapFitMin))
  {
    std::cerr << "Invalid map fit window: MapFitMax must be > MapFitMin." << std::endl;
    return -1;
  }


  std::unique_ptr<TFile, TFileCloser> targetFile(
    RootIOUtils::OpenFileOrNull(targetFileName, "READ", "target file"));
  std::unique_ptr<TFile, TFileCloser> fitFile(
    RootIOUtils::OpenFileOrNull(fileToFitName, "READ", "fit file"));
  if(targetFile == nullptr || fitFile == nullptr)
    return -1;

  TTree *targetTree = RootIOUtils::GetTreeOrNull(targetFile.get(), treeName, "target file: " + targetFileName);
  TTree *fitTree = RootIOUtils::GetTreeOrNull(fitFile.get(), treeName, "fit file: " + fileToFitName);
  if(targetTree == nullptr || fitTree == nullptr)
    return -1;

  if(RootIOUtils::RequireBranchOrNull(targetTree, varTargetName, "target tree") == nullptr ||
     RootIOUtils::RequireBranchOrNull(fitTree, varFitName, "fit tree") == nullptr)
    return -1;

  std::unique_ptr<TH1D> targetHistROI(
    BuildROIHistogram(targetTree, varTargetName, "targetHistROI", nBinsTarget, eMinTarget, eMaxTargetROI));
  std::unique_ptr<TH1D> fitHistROI(
    BuildROIHistogram(fitTree, varFitName, "fitHistROI", nBinsFit, eMinFit, eMaxFitROI));
  if(targetHistROI == nullptr || fitHistROI == nullptr)
  {
    std::cerr << "Failed to build ROI histograms." << std::endl;
    return -1;
  }

  const double targetIntegralROI = targetHistROI->Integral(1, targetHistROI->GetNbinsX());
  const double fitIntegralROI = fitHistROI->Integral(1, fitHistROI->GetNbinsX());
  if(targetIntegralROI <= 0.0 || fitIntegralROI <= 0.0)
  {
    std::cerr << "ROI histogram integral is non-positive. Check ROI and inputs." << std::endl;
    return -1;
  }
  targetHistROI->Scale(1.0 / targetIntegralROI);
  fitHistROI->Scale(1.0 / fitIntegralROI);

  std::unique_ptr<TH1D> cdfFit(dynamic_cast<TH1D *>(fitHistROI->GetCumulative()));
  if(cdfFit == nullptr)
  {
    std::cerr << "Failed to build cumulative histogram for fit ROI." << std::endl;
    return -1;
  }
  cdfFit->SetName("cdfFitROI");

  std::unique_ptr<TGraph> cdfMapGraph(new TGraph());
  cdfMapGraph->SetName("cdfMapGraph");
  cdfMapGraph->SetTitle("CDF mapping graph;E_{fit} [GeV];E_{target} [GeV]");

  std::vector<double> probs(1);
  std::vector<double> quantiles(1);
  for(int i = 1; i <= cdfFit->GetNbinsX(); ++i)
  {
    const double eFit = cdfFit->GetBinCenter(i);
    if(eFit < mapFitMin || eFit > mapFitMax)
      continue;

    probs[0] = std::clamp(cdfFit->GetBinContent(i), 0.0, 1.0);
    targetHistROI->GetQuantiles(1, quantiles.data(), probs.data());
    cdfMapGraph->SetPoint(cdfMapGraph->GetN(), eFit, quantiles[0]);
  }

  const int nFreeParameters = forceZeroConstant ? polyOrder : (polyOrder + 1);
  if(cdfMapGraph->GetN() < nFreeParameters)
  {
    std::cerr << "Not enough CDF mapping points (" << cdfMapGraph->GetN()
              << ") for polynomial order " << polyOrder
              << " with ForceZeroConstant=" << (forceZeroConstant ? "true" : "false")
              << "." << std::endl;
    return -1;
  }

  std::unique_ptr<TF1> fitFunc(new TF1("fitFunc", BuildPolynomialFormula(polyOrder).c_str(), mapFitMin, mapFitMax));
  for(int i = 0; i <= polyOrder; ++i)
    fitFunc->SetParameter(i, 0.0);
  fitFunc->SetParameter(1, 1.0);
  if(forceZeroConstant)
    fitFunc->FixParameter(0, 0.0);
  cdfMapGraph->Fit(fitFunc.get(), "QR");

  const std::string outputDir = DirectoryFromPath(outputFileName);
  if(!outputDir.empty() && outputDir != ".")
    gSystem->mkdir(outputDir.c_str(), true);

  std::unique_ptr<TFile, TFileCloser> outputFile(
    RootIOUtils::OpenFileOrNull(outputFileName, "RECREATE", "output file"));
  if(outputFile == nullptr)
    return -1;

  TTimeStamp currentTime;
  GeneralInfoManager man(outputFile.get(), "InfoDir", false);
  man.AddSourceFile(targetFileName, &currentTime);
  man.AddSourceFile(fileToFitName, &currentTime);
  man.AddCutParameter("NBinsTarget", nBinsTarget, &currentTime);
  man.AddCutParameter("NBinsFit", nBinsFit, &currentTime);
  man.AddCutParameter("EMinTarget", eMinTarget, &currentTime);
  man.AddCutParameter("EMaxTargetROI", eMaxTargetROI, &currentTime);
  man.AddCutParameter("EMinFit", eMinFit, &currentTime);
  man.AddCutParameter("EMaxFitROI", eMaxFitROI, &currentTime);
  man.AddCutParameter("MapFitMin", mapFitMin, &currentTime);
  man.AddCutParameter("MapFitMax", mapFitMax, &currentTime);
  man.AddCutParameter("PolyOrder", polyOrder, &currentTime);
  man.AddCutParameter("ForceZeroConstant", forceZeroConstant ? 1.0 : 0.0, &currentTime);
  man.AddCutParameter("PlotXMax", plotXMax, &currentTime);
  man.AddCutParameter("PlotYMax", plotYMax, &currentTime);

  outputFile->cd();
  targetHistROI->Write("targetHistROI");
  fitHistROI->Write("fitHistROI");
  cdfFit->Write("cdfFitROI");
  cdfMapGraph->Write("cdfMapGraph");
  fitFunc->Write("fitFunc");
  man.SaveToFile();

  std::unique_ptr<TCanvas> canvas(new TCanvas("cCDFMap", "cCDFMap", 800, 800));
  canvas->cd();
  gStyle->SetOptStat(0);
  cdfMapGraph->SetMarkerStyle(20);
  cdfMapGraph->SetMarkerColor(kBlue + 1);
  cdfMapGraph->SetLineColor(kBlue + 1);
  cdfMapGraph->SetLineWidth(2);
  cdfMapGraph->GetXaxis()->SetRangeUser(eMinFit, plotXMax);
  cdfMapGraph->GetYaxis()->SetRangeUser(eMinTarget, plotYMax);
  cdfMapGraph->SetTitle("ROI CDF mapping;E_{fit} [GeV];E_{target} [GeV]");
  cdfMapGraph->Draw("ALP");
  fitFunc->SetLineColor(kRed + 1);
  fitFunc->SetLineWidth(2);
  fitFunc->Draw("SAME");

  canvas->SaveAs((StripExtension(outputFileName) + "_CDF_fit.pdf").c_str());
  return 0;
}

