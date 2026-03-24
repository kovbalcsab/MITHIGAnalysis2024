#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLeaf.h"
#include "TLatex.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "CommandLine.h"

static double ComputeChi2(const TH1D *target, const TH1D *rescaled, int &ndf)
{
  ndf = 0;
  if(target == nullptr || rescaled == nullptr)
    return 0.0;

  const int nBins = std::min(target->GetNbinsX(), rescaled->GetNbinsX());
  double chi2 = 0.0;
  for(int i = 1; i <= nBins; ++i)
  {
    const double t = target->GetBinContent(i);
    const double r = rescaled->GetBinContent(i);
    const double et = target->GetBinError(i);
    const double er = rescaled->GetBinError(i);
    const double variance = et * et + er * er;
    if(variance <= 0.0)
      continue;

    const double diff = t - r;
    chi2 += diff * diff / variance;
    ndf++;
  }
  return chi2;
}

static bool ReadAB(const std::string &parameterFile, double &a, double &b)
{
  std::ifstream in(parameterFile);
  if(!in.is_open())
  {
    std::cerr << "Cannot open parameter file: " << parameterFile << std::endl;
    return false;
  }

  bool foundA = false;
  bool foundB = false;
  std::string line;
  while(std::getline(in, line))
  {
    if((line.rfind("a:", 0) == 0 || line.rfind("A:", 0) == 0) && std::sscanf(line.c_str(), "%*[^:]: %lf", &a) == 1)
      foundA = true;
    if((line.rfind("b:", 0) == 0 || line.rfind("B:", 0) == 0) && std::sscanf(line.c_str(), "%*[^:]: %lf", &b) == 1)
      foundB = true;
  }

  if(!foundA || !foundB)
  {
    std::cerr << "Failed to read both a and b from: " << parameterFile << std::endl;
    return false;
  }
  return true;
}

static TH1D *FillRawTreeHistogram(TTree *tree, const std::string &varName, const std::string &histName,
                                  int nBins, double xMin, double xMax)
{
  if(tree == nullptr)
    return nullptr;
  if(tree->GetBranch(varName.c_str()) == nullptr)
    return nullptr;

  TH1D *h = new TH1D(histName.c_str(), ";E [GeV];1/N dN/dE", nBins, xMin, xMax);
  h->Sumw2();
  tree->Draw((varName + ">>" + histName).c_str(), "", "goff");
  return h;
}

static TH1D *FillRescaledEventByEventHistogram(TTree *tree, const std::string &varName, const std::string &histName,
                                               int nBins, double xMin, double xMax, double a, double b)
{
  if(tree == nullptr)
    return nullptr;
  TBranch *branch = tree->GetBranch(varName.c_str());
  if(branch == nullptr)
    return nullptr;

  TLeaf *leaf = branch->GetLeaf(varName.c_str());
  if(leaf == nullptr)
    leaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->First());
  if(leaf == nullptr)
    return nullptr;

  TH1D *h = new TH1D(histName.c_str(), ";E [GeV];1/N dN/dE", nBins, xMin, xMax);
  h->Sumw2();

  const std::string leafType = leaf->GetTypeName() ? std::string(leaf->GetTypeName()) : "";
  const bool isFloat = (leafType == "Float_t" || leafType == "float");

  double valueD = 0.0;
  float valueF = 0.0f;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(varName.c_str(), 1);
  if(isFloat)
    tree->SetBranchAddress(varName.c_str(), &valueF);
  else
    tree->SetBranchAddress(varName.c_str(), &valueD);

  const Long64_t nEntries = tree->GetEntries();
  for(Long64_t i = 0; i < nEntries; ++i)
  {
    tree->GetEntry(i);
    const double e = isFloat ? static_cast<double>(valueF) : valueD;
    const double mapped = a * e + b * e * e;
    h->Fill(mapped);
  }

  tree->SetBranchStatus("*", 1);
  tree->ResetBranchAddresses();
  return h;
}

int main(int argc, char **argv)
{
  CommandLine CL(argc, argv);

  const std::string targetFileName = CL.Get("TargetFile", "");
  const std::string fileToFitName = CL.Get("FileToFit", "");
  const std::string parameterFileName = CL.Get("ParameterFile", "");
  const std::string plotFileName = CL.Get("PlotFileName", "scaling_fit_result.pdf");
  const std::string treeName = CL.Get("TreeName", "OutputTree");
  const std::string varTargetName = CL.Get("VarTargetName", "HFEMaxPlus_forest");
  const std::string varFitName = CL.Get("VarFitName", "HFEMaxPlus_forest");
  const int nBins = CL.GetInt("NBins", 200);
  const double xMin = CL.GetDouble("XMin", 0.0);
  const double xMax = CL.GetDouble("XMax", 25.0);

  if(targetFileName.empty() || fileToFitName.empty() || parameterFileName.empty())
  {
    std::cerr << "Please provide --TargetFile, --FileToFit, and --ParameterFile." << std::endl;
    return -1;
  }
  if(!(xMax > xMin) || nBins <= 0)
  {
    std::cerr << "Invalid histogram settings: require XMax>XMin and NBins>0." << std::endl;
    return -1;
  }

  double a = 0.0;
  double b = 0.0;
  if(!ReadAB(parameterFileName, a, b))
    return -1;

  TFile targetFile(targetFileName.c_str(), "READ");
  if(targetFile.IsZombie())
  {
    std::cerr << "Cannot open target file: " << targetFileName << std::endl;
    return -1;
  }
  TTree *targetTree = dynamic_cast<TTree *>(targetFile.Get(treeName.c_str()));
  if(targetTree == nullptr)
  {
    std::cerr << "Cannot find tree " << treeName << " in target file." << std::endl;
    return -1;
  }

  TFile fitFile(fileToFitName.c_str(), "READ");
  if(fitFile.IsZombie())
  {
    std::cerr << "Cannot open fit file: " << fileToFitName << std::endl;
    return -1;
  }
  TTree *fitTree = dynamic_cast<TTree *>(fitFile.Get(treeName.c_str()));
  if(fitTree == nullptr)
  {
    std::cerr << "Cannot find tree " << treeName << " in fit file." << std::endl;
    return -1;
  }

  TH1D *hTargetRaw = FillRawTreeHistogram(targetTree, varTargetName, "hTargetRaw", nBins, xMin, xMax);
  TH1D *hFitRaw = FillRawTreeHistogram(fitTree, varFitName, "hFitRaw", nBins, xMin, xMax);
  TH1D *hFitRescaled = FillRescaledEventByEventHistogram(fitTree, varFitName, "hFitRescaled", nBins, xMin, xMax, a, b);

  if(hTargetRaw == nullptr || hFitRaw == nullptr || hFitRescaled == nullptr)
  {
    std::cerr << "Failed to construct one or more histograms. Check tree/branch names." << std::endl;
    delete hTargetRaw;
    delete hFitRaw;
    delete hFitRescaled;
    return -1;
  }

  const double targetInt = hTargetRaw->Integral("width");
  const double fitInt = hFitRaw->Integral("width");
  const double scaledInt = hFitRescaled->Integral("width");
  if(!(targetInt > 0.0) || !(fitInt > 0.0) || !(scaledInt > 0.0))
  {
    std::cerr << "At least one histogram integral is non-positive." << std::endl;
    delete hTargetRaw;
    delete hFitRaw;
    delete hFitRescaled;
    return -1;
  }

  hTargetRaw->Scale(1.0 / targetInt);
  hFitRaw->Scale(1.0 / fitInt);
  hFitRescaled->Scale(1.0 / scaledInt);

  int chi2NDF = 0;
  const double chi2 = ComputeChi2(hTargetRaw, hFitRescaled, chi2NDF);
  const double chi2OverNDF = (chi2NDF > 0) ? (chi2 / chi2NDF) : 0.0;

  hTargetRaw->SetLineColor(kRed + 1);
  hTargetRaw->SetMarkerColor(kRed + 1);
  hTargetRaw->SetMarkerStyle(20);
  hTargetRaw->SetLineWidth(2);

  hFitRaw->SetLineColor(kBlue + 1);
  hFitRaw->SetMarkerColor(kBlue + 1);
  hFitRaw->SetMarkerStyle(21);
  hFitRaw->SetLineWidth(2);

  hFitRescaled->SetLineColor(kBlack);
  hFitRescaled->SetLineWidth(3);

  TCanvas *c = new TCanvas("cScalingSimple", "Scaling comparison", 900, 700);
  c->SetLogy();
  hFitRescaled->SetTitle("Raw and event-rescaled distributions");
  hFitRescaled->GetYaxis()->SetRangeUser(1e-6, 1.0);
  hFitRescaled->Draw("E");
  hTargetRaw->Draw("E SAME");
  hFitRaw->Draw("E SAME");

  TLegend *leg = new TLegend(0.54, 0.66, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hTargetRaw, "Target raw", "lep");
  leg->AddEntry(hFitRaw, "FileToFit raw", "lep");
  leg->AddEntry(hFitRescaled, "FileToFit event-by-event rescaled", "l");
  leg->Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.03);
  label.DrawLatex(0.13, 0.84,
                  Form("#chi^{2}/NDF (Target vs Rescaled) = %.4g / %d = %.4g", chi2, chi2NDF, chi2OverNDF));

  std::cout << "Using parameters: a=" << a << ", b=" << b << std::endl;
  std::cout << "Target vs rescaled chi2 = " << chi2 << ", NDF = " << chi2NDF
            << ", chi2/NDF = " << chi2OverNDF << std::endl;
  c->SaveAs(plotFileName.c_str());

  delete leg;
  delete c;
  delete hTargetRaw;
  delete hFitRaw;
  delete hFitRescaled;
  return 0;
}
