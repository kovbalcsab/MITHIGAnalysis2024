#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "InfoManager.h"
#include "RootIOUtils.h"

static int extractIteration(const std::string &name, const std::string &prefix)
{
  if(name.rfind(prefix, 0) != 0)
    return -1;

  const std::string suffix = name.substr(prefix.size());
  if(suffix.empty())
    return -1;

  for(char c : suffix)
    if(!std::isdigit(static_cast<unsigned char>(c)))
      return -1;

  return std::stoi(suffix);
}

static bool isCompatible(const TH1 *a, const TH1 *b)
{
  if(a == nullptr || b == nullptr)
    return false;
  if(a->GetNbinsX() != b->GetNbinsX())
    return false;
  const double axMin = a->GetXaxis()->GetXmin();
  const double bxMin = b->GetXaxis()->GetXmin();
  const double axMax = a->GetXaxis()->GetXmax();
  const double bxMax = b->GetXaxis()->GetXmax();
  const double tolMin = 1e-9 * std::max(1.0, std::max(std::fabs(axMin), std::fabs(bxMin)));
  const double tolMax = 1e-9 * std::max(1.0, std::max(std::fabs(axMax), std::fabs(bxMax)));
  if(std::fabs(axMin - bxMin) > tolMin)
    return false;
  if(std::fabs(axMax - bxMax) > tolMax)
    return false;
  return true;
}

static double computeChi2(const TH1 *measured, const TH1 *refolded, int &ndf)
{
  ndf = 0;
  double chi2 = 0.0;

  for(int i = 1; i <= measured->GetNbinsX(); i++)
  {
    const double m = measured->GetBinContent(i);
    const double r = refolded->GetBinContent(i);
    const double em = measured->GetBinError(i);
    const double er = refolded->GetBinError(i);
    const double variance = em * em + er * er;

    if(variance <= 0)
      continue;

    const double diff = m - r;
    chi2 += diff * diff / variance;
    ndf++;
  }

  return chi2;
}

static std::string stripExtension(const std::string &path)
{
  const size_t slashPos = path.find_last_of('/');
  const size_t dotPos = path.find_last_of('.');
  if(dotPos == std::string::npos)
    return path;
  if(slashPos != std::string::npos && dotPos < slashPos)
    return path;
  return path.substr(0, dotPos);
}

static void updatePositiveMinAndMax(const TH1 *h, double &minPos, double &maxVal)
{
  if(!h)
    return;
  for(int i = 1; i <= h->GetNbinsX(); ++i)
  {
    const double y = h->GetBinContent(i);
    if(y > 0.0 && y < minPos)
      minPos = y;
    if(y > maxVal)
      maxVal = y;
  }
}

int main(int argc, char **argv)
{
  CommandLine CL(argc, argv);

  const std::string inputFileName = CL.Get("InputFile", "unfold_svd.root");
  const std::string outputFileName = CL.Get("OutputFile", "unfold_svd_validation.root");
  const std::string measuredName = CL.Get("MeasuredHist", "hMeasured");
  const std::string refoldedPrefix = CL.Get("RefoldedPrefix", "hRefolded_iter");
  const std::string unfoldedPrefix = CL.Get("UnfoldedPrefix", "hUnfolded_iter");
  const std::string validationDataFile = CL.Get("ValidationDataFile", "");
  const std::string validationTreeName = CL.Get("ValidationTreeName", "OutputTree");
  const std::string validationVarName = CL.Get("ValidationVarName", "");
  const int validationQuarter = CL.GetInt("ValidationQuarter", -1);
  int trainingQuarter = CL.GetInt("TrainingQuarter", -1);
  const std::string regularizationTreeName = CL.Get("RegularizationTree", "RegularizationTree");
  const int unfoldedRatioDenominatorIteration = CL.GetInt("UnfoldedRatioDenominatorIteration", 1);
  const int maxIterations = CL.GetInt("MaxIterations", -1);

  if(validationQuarter < -1 || validationQuarter > 3)
  {
    std::cerr << "ValidationQuarter must be -1 or in [0,3]." << std::endl;
    return -1;
  }
  if(!validationDataFile.empty() && validationVarName.empty())
  {
    std::cerr << "ValidationVarName must be provided when ValidationDataFile is set." << std::endl;
    return -1;
  }
  if(!validationDataFile.empty() && validationQuarter < 0)
  {
    std::cerr << "ValidationQuarter must be in [0,3] when ValidationDataFile is set." << std::endl;
    return -1;
  }

  TFile *inputFile = RootIOUtils::OpenFileOrNull(inputFileName, "READ", "input file");
  if(inputFile == nullptr)
  {
    return -1;
  }

  TH1 *measuredFromFile = RootIOUtils::GetObjectOrNull<TH1>(inputFile, measuredName, "input file: " + inputFileName);
  if(!measuredFromFile)
  {
    RootIOUtils::CloseAndDeleteFile(inputFile);
    return -1;
  }
  TH1 *measured = measuredFromFile;
  TH1D *validationMeasured = nullptr;
  std::string measuredLabel = measuredName;

  if(trainingQuarter < 0)
    trainingQuarter = RootIOUtils::InferCutParameterIntOrDefault(inputFile, "DataQuarter", -1);

  if(!validationDataFile.empty())
  {
    validationMeasured = RootIOUtils::LoadTreeBranchHistogramOrNull(
      validationDataFile, validationTreeName, validationVarName,
      "hMeasuredValidation",
      "hMeasuredValidation;" + validationVarName + ";Entries",
      measuredFromFile->GetNbinsX(),
      measuredFromFile->GetXaxis()->GetXmin(),
      measuredFromFile->GetXaxis()->GetXmax(),
      validationQuarter, true, true, "validation file");
    if(validationMeasured == nullptr)
    {
      RootIOUtils::CloseAndDeleteFile(inputFile);
      return -1;
    }
    measured = validationMeasured;
    measuredLabel = "ValidationData";
  }

  const int effectiveValidationQuarter = (!validationDataFile.empty()) ? validationQuarter : trainingQuarter;
  if(trainingQuarter >= 0 && effectiveValidationQuarter >= 0 && trainingQuarter == effectiveValidationQuarter)
  {
    std::cerr << "Warning: validation uses the same quarter (" << trainingQuarter
              << ") as the unfolding input." << std::endl;
  }

  const std::string unfoldedRatioDenName = unfoldedPrefix + std::to_string(unfoldedRatioDenominatorIteration);
  TH1 *unfoldedRatioDenominator = RootIOUtils::GetObjectOrNull<TH1>(inputFile, unfoldedRatioDenName, "input file");
  if(!unfoldedRatioDenominator)
  {
    RootIOUtils::CloseAndDeleteFile(inputFile);
    return -1;
  }

  std::vector<int> kTermForIteration;
  std::vector<double> covTraceForIteration;
  std::vector<double> covMinForIteration;
  std::vector<double> covMaxForIteration;
  TTree *regTree = RootIOUtils::GetObjectOrNull<TTree>(inputFile, regularizationTreeName, "input file");
  if(regTree)
  {
    int regIteration = 0;
    int regKTerm = 0;
    double regCovTrace = std::numeric_limits<double>::quiet_NaN();
    double regCovDiagMin = std::numeric_limits<double>::quiet_NaN();
    double regCovDiagMax = std::numeric_limits<double>::quiet_NaN();

    if(regTree->GetBranch("Iteration"))
      regTree->SetBranchAddress("Iteration", &regIteration);
    if(regTree->GetBranch("KTerm"))
      regTree->SetBranchAddress("KTerm", &regKTerm);
    if(regTree->GetBranch("CovTrace"))
      regTree->SetBranchAddress("CovTrace", &regCovTrace);
    if(regTree->GetBranch("CovDiagMin"))
      regTree->SetBranchAddress("CovDiagMin", &regCovDiagMin);
    if(regTree->GetBranch("CovDiagMax"))
      regTree->SetBranchAddress("CovDiagMax", &regCovDiagMax);

    const Long64_t nEntries = regTree->GetEntries();
    const size_t allocSize = static_cast<size_t>(nEntries) + 1;
    kTermForIteration.assign(allocSize, -1);
    covTraceForIteration.assign(allocSize, std::numeric_limits<double>::quiet_NaN());
    covMinForIteration.assign(allocSize, std::numeric_limits<double>::quiet_NaN());
    covMaxForIteration.assign(allocSize, std::numeric_limits<double>::quiet_NaN());

    for(Long64_t i = 0; i < nEntries; ++i)
    {
      regTree->GetEntry(i);
      if(regIteration >= 0 && static_cast<size_t>(regIteration) < allocSize)
      {
        kTermForIteration[static_cast<size_t>(regIteration)] = regKTerm;
        covTraceForIteration[static_cast<size_t>(regIteration)] = regCovTrace;
        covMinForIteration[static_cast<size_t>(regIteration)] = regCovDiagMin;
        covMaxForIteration[static_cast<size_t>(regIteration)] = regCovDiagMax;
      }
    }
  }

  TH2 *responseMatrix = RootIOUtils::GetObjectOrNull<TH2>(inputFile, "hResponseMatrix", "input file");

  std::vector<int> iterations;
  TIter nextKey(inputFile->GetListOfKeys());
  TKey *key = nullptr;
  while((key = static_cast<TKey *>(nextKey())) != nullptr)
  {
    const std::string objName = key->GetName();
    const int iter = extractIteration(objName, refoldedPrefix);
    if(iter < 1)
      continue;
    if(maxIterations > 0 && iter > maxIterations)
      continue;
    iterations.push_back(iter);
  }

  if(iterations.empty())
  {
    std::cerr << "No refolded histograms found with prefix: " << refoldedPrefix << std::endl;
    RootIOUtils::CloseAndDeleteFile(inputFile);
    return -1;
  }

  std::sort(iterations.begin(), iterations.end());
  iterations.erase(std::unique(iterations.begin(), iterations.end()), iterations.end());

  TFile *outputFile = RootIOUtils::OpenFileOrNull(outputFileName, "RECREATE", "output file");
  if(outputFile == nullptr)
  {
    RootIOUtils::CloseAndDeleteFile(inputFile);
    return -1;
  }
  TTimeStamp *currentTime = new TTimeStamp();
  GeneralInfoManager man(outputFile, "InfoDir", false);
  man.AddSourceFile(inputFileName, currentTime);
  if(!validationDataFile.empty())
    man.AddSourceFile(validationDataFile, currentTime);
  man.AddCutParameter("UnfoldedRatioDenominatorIteration", unfoldedRatioDenominatorIteration, currentTime);
  man.AddCutParameter("MaxIterations", maxIterations, currentTime);
  man.AddCutParameter("TrainingQuarter", trainingQuarter, currentTime);
  man.AddCutParameter("ValidationQuarter", effectiveValidationQuarter, currentTime);

  const int nIter = static_cast<int>(iterations.size());
  TH1D *hChi2 = new TH1D("hChi2VsIter", ";Iteration;#chi^{2}", nIter, 0.5, nIter + 0.5);
  TH1D *hChi2NDF = new TH1D("hChi2NDFVsIter", ";Iteration;#chi^{2}/NDF", nIter, 0.5, nIter + 0.5);
  TH1D *hPValue = new TH1D("hPValueVsIter", ";Iteration;p-value", nIter, 0.5, nIter + 0.5);

  TTree *resultTree = new TTree("Chi2Tree", "Diagnostic chi2 comparison between measured and refolded distributions");
  int iteration = 0;
  int ndf = 0;
  int kTerm = -1;
  double chi2 = 0.0;
  double chi2NDF = 0.0;
  double pvalue = 0.0;
  double covTrace = std::numeric_limits<double>::quiet_NaN();
  double covDiagMin = std::numeric_limits<double>::quiet_NaN();
  double covDiagMax = std::numeric_limits<double>::quiet_NaN();
  resultTree->Branch("Iteration", &iteration, "Iteration/I");
  resultTree->Branch("KTerm", &kTerm, "KTerm/I");
  resultTree->Branch("NDF", &ndf, "NDF/I");
  resultTree->Branch("Chi2", &chi2, "Chi2/D");
  resultTree->Branch("Chi2NDF", &chi2NDF, "Chi2NDF/D");
  resultTree->Branch("PValue", &pvalue, "PValue/D");
  resultTree->Branch("CovTrace", &covTrace, "CovTrace/D");
  resultTree->Branch("CovDiagMin", &covDiagMin, "CovDiagMin/D");
  resultTree->Branch("CovDiagMax", &covDiagMax, "CovDiagMax/D");

  double bestChi2NDF = std::numeric_limits<double>::infinity();
  int bestIteration = -1;
  std::vector<int> storedIterations;
  std::vector<int> storedKTerms;
  std::vector<TH1D *> unfoldedAbsHists;
  std::vector<TH1D *> refoldedAbsHists;
  std::vector<TH1D *> unfoldedRatioHists;
  std::vector<TH1D *> refoldedRatioHists;

  for(int idx = 0; idx < nIter; idx++)
  {
    iteration = iterations[idx];
    const std::string refoldedName = refoldedPrefix + std::to_string(iteration);
    const std::string unfoldedName = unfoldedPrefix + std::to_string(iteration);
    TH1 *refolded = dynamic_cast<TH1 *>(inputFile->Get(refoldedName.c_str()));
    TH1 *unfolded = dynamic_cast<TH1 *>(inputFile->Get(unfoldedName.c_str()));
    if(!refolded)
    {
      std::cerr << "Missing histogram: " << refoldedName << std::endl;
      continue;
    }
    if(!unfolded)
    {
      std::cerr << "Missing histogram: " << unfoldedName << std::endl;
      continue;
    }
    if(!isCompatible(measured, refolded))
    {
      std::cerr << "Incompatible binning between " << measuredLabel << " and " << refoldedName << std::endl;
      continue;
    }
    if(!isCompatible(measured, unfolded))
    {
      std::cerr << "Incompatible binning between " << measuredLabel << " and " << unfoldedName << std::endl;
      continue;
    }
    if(!isCompatible(unfoldedRatioDenominator, unfolded))
    {
      std::cerr << "Incompatible binning between " << unfoldedRatioDenName << " and " << unfoldedName << std::endl;
      continue;
    }

    chi2 = computeChi2(measured, refolded, ndf);
    chi2NDF = (ndf > 0) ? chi2 / ndf : 0.0;
    pvalue = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0.0;
    kTerm = (iteration >= 0 && static_cast<size_t>(iteration) < kTermForIteration.size()) ?
      kTermForIteration[static_cast<size_t>(iteration)] : -1;
    covTrace = (iteration >= 0 && static_cast<size_t>(iteration) < covTraceForIteration.size()) ?
      covTraceForIteration[static_cast<size_t>(iteration)] : std::numeric_limits<double>::quiet_NaN();
    covDiagMin = (iteration >= 0 && static_cast<size_t>(iteration) < covMinForIteration.size()) ?
      covMinForIteration[static_cast<size_t>(iteration)] : std::numeric_limits<double>::quiet_NaN();
    covDiagMax = (iteration >= 0 && static_cast<size_t>(iteration) < covMaxForIteration.size()) ?
      covMaxForIteration[static_cast<size_t>(iteration)] : std::numeric_limits<double>::quiet_NaN();

    hChi2->SetBinContent(idx + 1, chi2);
    hChi2NDF->SetBinContent(idx + 1, chi2NDF);
    hPValue->SetBinContent(idx + 1, pvalue);

    hChi2->GetXaxis()->SetBinLabel(idx + 1, std::to_string(iteration).c_str());
    hChi2NDF->GetXaxis()->SetBinLabel(idx + 1, std::to_string(iteration).c_str());
    hPValue->GetXaxis()->SetBinLabel(idx + 1, std::to_string(iteration).c_str());

    resultTree->Fill();

    std::cout << "Iteration " << iteration
              << " : chi2 = " << chi2
              << ", NDF = " << ndf
              << ", chi2/NDF = " << chi2NDF
              << ", p-value = " << pvalue;
    if(kTerm > 0)
      std::cout << ", k-term = " << kTerm;
    std::cout << std::endl;

    if(ndf > 0 && chi2NDF < bestChi2NDF)
    {
      bestChi2NDF = chi2NDF;
      bestIteration = iteration;
    }

    const std::string unfoldedRatioName = "hUnfoldedToMeasuredRatio_iter" + std::to_string(iteration);
    TH1D *unfoldedRatio = dynamic_cast<TH1D *>(unfolded->Clone(unfoldedRatioName.c_str()));
    if(!unfoldedRatio)
    {
      std::cerr << "Failed to build unfolded ratio histogram for iteration " << iteration << std::endl;
      continue;
    }

    const std::string unfoldedAbsName = "hUnfoldedAbs_iter" + std::to_string(iteration);
    TH1D *unfoldedAbs = dynamic_cast<TH1D *>(unfolded->Clone(unfoldedAbsName.c_str()));
    if(!unfoldedAbs)
    {
      std::cerr << "Failed to build unfolded absolute histogram for iteration " << iteration << std::endl;
      delete unfoldedRatio;
      continue;
    }

    const std::string refoldedRatioName = "hRefoldedToMeasuredRatio_iter" + std::to_string(iteration);
    TH1D *refoldedRatio = dynamic_cast<TH1D *>(refolded->Clone(refoldedRatioName.c_str()));
    if(!refoldedRatio)
    {
      std::cerr << "Failed to build refolded ratio histogram for iteration " << iteration << std::endl;
      delete unfoldedRatio;
      delete unfoldedAbs;
      continue;
    }

    const std::string refoldedAbsName = "hRefoldedAbs_iter" + std::to_string(iteration);
    TH1D *refoldedAbs = dynamic_cast<TH1D *>(refolded->Clone(refoldedAbsName.c_str()));
    if(!refoldedAbs)
    {
      std::cerr << "Failed to build refolded absolute histogram for iteration " << iteration << std::endl;
      delete unfoldedRatio;
      delete unfoldedAbs;
      delete refoldedRatio;
      continue;
    }

    unfoldedAbs->SetDirectory(nullptr);
    refoldedAbs->SetDirectory(nullptr);
    unfoldedRatio->SetDirectory(nullptr);
    refoldedRatio->SetDirectory(nullptr);
    unfoldedAbs->SetTitle(("Unfolded absolute overlaps;" + std::string(measured->GetXaxis()->GetTitle()) + ";Entries").c_str());
    refoldedAbs->SetTitle(("Refolded absolute overlaps;" + std::string(measured->GetXaxis()->GetTitle()) + ";Entries").c_str());
    unfoldedRatio->SetTitle(("Unfolded/iter" + std::to_string(unfoldedRatioDenominatorIteration) + " ratio;" + std::string(measured->GetXaxis()->GetTitle()) + ";Unfolded / Unfolded_{den}").c_str());
    refoldedRatio->SetTitle(("Refolded/Measured ratio;" + std::string(measured->GetXaxis()->GetTitle()) + ";Refolded / Measured").c_str());
    for(int bin = 1; bin <= unfoldedRatio->GetNbinsX(); ++bin)
    {
      const double u = unfolded->GetBinContent(bin);
      const double du = unfoldedRatioDenominator->GetBinContent(bin);
      const double f = refolded->GetBinContent(bin);
      const double m = measured->GetBinContent(bin);
      const double eu = unfolded->GetBinError(bin);
      const double edu = unfoldedRatioDenominator->GetBinError(bin);
      const double ef = refolded->GetBinError(bin);
      const double em = measured->GetBinError(bin);

      if(du <= 0.0)
      {
        unfoldedRatio->SetBinContent(bin, 0.0);
        unfoldedRatio->SetBinError(bin, 0.0);
      }
      else
      {
        const double ru = u / du;
        double relU2 = 0.0;
        double relDU2 = 0.0;
        if(u > 0)
          relU2 = (eu / u) * (eu / u);
        if(du > 0)
          relDU2 = (edu / du) * (edu / du);
        const double eru = ru * std::sqrt(relU2 + relDU2);

        unfoldedRatio->SetBinContent(bin, ru);
        unfoldedRatio->SetBinError(bin, eru);
      }

      if(m <= 0.0)
      {
        refoldedRatio->SetBinContent(bin, 0.0);
        refoldedRatio->SetBinError(bin, 0.0);
        continue;
      }

      const double rf = f / m;
      double relF2 = 0.0;
      double relM2 = 0.0;
      if(f > 0)
        relF2 = (ef / f) * (ef / f);
      relM2 = (em / m) * (em / m);
      const double erf = rf * std::sqrt(relF2 + relM2);

      refoldedRatio->SetBinContent(bin, rf);
      refoldedRatio->SetBinError(bin, erf);
    }

    storedIterations.push_back(iteration);
    storedKTerms.push_back(kTerm);
    unfoldedAbsHists.push_back(unfoldedAbs);
    refoldedAbsHists.push_back(refoldedAbs);
    unfoldedRatioHists.push_back(unfoldedRatio);
    refoldedRatioHists.push_back(refoldedRatio);
  }

  TH1D *measuredForPlot = dynamic_cast<TH1D *>(measured->Clone("hMeasuredForOverlay"));
  if(measuredForPlot)
  {
    measuredForPlot->SetDirectory(nullptr);
    measuredForPlot->SetLineColor(kAzure + 2);
    measuredForPlot->SetMarkerColor(kAzure + 2);
    measuredForPlot->SetMarkerStyle(20);
    measuredForPlot->SetLineWidth(2);
    measuredForPlot->SetStats(0);
  }

  gStyle->SetOptStat(0);

  const double topToBottomRatio = 2.2;
  const double lowerFrac = 1.0 / (1.0 + topToBottomRatio);
  const int palette[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1, kSpring + 5, kPink + 7, kTeal + 3};
  const int nPalette = static_cast<int>(sizeof(palette) / sizeof(int));

  auto makeComparisonCanvas = [&](const std::string &canvasName,
                                  const std::string &canvasTitle,
                                  const std::vector<TH1D *> &absHists,
                                  const std::vector<TH1D *> &ratioHists,
                                  const std::string &legendPrefix,
                                  bool zoomToLowX) -> TCanvas *
  {
    if(ratioHists.empty() || measuredForPlot == nullptr)
      return nullptr;

    TH1D *measuredLocal = dynamic_cast<TH1D *>(measuredForPlot->Clone((canvasName + "_measured").c_str()));
    if(!measuredLocal)
      return nullptr;
    measuredLocal->SetDirectory(nullptr);
    measuredLocal->SetStats(0);

    std::vector<TH1D *> absLocal(absHists.size(), nullptr);
    std::vector<TH1D *> ratioLocal(ratioHists.size(), nullptr);
    for(size_t i = 0; i < absHists.size(); ++i)
    {
      if(absHists[i])
      {
        absLocal[i] = dynamic_cast<TH1D *>(absHists[i]->Clone((canvasName + "_abs_" + std::to_string(i)).c_str()));
        if(absLocal[i])
        {
          absLocal[i]->SetDirectory(nullptr);
          absLocal[i]->SetStats(0);
        }
      }
      if(ratioHists[i])
      {
        ratioLocal[i] = dynamic_cast<TH1D *>(ratioHists[i]->Clone((canvasName + "_ratio_" + std::to_string(i)).c_str()));
        if(ratioLocal[i])
        {
          ratioLocal[i]->SetDirectory(nullptr);
          ratioLocal[i]->SetStats(0);
        }
      }
    }

    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 900, 900);

    TPad *padTop = new TPad((canvasName + "_top").c_str(), "", 0.0, lowerFrac, 1.0, 1.0);
    TPad *padBottom = new TPad((canvasName + "_bottom").c_str(), "", 0.0, 0.0, 1.0, lowerFrac);

    padTop->SetLeftMargin(0.12);
    padTop->SetRightMargin(0.04);
    padTop->SetTopMargin(0.08);
    padTop->SetBottomMargin(0.0);
    padTop->SetLogy();

    padBottom->SetLeftMargin(0.12);
    padBottom->SetRightMargin(0.04);
    padBottom->SetTopMargin(0.0);
    padBottom->SetBottomMargin(0.30);

    canvas->cd();
    padTop->Draw();
    padBottom->Draw();

    double minPos = std::numeric_limits<double>::infinity();
    double maxVal = 0.0;
    updatePositiveMinAndMax(measuredLocal, minPos, maxVal);
    for(size_t i = 0; i < absHists.size(); ++i)
      updatePositiveMinAndMax(absLocal[i], minPos, maxVal);

    if(zoomToLowX)
    {
      measuredLocal->GetXaxis()->SetRangeUser(0.0, 25.0);
      for(TH1D *h : absLocal)
        h->GetXaxis()->SetRangeUser(0.0, 25.0);
      for(TH1D *h : ratioLocal)
        h->GetXaxis()->SetRangeUser(0.0, 25.0);
    }
    else
    {
      measuredLocal->GetXaxis()->UnZoom();
      for(TH1D *h : absLocal)
        h->GetXaxis()->UnZoom();
      for(TH1D *h : ratioLocal)
        h->GetXaxis()->UnZoom();
    }

    if(minPos < std::numeric_limits<double>::infinity())
      measuredLocal->SetMinimum(std::max(1e-2, minPos * 0.5));
    if(maxVal > 0.0)
      measuredLocal->SetMaximum(maxVal * 3.0);
    measuredLocal->GetXaxis()->SetLabelSize(0.0);
    measuredLocal->GetXaxis()->SetTitleSize(0.0);

    const double topHeight = 1.0 - lowerFrac;
    const double bottomHeight = lowerFrac;
    const double topText = 0.045;
    const double bottomText = topText * (topHeight / bottomHeight);
    measuredLocal->GetYaxis()->SetLabelSize(topText);
    measuredLocal->GetYaxis()->SetTitleSize(topText);

    padTop->cd();
    measuredLocal->Draw("E");
    TLegend *absLegend = new TLegend(0.52, 0.55, 0.88, 0.88);
    absLegend->SetBorderSize(0);
    absLegend->SetFillStyle(0);
    absLegend->AddEntry(measuredLocal, "Measured", "lep");

    padBottom->cd();

    bool drawnRatio = false;
    for(size_t i = 0; i < ratioLocal.size(); ++i)
    {
      TH1D *absHist = absLocal[i];
      TH1D *ratio = ratioLocal[i];
      if(!absHist || !ratio)
        continue;
      const int color = palette[static_cast<int>(i % nPalette)];

      padTop->cd();
      absHist->SetLineColor(color);
      absHist->SetMarkerColor(color);
      absHist->SetMarkerStyle(24 + static_cast<int>(i % 5));
      absHist->SetLineWidth(2);
      absHist->Draw("HIST SAME");
      std::string label = legendPrefix + " iter " + std::to_string(storedIterations[i]);
      if(storedKTerms[i] > 0)
        label += " (k=" + std::to_string(storedKTerms[i]) + ")";
      absLegend->AddEntry(absHist, label.c_str(), "l");

      padBottom->cd();
      ratio->SetLineColor(color);
      ratio->SetMarkerColor(color);
      ratio->SetMarkerStyle(20 + static_cast<int>(i % 5));
      ratio->SetLineWidth(2);
      ratio->GetYaxis()->SetRangeUser(0.0, 2.5);
      ratio->GetYaxis()->SetLabelSize(bottomText);
      ratio->GetYaxis()->SetTitleSize(bottomText);
      ratio->GetXaxis()->SetLabelSize(bottomText);
      ratio->GetXaxis()->SetTitleSize(bottomText);

      if(!drawnRatio)
        ratio->Draw("HIST");
      else
        ratio->Draw("HIST SAME");
      drawnRatio = true;
      ratio->SetTitle("");
      ratio->GetXaxis()->SetTitle("");
      ratio->GetYaxis()->SetTitle("");
      ratio->GetYaxis()->ChangeLabel(-1, -1, 0);
    }

    padTop->cd();
    absLegend->Draw();
    padBottom->cd();

    return canvas;
  };

  TCanvas *cUnfoldedComparison = makeComparisonCanvas("cUnfoldedComparison", "Unfolded Comparison", unfoldedAbsHists, unfoldedRatioHists, "Unfolded", false);
  TCanvas *cRefoldedComparison = makeComparisonCanvas("cRefoldedComparison", "Refolded Comparison", refoldedAbsHists, refoldedRatioHists, "Refolded", false);
  TCanvas *cUnfoldedComparisonX0to25 = makeComparisonCanvas("cUnfoldedComparison_X0to25", "Unfolded Comparison (x=[0,25])", unfoldedAbsHists, unfoldedRatioHists, "Unfolded", true);
  TCanvas *cRefoldedComparisonX0to25 = makeComparisonCanvas("cRefoldedComparison_X0to25", "Refolded Comparison (x=[0,25])", refoldedAbsHists, refoldedRatioHists, "Refolded", true);

  TCanvas *cChi2NDF = new TCanvas("cChi2NDF", "Chi2/NDF vs Iteration", 800, 600);
  hChi2NDF->SetStats(0);
  hChi2NDF->SetLineColor(kBlue + 1);
  hChi2NDF->SetMarkerColor(kBlue + 1);
  hChi2NDF->SetMarkerStyle(20);
  hChi2NDF->SetLineWidth(2);
  hChi2NDF->Draw("HIST");
  hChi2NDF->Draw("P SAME");

  TCanvas *cResponseMatrix = nullptr;
  if(responseMatrix)
  {
    cResponseMatrix = new TCanvas("cResponseMatrix", "Unfolding Response Matrix", 800, 700);
    responseMatrix->SetStats(0);
    responseMatrix->SetTitle("Unfolding response matrix;Measured;True");
    responseMatrix->Draw("COLZ");
  }

  const std::string outputBase = stripExtension(outputFileName);
  if(cUnfoldedComparison)
    cUnfoldedComparison->SaveAs((outputBase + "_unfolded_comparison.pdf").c_str());
  if(cRefoldedComparison)
    cRefoldedComparison->SaveAs((outputBase + "_refolded_comparison.pdf").c_str());
  if(cUnfoldedComparisonX0to25)
    cUnfoldedComparisonX0to25->SaveAs((outputBase + "_unfolded_comparison_x0_25.pdf").c_str());
  if(cRefoldedComparisonX0to25)
    cRefoldedComparisonX0to25->SaveAs((outputBase + "_refolded_comparison_x0_25.pdf").c_str());
  if(cChi2NDF)
    cChi2NDF->SaveAs((outputBase + "_chi2ndf.pdf").c_str());
  if(cResponseMatrix)
    cResponseMatrix->SaveAs((outputBase + "_response_matrix.pdf").c_str());

  outputFile->cd();
  hChi2->Write();
  hChi2NDF->Write();
  hPValue->Write();
  resultTree->Write();
  if(measuredForPlot)
    measuredForPlot->Write();
  for(TH1D *hist : unfoldedAbsHists)
    hist->Write();
  for(TH1D *hist : refoldedAbsHists)
    hist->Write();
  for(TH1D *ratio : unfoldedRatioHists)
    ratio->Write();
  for(TH1D *ratio : refoldedRatioHists)
    ratio->Write();
  if(cUnfoldedComparison)
    cUnfoldedComparison->Write();
  if(cRefoldedComparison)
    cRefoldedComparison->Write();
  if(cUnfoldedComparisonX0to25)
    cUnfoldedComparisonX0to25->Write();
  if(cRefoldedComparisonX0to25)
    cRefoldedComparisonX0to25->Write();
  if(cChi2NDF)
    cChi2NDF->Write();
  if(cResponseMatrix)
    cResponseMatrix->Write();
  man.SaveToFile();
  RootIOUtils::CloseAndDeleteFile(outputFile);
  RootIOUtils::CloseAndDeleteFile(inputFile);

  if(bestIteration > 0)
    std::cout << "Lowest diagnostic chi2/NDF iteration (not an optimization recommendation): " << bestIteration << " (" << bestChi2NDF << ")" << std::endl;
  else
    std::cout << "No valid chi2 result was produced." << std::endl;

  delete currentTime;
  if(validationMeasured)
    delete validationMeasured;
  return 0;
}
