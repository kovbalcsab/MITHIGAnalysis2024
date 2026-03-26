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
#include <cctype>
#include <cstdio>
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

static double computeChi2(const TH1 *measured, const TH1 *reconv, int &ndf)
{
  ndf = 0;
  double chi2 = 0.0;

  for(int i = 1; i <= measured->GetNbinsX(); i++)
  {
    const double m = measured->GetBinContent(i);
    const double r = reconv->GetBinContent(i);
    const double em = measured->GetBinError(i);
    const double er = reconv->GetBinError(i);
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

static TH1D *buildReconvolvedFromDeconvolved(const TH1 *deconv, const TH1 *kernel, const std::string &name)
{
  if(deconv == nullptr || kernel == nullptr)
    return nullptr;
  if(deconv->GetNbinsX() != kernel->GetNbinsX())
    return nullptr;

  const int N = deconv->GetNbinsX();
  std::vector<double> s(N, 0.0), k(N, 0.0);
  for(int i = 0; i < N; ++i)
  {
    s[i] = std::max(0.0, deconv->GetBinContent(i + 1));
    k[i] = std::max(0.0, kernel->GetBinContent(i + 1));
  }

  double ksum = 0.0;
  for(double v : k)
    ksum += v;
  if(ksum <= 0.0)
    return nullptr;
  for(double &v : k)
    v /= ksum;

  const double binWidth = deconv->GetXaxis()->GetBinWidth(1);
  const double xMinVal = deconv->GetXaxis()->GetXmin();
  int kc = static_cast<int>(std::round(-xMinVal / binWidth));
  kc = std::max(0, std::min(kc, N - 1));

  TH1D *out = dynamic_cast<TH1D *>(deconv->Clone(name.c_str()));
  if(!out)
    return nullptr;
  out->SetDirectory(nullptr);
  out->Reset();

  for(int i = 0; i < N; ++i)
  {
    double sum = 0.0;
    for(int j = 0; j < N; ++j)
    {
      const int ki = i - j + kc;
      if(ki >= 0 && ki < N)
        sum += s[j] * k[ki];
    }
    out->SetBinContent(i + 1, std::max(0.0, sum));
    out->SetBinError(i + 1, 0.0);
  }

  const double intg = out->Integral();
  if(intg > 0.0)
    out->Scale(1.0 / intg);

  return out;
}

int main(int argc, char **argv)
{
  CommandLine CL(argc, argv);

  const std::string inputFileName = CL.Get("InputFile", "deconvolve_fft.root");
  const std::string outputFileName = CL.Get("OutputFile", "deconvolve_fft_validation.root");
  const std::string measuredName = CL.Get("MeasuredHist", "hMeasured");
  const std::string deconvolvedPrefix = CL.Get("DeconvolvedPrefix", "hDeconvolved_iter");
  const std::string reconvolvedPrefix = CL.Get("ReconvolvedPrefix", "hReconvolved_iter");
  const std::string validationDataFile = CL.Get("ValidationDataFile", "");
  const std::string validationTreeName = CL.Get("ValidationTreeName", "OutputTree");
  const std::string validationVarName = CL.Get("ValidationVarName", "");
  const int validationQuarter = CL.GetInt("ValidationQuarter", -1);
  int trainingQuarter = CL.GetInt("TrainingQuarter", -1);
  const std::string dampingTreeName = CL.Get("RegularizationTree", "RegularizationTree");
  const std::string dampingBranchName = CL.Get("DampingBranch", "DampingFactor");
  const int unfoldedRatioDenominatorIteration = CL.GetInt("UnfoldedRatioDenominatorIteration", 1);
  const int maxIterations = CL.GetInt("MaxIterations", -1);
  const bool allowReconvolvedFallback = CL.GetBool("AllowReconvolvedFallback", false);

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
              << ") as the deconvolution input." << std::endl;
  }

  const std::string unfoldedRatioDenName = deconvolvedPrefix + std::to_string(unfoldedRatioDenominatorIteration);
  TH1 *unfoldedRatioDenominator = RootIOUtils::GetObjectOrNull<TH1>(inputFile, unfoldedRatioDenName, "input file");
  if(!unfoldedRatioDenominator)
  {
    RootIOUtils::CloseAndDeleteFile(inputFile);
    return -1;
  }

  TTree *regTreeIn = RootIOUtils::GetObjectOrNull<TTree>(inputFile, dampingTreeName, "input file");
  std::vector<double> dampingValues;
  if(regTreeIn)
  {
    int regIteration = 0;
    double damping = 0.0;
    if(regTreeIn->GetBranch("Iteration"))
      regTreeIn->SetBranchAddress("Iteration", &regIteration);
    if(regTreeIn->GetBranch(dampingBranchName.c_str()))
      regTreeIn->SetBranchAddress(dampingBranchName.c_str(), &damping);

    const Long64_t nEntries = regTreeIn->GetEntries();
    dampingValues.assign(static_cast<size_t>(nEntries) + 1, std::numeric_limits<double>::quiet_NaN());
    for(Long64_t i = 0; i < nEntries; ++i)
    {
      regTreeIn->GetEntry(i);
      if(regIteration >= 0 && static_cast<size_t>(regIteration) < dampingValues.size())
        dampingValues[static_cast<size_t>(regIteration)] = damping;
    }
  }

  TH1 *kernel = RootIOUtils::GetObjectOrNull<TH1>(inputFile, "hKernel", "input file");

  std::vector<int> iterations;
  TIter nextKey(inputFile->GetListOfKeys());
  TKey *key = nullptr;
  while((key = static_cast<TKey *>(nextKey())) != nullptr)
  {
    const std::string objName = key->GetName();
    const int iter = extractIteration(objName, reconvolvedPrefix);
    if(iter < 1)
      continue;
    if(maxIterations > 0 && iter > maxIterations)
      continue;
    iterations.push_back(iter);
  }

  if(iterations.empty())
  {
    TIter nextKeyDeconv(inputFile->GetListOfKeys());
    while((key = static_cast<TKey *>(nextKeyDeconv())) != nullptr)
    {
      const std::string objName = key->GetName();
      const int iter = extractIteration(objName, deconvolvedPrefix);
      if(iter < 1)
        continue;
      if(maxIterations > 0 && iter > maxIterations)
        continue;
      iterations.push_back(iter);
    }
  }

  if(iterations.empty())
  {
    std::cerr << "No reconvolved/deconvolved histograms found with prefixes: "
              << reconvolvedPrefix << " and " << deconvolvedPrefix << std::endl;
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
  man.AddCutParameter("AllowReconvolvedFallback", allowReconvolvedFallback ? 1.0 : 0.0, currentTime);
  man.AddCutParameter("TrainingQuarter", trainingQuarter, currentTime);
  man.AddCutParameter("ValidationQuarter", effectiveValidationQuarter, currentTime);

  const int nIter = static_cast<int>(iterations.size());
  TH1D *hChi2 = new TH1D("hChi2VsIter", ";Iteration;#chi^{2}", nIter, 0.5, nIter + 0.5);
  TH1D *hChi2NDF = new TH1D("hChi2NDFVsIter", ";Iteration;#chi^{2}/NDF", nIter, 0.5, nIter + 0.5);
  TH1D *hPValue = new TH1D("hPValueVsIter", ";Iteration;p-value", nIter, 0.5, nIter + 0.5);

  TTree *resultTree = new TTree("Chi2Tree", "Diagnostic chi2 comparison between measured and reconvolved distributions");
  int iteration = 0;
  int ndf = 0;
  double chi2 = 0.0;
  double chi2NDF = 0.0;
  double pvalue = 0.0;
  double damping = std::numeric_limits<double>::quiet_NaN();
  resultTree->Branch("Iteration", &iteration, "Iteration/I");
  resultTree->Branch("NDF", &ndf, "NDF/I");
  resultTree->Branch("Chi2", &chi2, "Chi2/D");
  resultTree->Branch("Chi2NDF", &chi2NDF, "Chi2NDF/D");
  resultTree->Branch("PValue", &pvalue, "PValue/D");
  resultTree->Branch("DampingFactor", &damping, "DampingFactor/D");

  double bestChi2NDF = std::numeric_limits<double>::infinity();
  int bestIteration = -1;
  std::vector<int> storedIterations;
  std::vector<double> storedDamping;
  std::vector<TH1D *> unfoldedAbsHists;
  std::vector<TH1D *> refoldedAbsHists;
  std::vector<TH1D *> unfoldedRatioHists;
  std::vector<TH1D *> refoldedRatioHists;

  for(int idx = 0; idx < nIter; idx++)
  {
    iteration = iterations[idx];
    const std::string reconvolvedName = reconvolvedPrefix + std::to_string(iteration);
    const std::string deconvolvedName = deconvolvedPrefix + std::to_string(iteration);
    TH1 *deconv = dynamic_cast<TH1 *>(inputFile->Get(deconvolvedName.c_str()));
    if(!deconv)
    {
      std::cerr << "Missing histogram: " << deconvolvedName << std::endl;
      continue;
    }

    TH1 *reconv = dynamic_cast<TH1 *>(inputFile->Get(reconvolvedName.c_str()));
    TH1D *reconvComputed = nullptr;
    if(!reconv)
    {
      if(!allowReconvolvedFallback)
      {
        std::cerr << "Missing histogram: " << reconvolvedName
                  << ". Set --AllowReconvolvedFallback true to synthesize reconvolution from deconvolved+kernel."
                  << std::endl;
        continue;
      }
      reconvComputed = buildReconvolvedFromDeconvolved(deconv, kernel, reconvolvedName + "_computed");
      reconv = reconvComputed;
      if(!reconv)
      {
        std::cerr << "Missing histogram " << reconvolvedName
                  << " and failed to compute reconvolution from hKernel." << std::endl;
        continue;
      }
      std::cout << "Warning: using synthesized reconvolution for iteration " << iteration
                << " (fallback path)." << std::endl;
    }

    if(!isCompatible(measured, reconv))
    {
      std::cerr << "Incompatible binning between " << measuredLabel << " and " << reconvolvedName << std::endl;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }
    if(!isCompatible(measured, deconv))
    {
      std::cerr << "Incompatible binning between " << measuredLabel << " and " << deconvolvedName << std::endl;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }
    if(!isCompatible(unfoldedRatioDenominator, deconv))
    {
      std::cerr << "Incompatible binning between " << unfoldedRatioDenName << " and " << deconvolvedName << std::endl;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }

    chi2 = computeChi2(measured, reconv, ndf);
    chi2NDF = (ndf > 0) ? chi2 / ndf : 0.0;
    pvalue = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0.0;
    damping = (iteration >= 0 && static_cast<size_t>(iteration) < dampingValues.size()) ?
      dampingValues[static_cast<size_t>(iteration)] : std::numeric_limits<double>::quiet_NaN();

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
    if(std::isfinite(damping))
      std::cout << ", damping = " << damping;
    std::cout << std::endl;

    if(ndf > 0 && chi2NDF < bestChi2NDF)
    {
      bestChi2NDF = chi2NDF;
      bestIteration = iteration;
    }

    const std::string unfoldedRatioName = "hUnfoldedToMeasuredRatio_iter" + std::to_string(iteration);
    TH1D *unfoldedRatio = dynamic_cast<TH1D *>(deconv->Clone(unfoldedRatioName.c_str()));
    if(!unfoldedRatio)
    {
      std::cerr << "Failed to build deconvolved ratio histogram for iteration " << iteration << std::endl;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }

    const std::string unfoldedAbsName = "hUnfoldedAbs_iter" + std::to_string(iteration);
    TH1D *unfoldedAbs = dynamic_cast<TH1D *>(deconv->Clone(unfoldedAbsName.c_str()));
    if(!unfoldedAbs)
    {
      std::cerr << "Failed to build deconvolved absolute histogram for iteration " << iteration << std::endl;
      delete unfoldedRatio;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }

    const std::string refoldedRatioName = "hRefoldedToMeasuredRatio_iter" + std::to_string(iteration);
    TH1D *refoldedRatio = dynamic_cast<TH1D *>(reconv->Clone(refoldedRatioName.c_str()));
    if(!refoldedRatio)
    {
      std::cerr << "Failed to build reconvolved ratio histogram for iteration " << iteration << std::endl;
      delete unfoldedRatio;
      delete unfoldedAbs;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }

    const std::string refoldedAbsName = "hRefoldedAbs_iter" + std::to_string(iteration);
    TH1D *refoldedAbs = dynamic_cast<TH1D *>(reconv->Clone(refoldedAbsName.c_str()));
    if(!refoldedAbs)
    {
      std::cerr << "Failed to build reconvolved absolute histogram for iteration " << iteration << std::endl;
      delete unfoldedRatio;
      delete unfoldedAbs;
      delete refoldedRatio;
      if(reconvComputed)
        delete reconvComputed;
      continue;
    }

    unfoldedAbs->SetDirectory(nullptr);
    refoldedAbs->SetDirectory(nullptr);
    unfoldedRatio->SetDirectory(nullptr);
    refoldedRatio->SetDirectory(nullptr);
    unfoldedAbs->SetTitle(("Deconvolved absolute overlaps;" + std::string(measured->GetXaxis()->GetTitle()) + ";Entries").c_str());
    refoldedAbs->SetTitle(("Reconvolved absolute overlaps;" + std::string(measured->GetXaxis()->GetTitle()) + ";Entries").c_str());
    unfoldedRatio->SetTitle(("Deconvolved/iter" + std::to_string(unfoldedRatioDenominatorIteration) + " ratio;" + std::string(measured->GetXaxis()->GetTitle()) + ";Deconvolved / Deconvolved_{den}").c_str());
    refoldedRatio->SetTitle(("Reconvolved/Measured ratio;" + std::string(measured->GetXaxis()->GetTitle()) + ";Reconvolved / Measured").c_str());
    for(int bin = 1; bin <= unfoldedRatio->GetNbinsX(); ++bin)
    {
      const double u = deconv->GetBinContent(bin);
      const double du = unfoldedRatioDenominator->GetBinContent(bin);
      const double f = reconv->GetBinContent(bin);
      const double m = measured->GetBinContent(bin);
      const double eu = deconv->GetBinError(bin);
      const double edu = unfoldedRatioDenominator->GetBinError(bin);
      const double ef = reconv->GetBinError(bin);
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
    storedDamping.push_back(damping);
    unfoldedAbsHists.push_back(unfoldedAbs);
    refoldedAbsHists.push_back(refoldedAbs);
    unfoldedRatioHists.push_back(unfoldedRatio);
    refoldedRatioHists.push_back(refoldedRatio);
    if(reconvComputed)
      delete reconvComputed;
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
      measuredLocal->SetMinimum(std::max(1e-5, minPos * 0.5));
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
      if(std::isfinite(storedDamping[i]))
      {
        char dbuf[64];
        std::snprintf(dbuf, sizeof(dbuf), " (#alpha=%.2e)", storedDamping[i]);
        label += dbuf;
      }
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

  TCanvas *cUnfoldedComparison = makeComparisonCanvas("cUnfoldedComparison", "Deconvolved Comparison", unfoldedAbsHists, unfoldedRatioHists, "Deconvolved", false);
  TCanvas *cRefoldedComparison = makeComparisonCanvas("cRefoldedComparison", "Reconvolved Comparison", refoldedAbsHists, refoldedRatioHists, "Reconvolved", false);
  TCanvas *cUnfoldedComparisonX0to25 = makeComparisonCanvas("cUnfoldedComparison_X0to25", "Deconvolved Comparison (x=[0,25])", unfoldedAbsHists, unfoldedRatioHists, "Deconvolved", true);
  TCanvas *cRefoldedComparisonX0to25 = makeComparisonCanvas("cRefoldedComparison_X0to25", "Reconvolved Comparison (x=[0,25])", refoldedAbsHists, refoldedRatioHists, "Reconvolved", true);

  TCanvas *cChi2NDF = new TCanvas("cChi2NDF", "Chi2/NDF vs Iteration", 800, 600);
  hChi2NDF->SetStats(0);
  hChi2NDF->SetLineColor(kBlue + 1);
  hChi2NDF->SetMarkerColor(kBlue + 1);
  hChi2NDF->SetMarkerStyle(20);
  hChi2NDF->SetLineWidth(2);
  hChi2NDF->Draw("HIST");
  hChi2NDF->Draw("P SAME");

  const std::string outputBase = stripExtension(outputFileName);
  if(cUnfoldedComparison)
    cUnfoldedComparison->SaveAs((outputBase + "_deconvolved_comparison.pdf").c_str());
  if(cRefoldedComparison)
    cRefoldedComparison->SaveAs((outputBase + "_reconvolved_comparison.pdf").c_str());
  if(cUnfoldedComparisonX0to25)
    cUnfoldedComparisonX0to25->SaveAs((outputBase + "_deconvolved_comparison_x0_25.pdf").c_str());
  if(cRefoldedComparisonX0to25)
    cRefoldedComparisonX0to25->SaveAs((outputBase + "_reconvolved_comparison_x0_25.pdf").c_str());
  if(cChi2NDF)
    cChi2NDF->SaveAs((outputBase + "_chi2ndf.pdf").c_str());

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
