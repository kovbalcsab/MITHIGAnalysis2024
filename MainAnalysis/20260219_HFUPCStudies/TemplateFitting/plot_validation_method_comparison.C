#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TLine.h>
#include <TGraphErrors.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <sstream>

#include "CommandLine.h"
#include "RootIOUtils.h"

static std::string replaceTilde(std::string s)
{
  size_t pos = s.find('~');
  while(pos != std::string::npos)
  {
    s[pos] = ' ';
    pos = s.find('~');
  }
  return s;
}

static void styleUpAxis(TH1 *hist)
{
  hist->GetXaxis()->SetTitleSize(0.045);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->ChangeLabel(1, 0, 0);
}

static void styleDownAxis(TH1 *hist)
{
  hist->GetXaxis()->SetTitleSize(0.05 / 0.4);
  hist->GetXaxis()->SetLabelSize(0.045 / 0.4);
  hist->GetYaxis()->SetTitleSize(0.05 / 0.4);
  hist->GetYaxis()->SetLabelSize(0.045 / 0.4);
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetYaxis()->SetTitleOffset(0.3);
  hist->GetYaxis()->ChangeLabel(-1, 0, 0);
  hist->GetYaxis()->SetNdivisions(505, true);
}

static TH1D *loadHistogram(TFile *file, const std::string &name, const std::string &cloneName)
{
  if(file == nullptr)
    return nullptr;

  TH1 *h = dynamic_cast<TH1 *>(file->Get(name.c_str()));
  if(h == nullptr)
    return nullptr;

  TH1D *h1 = dynamic_cast<TH1D *>(h->Clone(cloneName.c_str()));
  if(h1 == nullptr)
    return nullptr;

  h1->SetDirectory(nullptr);
  return h1;
}

static std::string joinStrings(const std::vector<std::string> &items, const std::string &sep)
{
  std::ostringstream out;
  for(size_t i = 0; i < items.size(); ++i)
  {
    if(i > 0)
      out << sep;
    out << items[i];
  }
  return out.str();
}

static void normalizeHistogram(TH1D *hist)
{
  if(hist == nullptr)
    return;
  const double integral = hist->Integral("width");
  if(integral > 0.0)
    hist->Scale(1.0 / integral);
}

static bool almostEqual(double a, double b, double relTol = 1e-10, double absTol = 1e-14)
{
  if(std::fabs(a - b) <= absTol)
    return true;

  const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
  return std::fabs(a - b) <= relTol * scale;
}

static bool histogramsCompatibleUpToNormalization(const TH1D *a,
                                                  const TH1D *b,
                                                  double relTol = 1e-8,
                                                  double absTol = 1e-12)
{
  if(a == nullptr || b == nullptr)
    return false;
  if(a->GetNbinsX() != b->GetNbinsX())
    return false;
  if(!almostEqual(a->GetXaxis()->GetXmin(), b->GetXaxis()->GetXmin()))
    return false;
  if(!almostEqual(a->GetXaxis()->GetXmax(), b->GetXaxis()->GetXmax()))
    return false;

  for(int edge = 1; edge <= a->GetNbinsX() + 1; ++edge)
  {
    if(!almostEqual(a->GetXaxis()->GetBinLowEdge(edge), b->GetXaxis()->GetBinLowEdge(edge)))
      return false;
  }

  const double intA = a->Integral("width");
  const double intB = b->Integral("width");
  if(intA <= 0.0 || intB <= 0.0)
    return false;

  for(int bin = 0; bin <= a->GetNbinsX() + 1; ++bin)
  {
    const double aContent = a->GetBinContent(bin) / intA;
    const double bContent = b->GetBinContent(bin) / intB;
    if(!almostEqual(aContent, bContent, relTol, absTol))
      return false;

    const double aError = a->GetBinError(bin) / intA;
    const double bError = b->GetBinError(bin) / intB;
    if(!almostEqual(aError, bError, relTol, absTol))
      return false;
  }
  return true;
}

static TGraphErrors *makeRatioGraph(const TH1D *numerator,
                                    const TH1D *denominator,
                                    const std::string &name,
                                    double xMinUse,
                                    double xMaxUse,
                                    bool includeErrors)
{
  if(numerator == nullptr || denominator == nullptr)
    return nullptr;

  const int nBins = numerator->GetNbinsX();
  const double denXMin = denominator->GetXaxis()->GetXmin();
  const double denXMax = denominator->GetXaxis()->GetXmax();

  auto *graph = new TGraphErrors();
  graph->SetName(name.c_str());
  graph->SetTitle("");

  int point = 0;
  for(int i = 1; i <= nBins; ++i)
  {
    const double x = numerator->GetBinCenter(i);
    if(x < xMinUse || x > xMaxUse)
      continue;
    if(x < denXMin || x > denXMax)
      continue;

    const double ex = 0.5 * numerator->GetBinWidth(i);
    const double num = numerator->GetBinContent(i);
    const double enumv = numerator->GetBinError(i);

    int denBin = denominator->GetXaxis()->FindFixBin(x);
    if(denBin < 1)
      denBin = 1;
    if(denBin > denominator->GetNbinsX())
      denBin = denominator->GetNbinsX();
    const double den = denominator->GetBinContent(denBin);
    const double eden = denominator->GetBinError(denBin);
    if(den <= 0.0)
      continue;

    const double ratio = num / den;
    double eratio = 0.0;
    if(includeErrors)
    {
      double relNum2 = 0.0;
      double relDen2 = 0.0;
      if(num > 0.0)
        relNum2 = (enumv / num) * (enumv / num);
      relDen2 = (eden / den) * (eden / den);
      eratio = ratio * std::sqrt(relNum2 + relDen2);
    }

    graph->SetPoint(point, x, ratio);
    graph->SetPointError(point, ex, eratio);
    point++;
  }

  return graph;
}

static void updateMinMaxVisible(const TH1D *hist,
                                double xMinUse,
                                double xMaxUse,
                                double &minY,
                                double &maxY,
                                bool includeErrors,
                                bool positiveOnly)
{
  if(hist == nullptr)
    return;

  for(int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    const double x = hist->GetBinCenter(i);
    if(x < xMinUse || x > xMaxUse)
      continue;

    const double y = hist->GetBinContent(i);
    const double ey = includeErrors ? hist->GetBinError(i) : 0.0;
    const double low = y - ey;
    const double high = y + ey;

    if(positiveOnly)
    {
      if(high <= 0.0)
        continue;
      if(low > 0.0)
        minY = std::min(minY, low);
      else
        minY = std::min(minY, high * 0.2);
    }
    else
      minY = std::min(minY, low);

    maxY = std::max(maxY, high);
  }
}

struct MethodInput
{
  TH1D *unfoldedNorm = nullptr;
  TH1D *closureRatio = nullptr;
  TH1D *refoldedNorm = nullptr;
  TH1D *chi2Profile = nullptr;
  Color_t color = kBlack;
  std::string label;
  int iteration = -1;
};

int main(int argc, char **argv)
{
  CommandLine CL(argc, argv);

  std::vector<std::string> inputFileNames = CL.GetStringVector("InputFileNames", "");
  std::vector<std::string> labels = CL.GetStringVector("Labels", "");
  std::vector<int> colorsIn = CL.GetIntVector("Colors", "");
  std::vector<int> selectedIterations = CL.GetIntVector("SelectedIterations", "");

  std::string outputPrefix = CL.Get("OutputPrefix", "TemplateFitting/output_unfolding/validation_method_comparison_selected");
  std::string unfoldedOutput = CL.Get("OutputUnfoldedFileName", outputPrefix + "_unfolded.pdf");
  std::string closureRatioOutput = CL.Get("OutputClosureRatioFileName", outputPrefix + "_closure_ratio.pdf");
  std::string closureAbsOutput = CL.Get("OutputClosureAbsoluteFileName", outputPrefix + "_closure_absolute.pdf");
  std::string closureAbsZoomOutput = CL.Get("OutputClosureAbsoluteZoomFileName", outputPrefix + "_closure_absolute_x0_25.pdf");
  std::string chi2Output = CL.Get("OutputChi2FileName", outputPrefix + "_chi2ndf.pdf");

  std::string unfoldedPrefix = CL.Get("UnfoldedHistPrefix", "hUnfoldedAbs_iter");
  std::string closurePrefix = CL.Get("ClosureHistPrefix", "hRefoldedToMeasuredRatio_iter");
  std::string refoldedAbsPrefix = CL.Get("RefoldedAbsHistPrefix", "hRefoldedAbs_iter");
  std::string measuredHistName = CL.Get("MeasuredHistName", "hMeasuredForOverlay");
  std::string chi2HistName = CL.Get("Chi2HistName", "hChi2NDFVsIter");

  std::string unfoldedXTitle = replaceTilde(CL.Get("UnfoldedXTitle", "HF~energy~[GeV]"));
  std::string unfoldedYTitle = replaceTilde(CL.Get("UnfoldedYTitle", "1/N~dN/dE"));
  std::string closureXTitle = replaceTilde(CL.Get("ClosureXTitle", "HF~energy~[GeV]"));
  std::string closureYTitle = replaceTilde(CL.Get("ClosureYTitle", "Refolded/Measured"));
  std::string closureAbsXTitle = replaceTilde(CL.Get("ClosureAbsXTitle", "HF~energy~[GeV]"));
  std::string closureAbsYTitle = replaceTilde(CL.Get("ClosureAbsYTitle", "1/N~dN/dE"));
  std::string chi2XTitle = replaceTilde(CL.Get("Chi2XTitle", "Regularization~iteration"));
  std::string chi2YTitle = replaceTilde(CL.Get("Chi2YTitle", "#chi^{2}/NDF"));

  std::string unfoldedTitle = replaceTilde(CL.Get("UnfoldedTitle", "Selected~iteration~unfolded/deconvolved~comparison"));
  std::string closureTitle = replaceTilde(CL.Get("ClosureTitle", "Selected~iteration~closure~ratio~comparison"));
  std::string closureAbsTitle = replaceTilde(CL.Get("ClosureAbsTitle", "Selected~iteration~measured~vs~refolded/reconvolved"));
  std::string chi2Title = replaceTilde(CL.Get("Chi2Title", "Method~comparison~of~#chi^{2}/NDF~scan"));

  bool logYUnfolded = CL.GetBool("LogYUnfolded", true);
  bool doRatioToFirst = CL.GetBool("DoRatioToFirst", true);
  bool requireAllMethods = CL.GetBool("RequireAllMethods", true);
  bool normalizeAbsolute = CL.GetBool("NormalizeAbsolute", true);
  bool chi2RemoveErrorBars = CL.GetBool("Chi2RemoveErrorBars", true);
  double ratioYMin = CL.GetDouble("RatioYMin", 0.0);
  double ratioYMax = CL.GetDouble("RatioYMax", 0.0);

  double xMin = CL.GetDouble("XMin", 0.0);
  double xMax = CL.GetDouble("XMax", 0.0);
  double zoomXMin = CL.GetDouble("ZoomXMin", 0.0);
  double zoomXMax = CL.GetDouble("ZoomXMax", 25.0);
  double yMinUnfolded = CL.GetDouble("YMinUnfolded", 0.0);
  double yMaxUnfolded = CL.GetDouble("YMaxUnfolded", 0.0);
  double yMinClosure = CL.GetDouble("YMinClosure", 0.0);
  double yMaxClosure = CL.GetDouble("YMaxClosure", 0.0);
  double yMinClosureAbs = CL.GetDouble("YMinClosureAbs", 0.0);
  double yMaxClosureAbs = CL.GetDouble("YMaxClosureAbs", 0.0);
  double yMinChi2 = CL.GetDouble("YMinChi2", 0.0);
  double yMaxChi2 = CL.GetDouble("YMaxChi2", 0.0);

  if(inputFileNames.empty() || labels.empty() || selectedIterations.empty())
  {
    std::cerr << "InputFileNames, Labels and SelectedIterations must be provided and non-empty." << std::endl;
    return -1;
  }
  if(inputFileNames.size() != labels.size() || inputFileNames.size() != selectedIterations.size())
  {
    std::cerr << "InputFileNames, Labels and SelectedIterations must have the same size." << std::endl;
    return -1;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const std::vector<Color_t> defaultColors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2, kOrange + 7};

  std::vector<TFile *> inputFiles;
  std::vector<MethodInput> methods;
  TH1D *commonMeasuredNorm = nullptr;

  for(size_t i = 0; i < inputFileNames.size(); ++i)
  {
    TFile *file = RootIOUtils::OpenFileOrNull(inputFileNames[i], "READ", "input file");
    if(file == nullptr)
    {
      if(requireAllMethods)
      {
        for(TFile *f : inputFiles)
          RootIOUtils::CloseAndDeleteFile(f);
        delete commonMeasuredNorm;
        return -1;
      }
      continue;
    }
    inputFiles.push_back(file);

    const int iter = selectedIterations[i];
    const std::string unfoldedName = unfoldedPrefix + std::to_string(iter);
    const std::string closureName = closurePrefix + std::to_string(iter);
    const std::string refoldedAbsName = refoldedAbsPrefix + std::to_string(iter);

    MethodInput method;
    method.unfoldedNorm = loadHistogram(file, unfoldedName, Form("hUnfoldedNorm_%zu", i));
    method.closureRatio = loadHistogram(file, closureName, Form("hClosureRatio_%zu", i));
    method.refoldedNorm = loadHistogram(file, refoldedAbsName, Form("hRefoldedNorm_%zu", i));
    TH1D *measuredCandidate = loadHistogram(file, measuredHistName, Form("hMeasuredNorm_%zu", i));
    method.chi2Profile = loadHistogram(file, chi2HistName, Form("hChi2_%zu", i));
    method.iteration = iter;
    method.label = replaceTilde(labels[i]);
    method.color = (i < colorsIn.size()) ? static_cast<Color_t>(colorsIn[i]) : defaultColors[i % defaultColors.size()];

    if(method.unfoldedNorm == nullptr || method.closureRatio == nullptr || method.refoldedNorm == nullptr || measuredCandidate == nullptr || method.chi2Profile == nullptr)
    {
      std::cerr << "Missing required histograms in " << inputFileNames[i] << " for selected iteration " << iter << std::endl;
      delete method.unfoldedNorm;
      delete method.closureRatio;
      delete method.refoldedNorm;
      delete measuredCandidate;
      delete method.chi2Profile;
      if(requireAllMethods)
      {
        for(TFile *f : inputFiles)
          RootIOUtils::CloseAndDeleteFile(f);
        return -1;
      }
      continue;
    }

    if(commonMeasuredNorm == nullptr)
    {
      commonMeasuredNorm = dynamic_cast<TH1D *>(measuredCandidate->Clone("hMeasuredCommonNorm"));
      if(commonMeasuredNorm == nullptr)
      {
        std::cerr << "Failed to clone common measured histogram from " << inputFileNames[i] << std::endl;
        delete method.unfoldedNorm;
        delete method.closureRatio;
        delete method.refoldedNorm;
        delete measuredCandidate;
        delete method.chi2Profile;
        for(TFile *f : inputFiles)
          RootIOUtils::CloseAndDeleteFile(f);
        delete commonMeasuredNorm;
        return -1;
      }
      commonMeasuredNorm->SetDirectory(nullptr);
    }
    else
    {
      if(!histogramsCompatibleUpToNormalization(commonMeasuredNorm, measuredCandidate))
      {
        std::cerr << "Measured histogram mismatch detected in " << inputFileNames[i]
                  << ". All input files must contain a measured histogram compatible up to a global normalization factor ("
                  << measuredHistName << ")." << std::endl;
        delete method.unfoldedNorm;
        delete method.closureRatio;
        delete method.refoldedNorm;
        delete measuredCandidate;
        delete method.chi2Profile;
        for(TFile *f : inputFiles)
          RootIOUtils::CloseAndDeleteFile(f);
        delete commonMeasuredNorm;
        return -1;
      }
    }
    delete measuredCandidate;

    if(normalizeAbsolute)
    {
      normalizeHistogram(method.unfoldedNorm);
      normalizeHistogram(method.refoldedNorm);
    }

    if(chi2RemoveErrorBars)
      for(int b = 1; b <= method.chi2Profile->GetNbinsX(); ++b)
        method.chi2Profile->SetBinError(b, 0.0);

    method.unfoldedNorm->SetLineWidth(2);
    method.closureRatio->SetLineWidth(2);
    method.refoldedNorm->SetLineWidth(2);
    method.chi2Profile->SetLineWidth(2);

    method.unfoldedNorm->SetLineColor(method.color);
    method.unfoldedNorm->SetMarkerColor(method.color);
    method.unfoldedNorm->SetMarkerStyle(20 + static_cast<int>(i % 5));
    method.closureRatio->SetLineColor(method.color);
    method.closureRatio->SetMarkerColor(method.color);
    method.closureRatio->SetMarkerStyle(20 + static_cast<int>(i % 5));
    method.refoldedNorm->SetLineColor(method.color);
    method.refoldedNorm->SetMarkerColor(method.color);
    method.refoldedNorm->SetMarkerStyle(24 + static_cast<int>(i % 5));
    method.chi2Profile->SetLineColor(method.color);
    method.chi2Profile->SetMarkerColor(method.color);
    method.chi2Profile->SetMarkerStyle(20 + static_cast<int>(i % 5));

    methods.push_back(method);
  }

  if(methods.empty())
  {
    std::cerr << "No valid method inputs loaded." << std::endl;
    for(TFile *f : inputFiles)
      RootIOUtils::CloseAndDeleteFile(f);
    delete commonMeasuredNorm;
    return -1;
  }
  if(commonMeasuredNorm == nullptr)
  {
    std::cerr << "Failed to load common measured histogram." << std::endl;
    for(TFile *f : inputFiles)
      RootIOUtils::CloseAndDeleteFile(f);
    return -1;
  }
  if(normalizeAbsolute)
    normalizeHistogram(commonMeasuredNorm);
  commonMeasuredNorm->SetLineColor(kGray + 2);
  commonMeasuredNorm->SetMarkerColor(kGray + 2);
  commonMeasuredNorm->SetMarkerStyle(20);
  commonMeasuredNorm->SetLineWidth(2);

  if(methods.size() < 2)
    doRatioToFirst = false;

  auto inferXRange = [&](const std::vector<TH1D *> &hists, double &xLo, double &xHi)
  {
    xLo = std::numeric_limits<double>::infinity();
    xHi = -std::numeric_limits<double>::infinity();
    for(const TH1D *h : hists)
    {
      if(h == nullptr)
        continue;
      xLo = std::min(xLo, h->GetXaxis()->GetXmin());
      xHi = std::max(xHi, h->GetXaxis()->GetXmax());
    }
    if(!std::isfinite(xLo) || !std::isfinite(xHi) || xLo >= xHi)
    {
      xLo = 0.0;
      xHi = 1.0;
    }
  };

  auto drawGenericComparison = [&](const std::vector<TH1D *> &hists,
                                   const std::vector<std::string> &legendLabels,
                                   const std::vector<Color_t> &plotColors,
                                   const std::string &canvasName,
                                   const std::string &title,
                                   const std::string &xTitle,
                                   const std::string &yTitle,
                                   const std::string &outputName,
                                   bool useLogY,
                                   bool withRatio,
                                   bool withErrors,
                                   double yMinUser,
                                   double yMaxUser,
                                   double xMinUser,
                                   double xMaxUser,
                                   bool forceAbsoluteYMin,
                                   bool dynamicRatioRange)
  {
    TCanvas *canvas = new TCanvas(canvasName.c_str(), title.c_str(), 900, 900);
    TPad *padUp = nullptr;
    TPad *padDown = nullptr;
    if(withRatio)
    {
      padUp = new TPad((canvasName + "_up").c_str(), "", 0.0, 0.30, 1.0, 1.0);
      padDown = new TPad((canvasName + "_down").c_str(), "", 0.0, 0.05, 1.0, 0.30);
      padUp->SetBottomMargin(0.0);
      padDown->SetTopMargin(0.0);
      padDown->SetBottomMargin(0.40);
      padUp->SetLeftMargin(0.13);
      padDown->SetLeftMargin(0.13);
      padUp->SetRightMargin(0.04);
      padDown->SetRightMargin(0.04);
      canvas->cd();
      padUp->Draw();
      padDown->Draw();
      padUp->cd();
      if(useLogY)
        padUp->SetLogy();
    }
    else if(useLogY)
      canvas->SetLogy();

    double xMinAuto = 0.0;
    double xMaxAuto = 1.0;
    inferXRange(hists, xMinAuto, xMaxAuto);
    const double xMinUse = (xMinUser < xMaxUser) ? xMinUser : xMinAuto;
    const double xMaxUse = (xMinUser < xMaxUser) ? xMaxUser : xMaxAuto;

    double yMinAuto = std::numeric_limits<double>::infinity();
    double yMaxAuto = -std::numeric_limits<double>::infinity();
    for(const TH1D *h : hists)
      updateMinMaxVisible(h, xMinUse, xMaxUse, yMinAuto, yMaxAuto, withErrors, useLogY);

    double yMinUse = yMinUser;
    double yMaxUse = yMaxUser;
    if(yMinUse >= yMaxUse)
    {
      if(useLogY)
      {
        const double minPos = (std::isfinite(yMinAuto) && yMinAuto > 0.0) ? yMinAuto : 1e-7;
        const double maxPos = (std::isfinite(yMaxAuto) && yMaxAuto > minPos) ? yMaxAuto : (minPos * 100.0);
        yMinUse = std::max(1e-7, minPos * 0.5);
        yMaxUse = maxPos * 2.2;
      }
      else
      {
        const double minVal = std::isfinite(yMinAuto) ? yMinAuto : 0.0;
        const double maxVal = (std::isfinite(yMaxAuto) && yMaxAuto > minVal) ? yMaxAuto : (minVal + 1.0);
        yMinUse = minVal - 0.10 * (maxVal - minVal);
        yMaxUse = maxVal + 0.20 * (maxVal - minVal);
      }
    }
    if(forceAbsoluteYMin)
    {
      if(yMaxUse <= 1e-5)
        yMaxUse = 10e-5;
      yMinUse = 1e-5;
    }

    TLegend *legend = new TLegend(0.52, 0.62, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.038);

    for(size_t i = 0; i < hists.size(); ++i)
    {
      TH1D *h = hists[i];
      const double histXMin = h->GetXaxis()->GetXmin();
      const double histXMax = h->GetXaxis()->GetXmax();
      const double histRangeMin = std::max(xMinUse, histXMin);
      const double histRangeMax = std::min(xMaxUse, histXMax);

      h->SetTitle(title.c_str());
      h->GetXaxis()->SetTitle(xTitle.c_str());
      h->GetYaxis()->SetTitle(yTitle.c_str());
      if(histRangeMin < histRangeMax)
        h->GetXaxis()->SetRangeUser(histRangeMin, histRangeMax);
      h->GetYaxis()->SetRangeUser(yMinUse, yMaxUse);
      h->SetLineColor(plotColors[i]);
      h->SetMarkerColor(plotColors[i]);

      if(withRatio)
      {
        styleUpAxis(h);
        h->GetXaxis()->SetLabelSize(0.0);
        h->GetXaxis()->SetTitleSize(0.0);
      }

      const char *drawOpt = withErrors ? "EP" : "P";
      const std::string opt = (i == 0) ? drawOpt : std::string(drawOpt) + " SAME";
      h->Draw(opt.c_str());
      legend->AddEntry(h, legendLabels[i].c_str(), "lp");
    }
    legend->Draw();

    if(withRatio)
    {
      canvas->cd();
      padDown->cd();
      TH1D *axisHist = dynamic_cast<TH1D *>(hists.front()->Clone((canvasName + "_ratioAxis").c_str()));
      axisHist->Reset("ICES");
      axisHist->SetDirectory(nullptr);
      axisHist->SetTitle("");
      axisHist->GetXaxis()->SetTitle(xTitle.c_str());
      std::string ratioTitle = "Ratio to first";
      if(!legendLabels.empty())
        ratioTitle = "Ratio to " + legendLabels.front();
      axisHist->GetYaxis()->SetTitle(ratioTitle.c_str());
      axisHist->GetXaxis()->SetRangeUser(xMinUse, xMaxUse);
      double ratioMinUse = ratioYMin;
      double ratioMaxUse = ratioYMax;
      if(!(ratioMinUse < ratioMaxUse))
      {
        ratioMinUse = 0.0;
        ratioMaxUse = 2.0;
      }
      if(dynamicRatioRange && !(ratioYMin < ratioYMax))
      {
        double ratioMinAuto = std::numeric_limits<double>::infinity();
        double ratioMaxAuto = -std::numeric_limits<double>::infinity();
        for(size_t i = 1; i < hists.size(); ++i)
        {
          std::unique_ptr<TGraphErrors> ratioProbe(
            makeRatioGraph(hists[i], hists.front(), Form("%s_ratio_probe_%zu", canvasName.c_str(), i), xMinUse, xMaxUse, withErrors));
          if(!ratioProbe || ratioProbe->GetN() == 0)
            continue;
          for(int p = 0; p < ratioProbe->GetN(); ++p)
          {
            double x = 0.0;
            double y = 0.0;
            ratioProbe->GetPoint(p, x, y);
            const double ey = withErrors ? ratioProbe->GetErrorY(p) : 0.0;
            ratioMinAuto = std::min(ratioMinAuto, y - ey);
            ratioMaxAuto = std::max(ratioMaxAuto, y + ey);
          }
        }
        if(std::isfinite(ratioMinAuto) && std::isfinite(ratioMaxAuto) && ratioMinAuto < ratioMaxAuto)
        {
          const double center = 1.0;
          const double span = std::max(std::fabs(ratioMaxAuto - center), std::fabs(center - ratioMinAuto));
          const double paddedSpan = std::max(0.15, span * 1.20);
          ratioMinUse = center - paddedSpan;
          ratioMaxUse = center + paddedSpan;
        }
      }
      axisHist->GetYaxis()->SetRangeUser(ratioMinUse, ratioMaxUse);
      styleDownAxis(axisHist);
      axisHist->Draw("AXIS");

      TLine *lineOne = new TLine(xMinUse, 1.0, xMaxUse, 1.0);
      lineOne->SetLineStyle(2);
      lineOne->SetLineColor(kGray + 2);
      lineOne->Draw("SAME");

      std::vector<std::unique_ptr<TGraphErrors>> ratioGraphs;
      for(size_t i = 1; i < hists.size(); ++i)
      {
        std::unique_ptr<TGraphErrors> ratio(makeRatioGraph(hists[i], hists.front(), Form("%s_ratio_%zu", canvasName.c_str(), i), xMinUse, xMaxUse, withErrors));
        if(!ratio || ratio->GetN() == 0)
          continue;
        ratio->SetLineColor(plotColors[i]);
        ratio->SetMarkerColor(plotColors[i]);
        ratio->SetMarkerStyle(20 + static_cast<int>(i % 5));
        ratio->SetLineWidth(2);
        ratio->Draw((withErrors ? "EP SAME" : "P SAME"));
        ratioGraphs.push_back(std::move(ratio));
      }
    }

    canvas->SaveAs(outputName.c_str());
    delete canvas;
  };

  std::vector<TH1D *> unfoldedHists;
  std::vector<TH1D *> closureRatioHists;
  std::vector<TH1D *> chi2Hists;
  std::vector<Color_t> methodColors;
  std::vector<std::string> methodLabels;
  std::vector<std::string> methodLabelsWithIter;
  for(const MethodInput &m : methods)
  {
    unfoldedHists.push_back(m.unfoldedNorm);
    closureRatioHists.push_back(m.closureRatio);
    chi2Hists.push_back(m.chi2Profile);
    methodColors.push_back(m.color);
    methodLabels.push_back(m.label);
    methodLabelsWithIter.push_back(m.label + Form(" (iter %d)", m.iteration));
  }

  if(normalizeAbsolute)
  {
    unfoldedTitle += " [normalized]";
    closureAbsTitle += " [normalized]";
  }
  if(chi2RemoveErrorBars)
    chi2Title += " [errors hidden]";

  drawGenericComparison(unfoldedHists,
                        methodLabelsWithIter,
                        methodColors,
                        "cValidationMethodUnfoldedComparison",
                        unfoldedTitle,
                        unfoldedXTitle,
                        unfoldedYTitle,
                        unfoldedOutput,
                        logYUnfolded,
                        doRatioToFirst,
                        true,
                        yMinUnfolded,
                        yMaxUnfolded,
                        xMin,
                        xMax,
                        true,
                        true);

  drawGenericComparison(closureRatioHists,
                        methodLabelsWithIter,
                        methodColors,
                        "cValidationMethodClosureRatioComparison",
                        closureTitle,
                        closureXTitle,
                        closureYTitle,
                        closureRatioOutput,
                        false,
                        doRatioToFirst,
                        true,
                        yMinClosure,
                        yMaxClosure,
                        xMin,
                        xMax,
                        false,
                        true);

  {
    std::vector<TH1D *> closureAbs;
    std::vector<Color_t> closureAbsColors;
    std::vector<std::string> closureAbsLabels;
    TH1D *meas = dynamic_cast<TH1D *>(commonMeasuredNorm->Clone("hMeasuredAbs_common"));
    if(meas)
    {
      meas->SetDirectory(nullptr);
      meas->SetLineColor(kGray + 2);
      meas->SetMarkerColor(kGray + 2);
      meas->SetMarkerStyle(20);
      meas->SetLineWidth(2);
      closureAbs.push_back(meas);
      closureAbsColors.push_back(kGray + 2);
      closureAbsLabels.push_back("Measured");
    }
    for(size_t i = 0; i < methods.size(); ++i)
    {
      TH1D *ref = dynamic_cast<TH1D *>(methods[i].refoldedNorm->Clone(Form("hRefoldedAbs_%zu", i)));
      if(!ref)
      {
        delete ref;
        continue;
      }
      ref->SetDirectory(nullptr);
      ref->SetLineColor(methods[i].color);
      ref->SetMarkerColor(methods[i].color);
      ref->SetMarkerStyle(24 + static_cast<int>(i % 5));
      ref->SetLineWidth(2);

      closureAbs.push_back(ref);
      closureAbsColors.push_back(methods[i].color);
      closureAbsLabels.push_back(methods[i].label + " refolded");
    }

    drawGenericComparison(closureAbs,
                          closureAbsLabels,
                          closureAbsColors,
                          "cValidationMethodClosureAbsoluteComparison",
                          closureAbsTitle,
                          closureAbsXTitle,
                          closureAbsYTitle,
                          closureAbsOutput,
                          true,
                          false,
                          true,
                           yMinClosureAbs,
                           yMaxClosureAbs,
                           xMin,
                           xMax,
                           true,
                           false);

    drawGenericComparison(closureAbs,
                          closureAbsLabels,
                          closureAbsColors,
                          "cValidationMethodClosureAbsoluteComparisonZoom",
                          closureAbsTitle + " (x < 25 GeV)",
                          closureAbsXTitle,
                          closureAbsYTitle,
                          closureAbsZoomOutput,
                          true,
                          false,
                          true,
                           yMinClosureAbs,
                           yMaxClosureAbs,
                           zoomXMin,
                           zoomXMax,
                           true,
                           false);

    for(TH1D *h : closureAbs)
      delete h;
  }

  drawGenericComparison(chi2Hists,
                        methodLabels,
                        methodColors,
                        "cValidationMethodChi2Comparison",
                        chi2Title,
                        chi2XTitle,
                        chi2YTitle,
                        chi2Output,
                        false,
                        doRatioToFirst,
                        false,
                        yMinChi2,
                        yMaxChi2,
                        xMin,
                        xMax,
                        false,
                        true);

  for(TFile *f : inputFiles)
    RootIOUtils::CloseAndDeleteFile(f);
  delete commonMeasuredNorm;

  std::cout << "Wrote selected-iteration comparison plots:" << std::endl;
  std::cout << "  " << unfoldedOutput << std::endl;
  std::cout << "  " << closureRatioOutput << std::endl;
  std::cout << "  " << closureAbsOutput << std::endl;
  std::cout << "  " << closureAbsZoomOutput << std::endl;
  std::cout << "  " << chi2Output << std::endl;
  std::cout << "Transforms: NormalizeAbsolute=" << (normalizeAbsolute ? "true" : "false")
            << ", Chi2RemoveErrorBars=" << (chi2RemoveErrorBars ? "true" : "false")
            << ", DoRatioToFirst=" << (doRatioToFirst ? "true" : "false");
  if(doRatioToFirst && !methods.empty())
    std::cout << " (baseline=" << methods.front().label << ")";
  std::cout << std::endl;
  std::cout << "Selected methods: " << joinStrings(methodLabelsWithIter, ", ") << std::endl;

  return 0;
}
