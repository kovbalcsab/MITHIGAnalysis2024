#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "RootIOUtils.h"

namespace
{
struct QuarterCurve
{
  int quarter = -1;
  TH1 *hist = nullptr;
  std::string fileName;
};

std::string replaceTilde(std::string text)
{
  std::size_t pos = text.find('~');
  while(pos != std::string::npos)
  {
    text[pos] = ' ';
    pos = text.find('~');
  }
  return text;
}

int inferQuarter(const std::string &fileName, int fallback)
{
  const std::string key = "_vs_q";
  const std::size_t pos = fileName.find(key);
  if(pos == std::string::npos || pos + key.size() >= fileName.size())
    return fallback;

  int quarter = 0;
  bool foundDigit = false;
  for(std::size_t i = pos + key.size(); i < fileName.size(); ++i)
  {
    const char c = fileName[i];
    if(!std::isdigit(static_cast<unsigned char>(c)))
      break;
    foundDigit = true;
    quarter = quarter * 10 + (c - '0');
  }

  if(!foundDigit)
    return fallback;
  return quarter;
}

bool sameBinning(const TH1 *a, const TH1 *b, double tolerance = 1e-12)
{
  if(a == nullptr || b == nullptr)
    return false;
  if(a->GetNbinsX() != b->GetNbinsX())
    return false;
  for(int edge = 1; edge <= a->GetNbinsX() + 1; ++edge)
  {
    const double edgeA = a->GetXaxis()->GetBinLowEdge(edge);
    const double edgeB = b->GetXaxis()->GetBinLowEdge(edge);
    if(std::fabs(edgeA - edgeB) > tolerance)
      return false;
  }
  return true;
}

Color_t quarterColor(int quarter, std::size_t index)
{
  if(quarter == 0) return kBlack;
  if(quarter == 1) return kRed + 1;
  if(quarter == 2) return kBlue + 1;
  if(quarter == 3) return kGreen + 2;

  const std::vector<Color_t> fallbackColors = {
    kMagenta + 1, kCyan + 2, kOrange + 7, kViolet + 1, kGray + 2
  };
  return fallbackColors[index % fallbackColors.size()];
}

bool loadQuarterCurves(const std::vector<std::string> &fileNames,
                       const std::string &histName,
                       const std::string &fallbackHistName,
                       const std::string &methodLabel,
                       bool requireAllQuarters,
                       std::vector<TFile *> &ownedFiles,
                       std::vector<QuarterCurve> &curves)
{
  curves.clear();

  for(std::size_t i = 0; i < fileNames.size(); ++i)
  {
    if(fileNames[i].empty())
      continue;

    TFile *file = RootIOUtils::OpenFileOrNull(fileNames[i], "READ", methodLabel + " input");
    if(file == nullptr)
    {
      if(requireAllQuarters)
        return false;
      continue;
    }

    TH1 *hist = dynamic_cast<TH1 *>(file->Get(histName.c_str()));
    if(hist == nullptr && !fallbackHistName.empty() && fallbackHistName != histName)
    {
      hist = dynamic_cast<TH1 *>(file->Get(fallbackHistName.c_str()));
      if(hist != nullptr)
      {
        std::cerr << "Warning: " << methodLabel << " file " << fileNames[i]
                  << " does not contain " << histName
                  << "; falling back to " << fallbackHistName << std::endl;
      }
    }
    if(hist == nullptr)
    {
      std::cerr << "Missing histogram " << histName
                << " in " << fileNames[i]
                << " for " << methodLabel << std::endl;
      RootIOUtils::CloseAndDeleteFile(file);
      if(requireAllQuarters)
        return false;
      continue;
    }

    QuarterCurve curve;
    curve.quarter = inferQuarter(fileNames[i], static_cast<int>(i));
    curve.hist = hist;
    curve.fileName = fileNames[i];
    curves.push_back(curve);
    ownedFiles.push_back(file);
  }

  if(curves.empty())
  {
    std::cerr << "No valid curves loaded for " << methodLabel << std::endl;
    return false;
  }

  std::sort(curves.begin(), curves.end(),
            [](const QuarterCurve &a, const QuarterCurve &b){ return a.quarter < b.quarter; });

  for(std::size_t i = 1; i < curves.size(); ++i)
  {
    if(!sameBinning(curves[0].hist, curves[i].hist))
    {
      std::cerr << "Incompatible iteration-axis binning in " << methodLabel
                << " between files:\n  " << curves[0].fileName
                << "\n  " << curves[i].fileName << std::endl;
      return false;
    }
  }

  if(requireAllQuarters)
  {
    std::set<int> loadedQuarters;
    for(const QuarterCurve &curve : curves)
      loadedQuarters.insert(curve.quarter);
    for(int quarter = 0; quarter < 4; ++quarter)
    {
      if(loadedQuarters.count(quarter) == 0)
      {
        std::cerr << "Missing quarter " << quarter << " in " << methodLabel << std::endl;
        return false;
      }
    }
  }

  return true;
}

void drawQuarterCurves(const std::vector<QuarterCurve> &curves,
                       const std::string &canvasName,
                       const std::string &title,
                       const std::string &xTitle,
                       const std::string &yTitle,
                       const std::string &outputFileName,
                       bool logY,
                       double yMinUser,
                       double yMaxUser)
{
  if(curves.empty())
    return;

  double yMin = std::numeric_limits<double>::infinity();
  double yMax = -std::numeric_limits<double>::infinity();
  for(const QuarterCurve &curve : curves)
  {
    TH1 *hist = curve.hist;
    if(hist == nullptr)
      continue;

    for(int bin = 1; bin <= hist->GetNbinsX(); ++bin)
    {
      const double value = hist->GetBinContent(bin);
      if(logY && value <= 0.0)
        continue;
      yMin = std::min(yMin, value);
      yMax = std::max(yMax, value);
    }
  }

  if(!std::isfinite(yMin) || !std::isfinite(yMax) || yMin >= yMax)
  {
    yMin = logY ? 1e-6 : 0.0;
    yMax = logY ? 1.0 : 1.0;
  }

  double yMinUse = yMinUser;
  double yMaxUse = yMaxUser;
  if(!(yMinUse < yMaxUse))
  {
    if(logY)
    {
      yMinUse = std::max(1e-8, yMin * 0.7);
      yMaxUse = yMax * 1.8;
    }
    else
    {
      yMinUse = yMin - 0.10 * (yMax - yMin);
      yMaxUse = yMax + 0.20 * (yMax - yMin);
    }
  }
  if(logY && yMinUse <= 0.0)
    yMinUse = std::max(1e-8, yMin * 0.7);

  TCanvas canvas(canvasName.c_str(), title.c_str(), 900, 700);
  if(logY)
    canvas.SetLogy();

  TLegend legend(0.56, 0.64, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.035);

  for(std::size_t i = 0; i < curves.size(); ++i)
  {
    TH1 *hist = curves[i].hist;
    if(hist == nullptr)
      continue;

    const Color_t color = quarterColor(curves[i].quarter, i);
    hist->SetStats(0);
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20 + static_cast<int>(i % 5));
    hist->SetLineWidth(2);

    if(i == 0)
    {
      hist->SetTitle((title + ";" + xTitle + ";" + yTitle).c_str());
      hist->GetYaxis()->SetRangeUser(yMinUse, yMaxUse);
      hist->Draw("HIST");
      hist->Draw("P SAME");
    }
    else
    {
      hist->Draw("HIST SAME");
      hist->Draw("P SAME");
    }

    legend.AddEntry(hist, ("Validation quarter " + std::to_string(curves[i].quarter)).c_str(), "lp");
  }

  std::vector<TLine> minLines;
  minLines.reserve(curves.size());
  for(std::size_t i = 0; i < curves.size(); ++i)
  {
    TH1 *hist = curves[i].hist;
    if(hist == nullptr)
      continue;

    int minBin = -1;
    double minValue = std::numeric_limits<double>::infinity();
    for(int bin = 1; bin <= hist->GetNbinsX(); ++bin)
    {
      const double value = hist->GetBinContent(bin);
      if(!std::isfinite(value))
        continue;
      if(logY && value <= 0.0)
        continue;
      if(value < minValue)
      {
        minValue = value;
        minBin = bin;
      }
    }
    if(minBin < 1)
      continue;

    const double xAtMin = hist->GetBinCenter(minBin);
    const Color_t color = quarterColor(curves[i].quarter, i);
    minLines.emplace_back(xAtMin, yMinUse, xAtMin, yMaxUse);
    TLine &line = minLines.back();
    line.SetLineColor(color);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.Draw("SAME");
  }

  legend.Draw();
  canvas.SaveAs(outputFileName.c_str());
}
} // namespace

int main(int argc, char **argv)
{
  CommandLine CL(argc, argv);

  const std::vector<std::string> bayesFiles = CL.GetStringVector("BayesFiles", "");
  const std::vector<std::string> svdFiles = CL.GetStringVector("SVDFiles", "");
  const std::vector<std::string> fftFiles = CL.GetStringVector("FFTFiles", "");

  if(bayesFiles.empty() || svdFiles.empty() || fftFiles.empty())
  {
    std::cerr << "Provide BayesFiles, SVDFiles and FFTFiles (comma-separated)." << std::endl;
    return -1;
  }

  const std::string chi2HistName = CL.Get("Chi2HistName", "hChi2VsIter");
  const bool makeMeasuredOnlyPlots = CL.GetBool("MakeMeasuredOnlyPlots", true);
  const std::string bayesMeasuredOnlyHistName = CL.Get("BayesMeasuredOnlyHistName", "hChi2NDFMeasuredOnlyVsIter");
  const std::string svdMeasuredOnlyHistName = CL.Get("SVDMeasuredOnlyHistName", "hChi2NDFMeasuredOnlyVsIter");
  const std::string fftMeasuredOnlyHistName = CL.Get("FFTMeasuredOnlyHistName", chi2HistName);
  const std::string outputPrefix = CL.Get("OutputPrefix", "TemplateFitting/output_unfolding/validation_quarter_chi2");
  const bool requireAllQuarters = CL.GetBool("RequireAllQuarters", true);
  const bool logY = CL.GetBool("LogY", false);
  const double yMinUser = CL.GetDouble("YMin", 0.0);
  const double yMaxUser = CL.GetDouble("YMax", 0.0);

  const std::string xTitle = replaceTilde(CL.Get("XTitle", "Iteration"));
  const std::string yTitle = replaceTilde(CL.Get("YTitle", "#chi^{2}"));
  const std::string bayesTitle = replaceTilde(CL.Get("BayesTitle", "Bayes~unfolding:~#chi^{2}~scan~by~validation~quarter"));
  const std::string svdTitle = replaceTilde(CL.Get("SVDTitle", "SVD~unfolding:~#chi^{2}~scan~by~validation~quarter"));
  const std::string fftTitle = replaceTilde(CL.Get("FFTTitle", "FFT~deconvolution:~#chi^{2}~scan~by~validation~quarter"));
  const std::string yTitleMeasuredOnly = replaceTilde(CL.Get("YTitleMeasuredOnly", "#chi^{2}/NDF~(measured~errors~only)"));
  const std::string bayesTitleMeasuredOnly = replaceTilde(CL.Get("BayesTitleMeasuredOnly", "Bayes~unfolding:~#chi^{2}/NDF~(measured~errors~only)~vs~iteration~for~each~validation~quarter"));
  const std::string svdTitleMeasuredOnly = replaceTilde(CL.Get("SVDTitleMeasuredOnly", "SVD~unfolding:~#chi^{2}/NDF~(measured~errors~only)~vs~iteration~for~each~validation~quarter"));
  const std::string fftTitleMeasuredOnly = replaceTilde(CL.Get("FFTTitleMeasuredOnly", "FFT~deconvolution:~#chi^{2}/NDF~(measured~errors~only)~vs~iteration~for~each~validation~quarter"));

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::vector<TFile *> ownedFiles;
  std::vector<QuarterCurve> bayesCurves;
  std::vector<QuarterCurve> svdCurves;
  std::vector<QuarterCurve> fftCurves;

  auto cleanup = [&]()
  {
    for(TFile *file : ownedFiles)
      RootIOUtils::CloseAndDeleteFile(file);
    ownedFiles.clear();
  };

  if(!loadQuarterCurves(bayesFiles, chi2HistName, "", "Bayes method", requireAllQuarters, ownedFiles, bayesCurves))
  {
    cleanup();
    return -1;
  }
  if(!loadQuarterCurves(svdFiles, chi2HistName, "", "SVD method", requireAllQuarters, ownedFiles, svdCurves))
  {
    cleanup();
    return -1;
  }
  if(!loadQuarterCurves(fftFiles, chi2HistName, "", "FFT method", requireAllQuarters, ownedFiles, fftCurves))
  {
    cleanup();
    return -1;
  }

  const std::string bayesOutput = outputPrefix + "_bayes.pdf";
  const std::string svdOutput = outputPrefix + "_svd.pdf";
  const std::string fftOutput = outputPrefix + "_fft.pdf";

  drawQuarterCurves(bayesCurves, "cBayesQuarterChi2", bayesTitle, xTitle, yTitle, bayesOutput, logY, yMinUser, yMaxUser);
  drawQuarterCurves(svdCurves, "cSVDQuarterChi2", svdTitle, xTitle, yTitle, svdOutput, logY, yMinUser, yMaxUser);
  drawQuarterCurves(fftCurves, "cFFTQuarterChi2", fftTitle, xTitle, yTitle, fftOutput, logY, yMinUser, yMaxUser);

  std::string bayesMeasuredOnlyOutput;
  std::string svdMeasuredOnlyOutput;
  std::string fftMeasuredOnlyOutput;
  if(makeMeasuredOnlyPlots)
  {
    std::vector<QuarterCurve> bayesCurvesMeasuredOnly;
    std::vector<QuarterCurve> svdCurvesMeasuredOnly;
    std::vector<QuarterCurve> fftCurvesMeasuredOnly;

    if(!loadQuarterCurves(bayesFiles, bayesMeasuredOnlyHistName, chi2HistName, "Bayes method (measured-only chi2)", requireAllQuarters, ownedFiles, bayesCurvesMeasuredOnly))
    {
      cleanup();
      return -1;
    }
    if(!loadQuarterCurves(svdFiles, svdMeasuredOnlyHistName, chi2HistName, "SVD method (measured-only chi2)", requireAllQuarters, ownedFiles, svdCurvesMeasuredOnly))
    {
      cleanup();
      return -1;
    }
    if(!loadQuarterCurves(fftFiles, fftMeasuredOnlyHistName, chi2HistName, "FFT method (measured-only chi2)", requireAllQuarters, ownedFiles, fftCurvesMeasuredOnly))
    {
      cleanup();
      return -1;
    }

    bayesMeasuredOnlyOutput = outputPrefix + "_bayes_measured_only.pdf";
    svdMeasuredOnlyOutput = outputPrefix + "_svd_measured_only.pdf";
    fftMeasuredOnlyOutput = outputPrefix + "_fft_measured_only.pdf";

    drawQuarterCurves(bayesCurvesMeasuredOnly, "cBayesQuarterChi2MeasuredOnly", bayesTitleMeasuredOnly, xTitle, yTitleMeasuredOnly, bayesMeasuredOnlyOutput, logY, yMinUser, yMaxUser);
    drawQuarterCurves(svdCurvesMeasuredOnly, "cSVDQuarterChi2MeasuredOnly", svdTitleMeasuredOnly, xTitle, yTitleMeasuredOnly, svdMeasuredOnlyOutput, logY, yMinUser, yMaxUser);
    drawQuarterCurves(fftCurvesMeasuredOnly, "cFFTQuarterChi2MeasuredOnly", fftTitleMeasuredOnly, xTitle, yTitleMeasuredOnly, fftMeasuredOnlyOutput, logY, yMinUser, yMaxUser);
  }

  cleanup();

  std::cout << "Wrote quarter-comparison chi2 plots:" << std::endl;
  std::cout << "  " << bayesOutput << std::endl;
  std::cout << "  " << svdOutput << std::endl;
  std::cout << "  " << fftOutput << std::endl;
  if(makeMeasuredOnlyPlots)
  {
    std::cout << "Wrote quarter-comparison chi2 plots (measured-only denominator):" << std::endl;
    std::cout << "  " << bayesMeasuredOnlyOutput << std::endl;
    std::cout << "  " << svdMeasuredOnlyOutput << std::endl;
    std::cout << "  " << fftMeasuredOnlyOutput << std::endl;
  }

  return 0;
}
