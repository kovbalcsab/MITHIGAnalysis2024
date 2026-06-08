#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CommandLine.h"

using namespace std;

namespace {
const int kOverlayColors[10] = {TColor::GetColor("#5790fc"), TColor::GetColor("#f89c20"), TColor::GetColor("#e42536"),
                                TColor::GetColor("#964a8b"), TColor::GetColor("#9c9ca1"), TColor::GetColor("#16a085"),
                                TColor::GetColor("#7a21dd"), TColor::GetColor("#80B5FF"), TColor::GetColor("#F05660"),
                                TColor::GetColor("#68686E")};

const int kCMSPaletteContours = 256;
int kCMSBluePalette[kCMSPaletteContours];

vector<string> SplitCommaSeparated(const vector<string> &input) {
  vector<string> result;
  for (const string &entry : input) {
    size_t start = 0;
    while (start <= entry.size()) {
      size_t end = entry.find(',', start);
      string piece = (end == string::npos) ? entry.substr(start) : entry.substr(start, end - start);
      if (piece.empty() == false)
        result.push_back(piece);
      if (end == string::npos)
        break;
      start = end + 1;
    }
  }
  return result;
}

string FormatLabel(const string &input) {
  string result = input;
  for (char &c : result)
    if (c == '_')
      c = ' ';
  return result;
}

vector<string> FormatLabels(const vector<string> &input) {
  vector<string> result;
  result.reserve(input.size());
  for (const string &entry : input)
    result.push_back(FormatLabel(entry));
  return result;
}

void BuildCMSPalettes() {
  static bool built = false;
  if (built == true)
    return;

  const int nrgbs = 3;
  Double_t stops[nrgbs] = {0.00, 0.45, 1.00};
  Double_t red[nrgbs] = {7. / 255., 51. / 255., 168. / 255.};
  Double_t green[nrgbs] = {59. / 255., 146. / 255., 233. / 255.};
  Double_t blue[nrgbs] = {136. / 255., 230. / 255., 255. / 255.};

  int firstColor = TColor::CreateGradientColorTable(nrgbs, stops, red, green, blue, kCMSPaletteContours);
  for (int i = 0; i < kCMSPaletteContours; ++i)
    kCMSBluePalette[i] = firstColor + i;

  built = true;
}
} // namespace

////////////////////////////////////////////////////////////
// CMS style (minimal clean version)
////////////////////////////////////////////////////////////
void SetCMSStyle() {
  BuildCMSPalettes();

  gStyle->SetOptStat(0);

  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetFrameLineWidth(1);

  gStyle->SetEndErrorSize(0);
  gStyle->SetNumberContours(kCMSPaletteContours);
  gStyle->SetPalette(kCMSPaletteContours, kCMSBluePalette);
}

////////////////////////////////////////////////////////////
// Histogram styling
////////////////////////////////////////////////////////////
void StyleHistogram(TH1D *hist, int color) {
  if (!hist)
    return;

  hist->SetLineColor(color);
  hist->SetLineWidth(1);

  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);

  hist->GetXaxis()->SetTitle("HF E_{max} (GeV)");
  hist->GetYaxis()->SetTitle("Entries");

  hist->SetTitle("");

  hist->GetXaxis()->SetRangeUser(0, 30);

  hist->GetYaxis()->SetRangeUser(1e-5, 1.);
}

////////////////////////////////////////////////////////////
// Draw histogram with smart axis visibility
////////////////////////////////////////////////////////////
void DrawHistogram(TH1D *hist, bool showX, bool showY, int color, const char *drawOption = "hist") {
  if (hist == nullptr)
    return;

  StyleHistogram(hist, color);

  if (!showX) {
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetXaxis()->SetTitleSize(0);
  }

  if (!showY) {
    hist->GetYaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0);
  }

  hist->Draw(drawOption);
}

////////////////////////////////////////////////////////////
// Draw style for eta-phi summary maps
////////////////////////////////////////////////////////////
void StyleEtaPhiMap(TH2D *hist) {
  if (hist == nullptr)
    return;

  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetZaxis()->SetTitleOffset(1.6);
}

void SetEtaPhiMapPadMargins(TVirtualPad *pad) {
  if (pad == nullptr)
    return;

  pad->SetLeftMargin(0.16);
  pad->SetRightMargin(0.20);
  pad->SetBottomMargin(0.12);
  pad->SetTopMargin(0.12);
}

void DrawCMSInternal(TPad *pad, bool insideFrame = true, double textSize = 0.06, double xShift = 0.0,
                     double yShift = 0.0) {
  if (pad == nullptr)
    return;

  pad->cd();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  latex.SetTextFont(42);
  latex.SetTextSize(textSize);

  const double left = pad->GetLeftMargin();
  const double top = pad->GetTopMargin();
  const double xCMS = (insideFrame == true ? (left + 0.02) : (left)) + xShift;
  const double yCMS = (insideFrame == true ? (1.0 - top - 0.04) : (1.0 - top + 0.007)) + yShift;
  const double xInternal = xCMS + 0.115;

  if (insideFrame == true) {
    latex.SetTextAlign(13);
    latex.DrawLatex(xCMS, yCMS, "#font[61]{#scale[1.25]{CMS}}");
    latex.DrawLatex(xInternal, yCMS, "#font[52]{Internal}");
  } else {
    latex.SetTextAlign(11);
    latex.DrawLatex(xCMS, yCMS, "#font[61]{#scale[1.25]{CMS}}");
    latex.DrawLatex(xInternal, yCMS, "#font[52]{Internal}");
  }
}

////////////////////////////////////////////////////////////

double gamma_pdf(double x, double k, double theta) {
  if (x < 0)
    return 0;
  return pow(x, k - 1) * exp(-x / theta) / (pow(theta, k) * TMath::Gamma(k));
}

Double_t exponential_pdf(Double_t *x, Double_t *par) {
  if (x[0] < 0)
    return 0;
  return par[1] / par[0] * exp(-x[0] / par[0]);
}

Double_t fitf(Double_t *x, Double_t *par) {
  return par[3] * gamma_pdf(x[0] - par[5], par[0], par[1]) + par[4] * exp(-x[0] / par[2]) / par[2];
}

// Calculate mean of histogram, neglecting the first bin (ibin = 1)
double calculateMean(const TH1D *hist) {
  if (hist == nullptr)
    return 0.0;

  double weightedSum = 0.0;
  double denominator = 0.0;

  for (int ibin = 2; ibin <= hist->GetNbinsX(); ibin++) {
    const double content = hist->GetBinContent(ibin);
    weightedSum += content * hist->GetBinCenter(ibin);
    denominator += content;
  }
  return denominator > 0.0 ? weightedSum / denominator : 0.0;
}

double calculateStdDev(const TH1D *hist, double mean) {
  if (hist == nullptr)
    return 0.0;

  double numerator = 0.0;
  double denominator = 0.0;

  for (int ibin = 2; ibin <= hist->GetNbinsX(); ibin++) {
    const double content = hist->GetBinContent(ibin);
    numerator += content * pow(hist->GetBinCenter(ibin) - mean, 2);
    denominator += content;
  }
  return denominator > 0.0 ? sqrt(numerator / denominator) : 0.0;
}

double calculateIntegral(const TH1D *hist) {
  if (hist == nullptr)
    return 0.0;
  return hist->Integral(2, hist->GetNbinsX());
}

void NormalizeHistogram(TH1D *hist) {
  if (hist == nullptr)
    return;

  const double integral = hist->Integral();
  if (integral > 0)
    hist->Scale(1.0 / integral);
}

TH1D *CloneHistogram(TH1D *hist, const string &name) {
  if (hist == nullptr)
    return nullptr;

  TH1D *clone = (TH1D *)hist->Clone(name.c_str());
  if (clone != nullptr)
    clone->SetDirectory(nullptr);
  return clone;
}

void StyleSimple1DHistogram(TH1D *hist, const string &xTitle, int color, bool doSelfNormalize) {
  if (hist == nullptr)
    return;

  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(xTitle.c_str());
  hist->GetYaxis()->SetTitle(doSelfNormalize == true ? "Normalized entries" : "Entries");
  hist->GetXaxis()->SetTitleSize(0.060);
  hist->GetXaxis()->SetLabelSize(0.050);
  hist->GetYaxis()->SetTitleSize(0.060);
  hist->GetYaxis()->SetLabelSize(0.050);
  hist->GetYaxis()->SetTitleOffset(1.10);
  hist->SetMinimum(1e-5);
}

void StyleNoRatioFullPadHistogram(TH1D *hist, const string &xTitle, int color, bool doSelfNormalize) {
  if (hist == nullptr)
    return;

  StyleSimple1DHistogram(hist, xTitle, color, doSelfNormalize);
  hist->GetXaxis()->SetTitleSize(0.042);
  hist->GetXaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTitleSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTitleOffset(1.20);
}

void StyleRatioHistogram(TH1D *hist, const string &xTitle, int color) {
  if (hist == nullptr)
    return;

  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(xTitle.c_str());
  hist->GetYaxis()->SetTitle("Ratio");
  hist->GetYaxis()->SetTitleSize(0.14);
  hist->GetYaxis()->SetTitleOffset(0.45);
  hist->GetYaxis()->SetLabelSize(0.117);
  hist->GetYaxis()->SetLabelOffset(0.012);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetTitleSize(0.14);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetXaxis()->SetLabelSize(0.117);
}

double GetDynamicXUpperEdge(const vector<TH1D *> &hists) {
  double maximumX = 0.0;
  double axisMaximum = 0.0;

  for (TH1D *hist : hists) {
    if (hist == nullptr)
      continue;

    axisMaximum = max(axisMaximum, hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()));

    for (int i = hist->GetNbinsX(); i >= 1; --i) {
      if (hist->GetBinContent(i) == 0.0)
        continue;

      maximumX = max(maximumX, hist->GetXaxis()->GetBinUpEdge(i));
      break;
    }
  }

  if (maximumX <= 0.0) {
    for (TH1D *hist : hists) {
      if (hist == nullptr)
        continue;
      maximumX = hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX());
      break;
    }
  }

  return min(maximumX * 1.1, axisMaximum);
}

void DrawMultiHistogramOverlay(const vector<TH1D *> &inputHistograms, const string &outputFileName,
                               const string &xTitle, const vector<string> &labels = {}, double fixedXUpper = -1.0,
                               bool doRatio = false, bool doSelfNormalize = true, bool doBinomialRatio = false,
                               double ratioDistanceAroundUnity = -1.0) {
  vector<TH1D *> histograms;
  vector<TH1D *> ratioSources;
  for (size_t i = 0; i < inputHistograms.size(); ++i) {
    if (inputHistograms[i] == nullptr)
      continue;

    TH1D *displayHistogram =
        CloneHistogram(inputHistograms[i], Form("%s_display_%zu", inputHistograms[i]->GetName(), i));
    TH1D *ratioHistogram =
        CloneHistogram(inputHistograms[i], Form("%s_ratioSource_%zu", inputHistograms[i]->GetName(), i));
    if (displayHistogram == nullptr || ratioHistogram == nullptr)
      continue;

    histograms.push_back(displayHistogram);
    ratioSources.push_back(ratioHistogram);
  }

  if (histograms.empty() == true)
    return;

  if (doSelfNormalize == true) {
    for (TH1D *hist : histograms)
      NormalizeHistogram(hist);
  }

  const double xUpperEdge = (fixedXUpper > 0.0) ? fixedXUpper : GetDynamicXUpperEdge(histograms);
  double maximumY = 0.0;
  for (TH1D *hist : histograms)
    maximumY = max(maximumY, hist->GetMaximum());

  for (size_t i = 0; i < histograms.size(); ++i) {
    StyleSimple1DHistogram(histograms[i], xTitle, kOverlayColors[min<size_t>(i, 9)], doSelfNormalize);
    histograms[i]->SetMaximum(maximumY * 1.15);
    histograms[i]->GetXaxis()->SetRangeUser(0.0, xUpperEdge);
  }

  TCanvas canvas("c1DOverlay", "", 900, 900);
  const bool drawRatio = (doRatio == true && histograms.size() > 1);

  if (drawRatio == true) {
    canvas.cd();
    TPad topPad("TopPad", "", 0.0, 0.30, 1.0, 1.0);
    TPad bottomPad("BottomPad", "", 0.0, 0.0, 1.0, 0.30);
    topPad.SetLeftMargin(0.14);
    topPad.SetRightMargin(0.04);
    topPad.SetBottomMargin(0.02);
    topPad.SetTopMargin(0.12);
    topPad.SetLogy();
    bottomPad.SetLeftMargin(0.14);
    bottomPad.SetRightMargin(0.04);
    bottomPad.SetBottomMargin(0.35);
    bottomPad.SetTopMargin(0.03);
    topPad.Draw();
    bottomPad.Draw();

    topPad.cd();
    for (TH1D *hist : histograms) {
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetTitleSize(0);
    }
    histograms[0]->Draw("hist");
    for (size_t i = 1; i < histograms.size(); ++i)
      histograms[i]->Draw("hist same");
    histograms[0]->GetYaxis()->ChangeLabel(1, -1, 0.0, -1, -1, -1, "");

    if (labels.empty() == false) {
      TLegend *legend = new TLegend(0.52, 0.64, 0.82, 0.85);
      legend->SetBorderSize(0);
      legend->SetLineColor(0);
      legend->SetFillStyle(0);
      legend->SetFillColorAlpha(0, 0.0);
      legend->SetTextSize(0.045);
      for (size_t i = 0; i < histograms.size() && i < labels.size(); ++i)
        legend->AddEntry(histograms[i], labels[i].c_str(), "l");
      legend->Draw();
    }

    bottomPad.cd();
    bool firstRatio = true;
    const double ratioMinimum = (ratioDistanceAroundUnity > 0.0) ? (1.0 - ratioDistanceAroundUnity) : 0.0;
    const double ratioMaximum = (ratioDistanceAroundUnity > 0.0) ? (1.0 + ratioDistanceAroundUnity) : 2.0;
    for (size_t i = 1; i < histograms.size(); ++i) {
      TH1D *ratioBaseNumerator = (doBinomialRatio == true) ? ratioSources[i] : histograms[i];
      TH1D *ratioBaseDenominator = (doBinomialRatio == true) ? ratioSources[0] : histograms[0];
      TH1D *ratio = (TH1D *)ratioBaseNumerator->Clone(Form("%s_ratio_%zu", ratioBaseNumerator->GetName(), i));
      if (ratio == nullptr)
        continue;
      ratio->SetDirectory(nullptr);
      ratio->Sumw2();
      ratio->Divide(ratioBaseNumerator, ratioBaseDenominator, 1.0, 1.0, doBinomialRatio == true ? "B" : "");
      StyleRatioHistogram(ratio, xTitle, kOverlayColors[min<size_t>(i, 9)]);
      ratio->GetXaxis()->SetRangeUser(0.0, xUpperEdge);
      ratio->SetMinimum(ratioMinimum);
      ratio->SetMaximum(ratioMaximum);
      ratio->GetYaxis()->SetRangeUser(ratioMinimum, ratioMaximum);
      ratio->Draw(firstRatio ? "E1" : "E1 same");
      firstRatio = false;
    }

    TLine unityLine(0.0, 1.0, xUpperEdge, 1.0);
    unityLine.SetLineStyle(2);
    unityLine.SetLineWidth(2);
    unityLine.Draw();

    TLatex ratioLabel;
    ratioLabel.SetNDC();
    ratioLabel.SetTextFont(42);
    ratioLabel.SetTextSize(0.09);
    ratioLabel.DrawLatex(0.16, 0.82,
                         Form("Denominator: %s", (labels.empty() == false ? labels[0].c_str() : "input 1")));
    topPad.Modified();
    bottomPad.Modified();
    canvas.cd();
    canvas.Modified();
    canvas.Update();
    DrawCMSInternal(&topPad, false);
    canvas.Modified();
    canvas.Update();
    canvas.SaveAs(outputFileName.c_str());
  } else {
    canvas.cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.12);

    for (size_t i = 0; i < histograms.size(); ++i)
      StyleNoRatioFullPadHistogram(histograms[i], xTitle, kOverlayColors[min<size_t>(i, 9)], doSelfNormalize);

    histograms[0]->Draw("hist");
    for (size_t i = 1; i < histograms.size(); ++i)
      histograms[i]->Draw("hist same");

    if (labels.empty() == false) {
      TLegend *legend = new TLegend(0.52, 0.60, 0.82, 0.78);
      legend->SetBorderSize(0);
      legend->SetLineColor(0);
      legend->SetFillStyle(0);
      legend->SetFillColorAlpha(0, 0.0);
      legend->SetTextSize(0.045);
      for (size_t i = 0; i < histograms.size() && i < labels.size(); ++i)
        legend->AddEntry(histograms[i], labels[i].c_str(), "l");
      legend->Draw();
    }
    DrawCMSInternal((TPad *)gPad, false, 0.042, 0.0, 0.005);
  }

  if (drawRatio == false) {
    canvas.Update();
    canvas.SaveAs(outputFileName.c_str());
  }
}

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  vector<string> InputFileNames = SplitCommaSeparated(CL.GetStringVector("Input", vector<string>()));
  if (InputFileNames.empty() == true) {
    const string InputFileName = CL.Get("Input", CL.Get("Input1", ""));
    const string InputFileName2 = CL.Get("Input2", "");
    const string InputFileName3 = CL.Get("Input3", "");
    if (InputFileName.empty() == false)
      InputFileNames.push_back(InputFileName);
    if (InputFileName2.empty() == false)
      InputFileNames.push_back(InputFileName2);
    if (InputFileName3.empty() == false)
      InputFileNames.push_back(InputFileName3);
  }
  if (InputFileNames.size() > 10)
    InputFileNames.resize(10);
  vector<string> Labels = SplitCommaSeparated(CL.GetStringVector("Labels", vector<string>()));
  Labels = FormatLabels(Labels);
  string OutputFileName = CL.Get("Output");
  bool doCumulative = CL.GetBool("DoCumulative", false);
  bool doRatio = CL.GetBool("DoRatio", false);
  bool doSelfNormalize = CL.GetBool("DoSelfNormalize", true);
  bool doBinomialRatio = CL.GetBool("DoBinomialRatio", false);
  double ratioDistanceAroundUnity = CL.GetDouble("RatioDistanceAroundUnity", -1.0);

  if (InputFileNames.empty() == true) {
    cerr << "No input files provided" << endl;
    return 1;
  }

  SetCMSStyle();
  gSystem->mkdir(OutputFileName.c_str(), true);

  vector<unique_ptr<TFile>> InputFiles;
  InputFiles.reserve(InputFileNames.size());
  for (const string &name : InputFileNames) {
    InputFiles.emplace_back(TFile::Open(name.c_str()));
    if (InputFiles.back() == nullptr || InputFiles.back()->IsZombie() == true) {
      cerr << "Failed to open input file " << name << endl;
      return 1;
    }
  }
  TFile *InputFile = InputFiles[0].get();

  auto LoadHistogramVector = [&](const string &name) {
    vector<TH1D *> result;
    result.reserve(InputFiles.size());
    for (const auto &file : InputFiles)
      result.push_back((file != nullptr) ? (TH1D *)file->Get(name.c_str()) : nullptr);
    return result;
  };

  vector<TH1D *> hPlusLeadingAll = LoadHistogramVector("hHFEMaxPlusLeading");
  vector<TH1D *> hMinusLeadingAll = LoadHistogramVector("hHFEMaxMinusLeading");
  TH1D *hPlusLeading = hPlusLeadingAll[0];
  TH1D *hMinusLeading = hMinusLeadingAll[0];
  TH1D *hHFnPF = (TH1D *)InputFile->Get("hHFnPF");
  TH1D *hHFnPF6p = (TH1D *)InputFile->Get("hHFnPF6p");
  TH1D *hHFnPF7p = (TH1D *)InputFile->Get("hHFnPF7p");
  TH1D *hHFnPF6m = (TH1D *)InputFile->Get("hHFnPF6m");
  TH1D *hHFnPF7m = (TH1D *)InputFile->Get("hHFnPF7m");
  TH1D *hAllEnergyp = (TH1D *)InputFile->Get("hAllEnergyp");
  TH1D *hAllEnergym = (TH1D *)InputFile->Get("hAllEnergym");
  TH1D *hAllEnergy6p = (TH1D *)InputFile->Get("hAllEnergy6p");
  TH1D *hAllEnergy7p = (TH1D *)InputFile->Get("hAllEnergy7p");
  TH1D *hAllEnergy6m = (TH1D *)InputFile->Get("hAllEnergy6m");
  TH1D *hAllEnergy7m = (TH1D *)InputFile->Get("hAllEnergy7m");

  ////////////////////////////////////////////////////////////
  // binning
  ////////////////////////////////////////////////////////////

  vector<float> etaBorders = {3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2};

  vector<float> phiBorders = {-M_PI, -2 * M_PI / 3, -M_PI / 3, 0, M_PI / 3, 2 * M_PI / 3, M_PI};

  const int nEta = etaBorders.size() - 1;
  const int nPhi = phiBorders.size() - 1;
  const TString firstPlusName =
      Form("hHFEMaxPlus_eta%.1f_%.1f_phi%.1f_%.1f", etaBorders[0], etaBorders[1], phiBorders[0], phiBorders[1]);
  const bool hasEtaPhiHistograms = (InputFile->Get(firstPlusName) != nullptr);

  if (hPlusLeading != nullptr && hMinusLeading != nullptr) {
    if (doCumulative == true) {
      for (size_t i = 0; i < hPlusLeadingAll.size(); ++i) {
        if (hPlusLeadingAll[i] != nullptr)
          hPlusLeadingAll[i] = (TH1D *)hPlusLeadingAll[i]->GetCumulative(kFALSE);
        if (hMinusLeadingAll[i] != nullptr)
          hMinusLeadingAll[i] = (TH1D *)hMinusLeadingAll[i]->GetCumulative(kFALSE);
      }
    }

    DrawMultiHistogramOverlay(hPlusLeadingAll, Form("%s/HFEMaxPlusLeadingOverlap.pdf", OutputFileName.c_str()),
                              "HF E^{+}_{max} [GeV]", Labels, 30.0, doRatio, doSelfNormalize, doBinomialRatio,
                              ratioDistanceAroundUnity);
    DrawMultiHistogramOverlay(hMinusLeadingAll, Form("%s/HFEMaxMinusLeadingOverlap.pdf", OutputFileName.c_str()),
                              "HF E^{-}_{max} [GeV]", Labels, 30.0, doRatio, doSelfNormalize, doBinomialRatio,
                              ratioDistanceAroundUnity);
    if (hasEtaPhiHistograms == false)
      return 0;
  }

  if (hHFnPF != nullptr) {
    DrawMultiHistogramOverlay(LoadHistogramVector("hHFnPF"), Form("%s/HFnPFOverlap.pdf", OutputFileName.c_str()),
                              "HF nPF", Labels, -1.0, doRatio, doSelfNormalize, doBinomialRatio,
                              ratioDistanceAroundUnity);
  }
  vector<pair<vector<TH1D *>, pair<string, string>>> extra1D = {
      {LoadHistogramVector("hHFnPF6p"), {"HFnPF6pOverlap.pdf", "HF nPF id6+"}},
      {LoadHistogramVector("hHFnPF7p"), {"HFnPF7pOverlap.pdf", "HF nPF id7+"}},
      {LoadHistogramVector("hHFnPF6m"), {"HFnPF6mOverlap.pdf", "HF nPF id6-"}},
      {LoadHistogramVector("hHFnPF7m"), {"HFnPF7mOverlap.pdf", "HF nPF id7-"}},
      {LoadHistogramVector("hAllEnergyp"), {"AllEnergypOverlap.pdf", "HF PF candidate energy + (GeV)"}},
      {LoadHistogramVector("hAllEnergym"), {"AllEnergymOverlap.pdf", "HF PF candidate energy - (GeV)"}},
      {LoadHistogramVector("hAllEnergy6p"), {"AllEnergy6pOverlap.pdf", "HF PF candidate energy id6+ (GeV)"}},
      {LoadHistogramVector("hAllEnergy7p"), {"AllEnergy7pOverlap.pdf", "HF PF candidate energy id7+ (GeV)"}},
      {LoadHistogramVector("hAllEnergy6m"), {"AllEnergy6mOverlap.pdf", "HF PF candidate energy id6- (GeV)"}},
      {LoadHistogramVector("hAllEnergy7m"), {"AllEnergy7mOverlap.pdf", "HF PF candidate energy id7- (GeV)"}}};
  for (const auto &entry : extra1D) {
    DrawMultiHistogramOverlay(entry.first, Form("%s/%s", OutputFileName.c_str(), entry.second.first.c_str()),
                              entry.second.second, Labels, -1.0, doRatio, doSelfNormalize, doBinomialRatio,
                              ratioDistanceAroundUnity);
  }
  vector<pair<pair<TH1D *, TH1D *>, pair<string, string>>> sideOverlaps = {
      {{hHFnPF6p, hHFnPF6m}, {"HFnPF6SideOverlap.pdf", "HF nPF id6"}},
      {{hHFnPF7p, hHFnPF7m}, {"HFnPF7SideOverlap.pdf", "HF nPF id7"}},
      {{hAllEnergyp, hAllEnergym}, {"AllEnergySideOverlap.pdf", "HF PF candidate energy (GeV)"}},
      {{hAllEnergy6p, hAllEnergy6m}, {"AllEnergy6SideOverlap.pdf", "HF PF candidate energy id6 (GeV)"}},
      {{hAllEnergy7p, hAllEnergy7m}, {"AllEnergy7SideOverlap.pdf", "HF PF candidate energy id7 (GeV)"}}};
  for (const auto &entry : sideOverlaps) {
    if (entry.first.first != nullptr && entry.first.second != nullptr)
      DrawMultiHistogramOverlay({entry.first.first, entry.first.second},
                                Form("%s/%s", OutputFileName.c_str(), entry.second.first.c_str()), entry.second.second,
                                {}, -1.0, doRatio, doSelfNormalize, doBinomialRatio, ratioDistanceAroundUnity);
  }

  ////////////////////////////////////////////////////////////
  // load histograms
  ////////////////////////////////////////////////////////////

  const size_t NInput = InputFiles.size();
  vector<vector<vector<TH1D *>>> hPlusAll(NInput, vector<vector<TH1D *>>(nEta, vector<TH1D *>(nPhi, nullptr)));
  vector<vector<vector<TH1D *>>> hMinusAll(NInput, vector<vector<TH1D *>>(nEta, vector<TH1D *>(nPhi, nullptr)));
  vector<vector<vector<TH1D *>>> hPlusCumAll(NInput, vector<vector<TH1D *>>(nEta, vector<TH1D *>(nPhi, nullptr)));
  vector<vector<vector<TH1D *>>> hMinusCumAll(NInput, vector<vector<TH1D *>>(nEta, vector<TH1D *>(nPhi, nullptr)));

  for (int iEta = 0; iEta < nEta; iEta++)
    for (int iPhi = 0; iPhi < nPhi; iPhi++) {
      const TString plusName = Form("hHFEMaxPlus_eta%.1f_%.1f_phi%.1f_%.1f", etaBorders[iEta], etaBorders[iEta + 1],
                                    phiBorders[iPhi], phiBorders[iPhi + 1]);
      const TString minusName = Form("hHFEMaxMinus_eta%.1f_%.1f_phi%.1f_%.1f", -etaBorders[iEta + 1], -etaBorders[iEta],
                                     phiBorders[iPhi], phiBorders[iPhi + 1]);

      for (size_t iInput = 0; iInput < NInput; ++iInput) {
        hPlusAll[iInput][iEta][iPhi] = (TH1D *)InputFiles[iInput]->Get(plusName);
        if (hPlusAll[iInput][iEta][iPhi] != nullptr && hPlusAll[iInput][iEta][iPhi]->Integral() > 0)
          hPlusAll[iInput][iEta][iPhi]->Scale(1.0 / hPlusAll[iInput][iEta][iPhi]->Integral());
        if (doCumulative && hPlusAll[iInput][iEta][iPhi] != nullptr)
          hPlusCumAll[iInput][iEta][iPhi] = (TH1D *)hPlusAll[iInput][iEta][iPhi]->GetCumulative(kFALSE);

        hMinusAll[iInput][iEta][iPhi] = (TH1D *)InputFiles[iInput]->Get(minusName);
        if (hMinusAll[iInput][iEta][iPhi] != nullptr && hMinusAll[iInput][iEta][iPhi]->Integral() > 0)
          hMinusAll[iInput][iEta][iPhi]->Scale(1.0 / hMinusAll[iInput][iEta][iPhi]->Integral());
        if (doCumulative && hMinusAll[iInput][iEta][iPhi] != nullptr)
          hMinusCumAll[iInput][iEta][iPhi] = (TH1D *)hMinusAll[iInput][iEta][iPhi]->GetCumulative(kFALSE);
      }
    }

  ////////////////////////////////////////////////////////////
  // pad geometry (fix squishing)
  ////////////////////////////////////////////////////////////

  int padW = 300;
  int padH = 300;

  double leftMargin = 0.16;
  double rightMargin = 0.02;
  double bottomMargin = 0.16;
  double topMargin = 0.10;

  int canvasW = padW * nEta / (1 - leftMargin - rightMargin);
  int canvasH = padH * nPhi / (1 - bottomMargin - topMargin);

  TCanvas *cPlus = new TCanvas("cPlus", "", canvasW, canvasH);
  TCanvas *cMinus = new TCanvas("cMinus", "", canvasW, canvasH);
  cPlus->SetTopMargin(0.10);
  cMinus->SetTopMargin(0.10);

  cPlus->Divide(nEta, nPhi, 0, 0);
  cMinus->Divide(nEta, nPhi, 0, 0);

  TLatex *label = new TLatex();
  label->SetNDC();
  label->SetTextSize(0.05);
  label->SetTextFont(42);

  std::vector<std::vector<TF1 *>> fitResultsPlus(nEta, std::vector<TF1 *>(nPhi, nullptr));
  std::vector<std::vector<TF1 *>> fitResultsMinus(nEta, std::vector<TF1 *>(nPhi, nullptr));
  TH2D *hEtaPhiMeanPlus =
      new TH2D("hEtaPhiMeanPlus", ";#eta;#phi;Mean HF E_{max} (GeV)", nEta, etaBorders.data(), nPhi, phiBorders.data());
  TH2D *hEtaPhiStdPlus = new TH2D("hEtaPhiStdPlus", ";#eta;#phi;Std Dev of HF E_{max} (GeV)", nEta, etaBorders.data(),
                                  nPhi, phiBorders.data());
  TH2D *hEtaPhiNumPlus = new TH2D("hEtaPhiNumPlus", ";#eta;#phi;Integral excluding first bin", nEta, etaBorders.data(),
                                  nPhi, phiBorders.data());
  TH2D *hEtaPhiMeanMinus = new TH2D("hEtaPhiMeanMinus", ";#eta;#phi;Mean HF E_{max} (GeV)", nEta, etaBorders.data(),
                                    nPhi, phiBorders.data());
  TH2D *hEtaPhiStdMinus = new TH2D("hEtaPhiStdMinus", ";#eta;#phi;Std Dev of HF E_{max} (GeV)", nEta, etaBorders.data(),
                                   nPhi, phiBorders.data());
  TH2D *hEtaPhiNumMinus = new TH2D("hEtaPhiNumMinus", ";#eta;#phi;Integral excluding first bin", nEta,
                                   etaBorders.data(), nPhi, phiBorders.data());

  ////////////////////////////////////////////////////////////
  // drawing
  ////////////////////////////////////////////////////////////

  for (int iEta = 0; iEta < nEta; iEta++)
    for (int iPhi = 0; iPhi < nPhi; iPhi++) {

      int padID = iPhi * nEta + iEta + 1;

      bool showY = (iEta == 0);
      bool showX = (iPhi == nPhi - 1);

      ///////////////////////
      // PLUS
      ///////////////////////

      cPlus->cd(padID);

      gPad->SetLogy();

      gPad->SetLeftMargin(showY ? leftMargin : 0.02);
      gPad->SetRightMargin(rightMargin);
      gPad->SetBottomMargin(showX ? bottomMargin : 0.02);
      gPad->SetTopMargin(topMargin);

      gPad->SetFrameLineWidth(1);

      TH1D *histToDrawPlus = doCumulative ? hPlusCumAll[0][iEta][iPhi] : hPlusAll[0][iEta][iPhi];
      DrawHistogram(histToDrawPlus, showX, showY, kBlue, "hist");
      for (size_t iInput = 1; iInput < NInput; ++iInput) {
        TH1D *histToDrawPlusN = doCumulative ? hPlusCumAll[iInput][iEta][iPhi] : hPlusAll[iInput][iEta][iPhi];
        DrawHistogram(histToDrawPlusN, showX, showY, kOverlayColors[min<size_t>(iInput, 9)], "hist same");
      }
      label->DrawLatex(0.15, 0.85,
                       Form("#eta: [%.1f,%.1f], #phi: [%.1f,%.1f]", etaBorders[iEta], etaBorders[iEta + 1],
                            phiBorders[iPhi], phiBorders[iPhi + 1]));
      const double plusMean = calculateMean(hPlusAll[0][iEta][iPhi]);
      hEtaPhiMeanPlus->SetBinContent(iEta + 1, iPhi + 1, plusMean);
      hEtaPhiStdPlus->SetBinContent(iEta + 1, iPhi + 1, calculateStdDev(hPlusAll[0][iEta][iPhi], plusMean));
      hEtaPhiNumPlus->SetBinContent(iEta + 1, iPhi + 1, calculateIntegral(hPlusAll[0][iEta][iPhi]));

      ///////////////////////
      // MINUS
      ///////////////////////

      cMinus->cd(padID);

      gPad->SetLogy();

      gPad->SetLeftMargin(showY ? leftMargin : 0.02);
      gPad->SetRightMargin(rightMargin);
      gPad->SetBottomMargin(showX ? bottomMargin : 0.02);
      gPad->SetTopMargin(topMargin);

      gPad->SetFrameLineWidth(1);

      TH1D *histToDrawMinus = doCumulative ? hMinusCumAll[0][iEta][iPhi] : hMinusAll[0][iEta][iPhi];
      DrawHistogram(histToDrawMinus, showX, showY, kBlue, "hist");
      for (size_t iInput = 1; iInput < NInput; ++iInput) {
        TH1D *histToDrawMinusN = doCumulative ? hMinusCumAll[iInput][iEta][iPhi] : hMinusAll[iInput][iEta][iPhi];
        DrawHistogram(histToDrawMinusN, showX, showY, kOverlayColors[min<size_t>(iInput, 9)], "hist same");
      }
      label->DrawLatex(0.15, 0.85,
                       Form("#eta: [%.1f,%.1f], #phi: [%.1f,%.1f]", -etaBorders[iEta + 1], -etaBorders[iEta],
                            phiBorders[iPhi], phiBorders[iPhi + 1]));
      const double minusMean = calculateMean(hMinusAll[0][iEta][iPhi]);
      hEtaPhiMeanMinus->SetBinContent(iEta + 1, iPhi + 1, minusMean);
      hEtaPhiStdMinus->SetBinContent(iEta + 1, iPhi + 1, calculateStdDev(hMinusAll[0][iEta][iPhi], minusMean));
      hEtaPhiNumMinus->SetBinContent(iEta + 1, iPhi + 1, calculateIntegral(hMinusAll[0][iEta][iPhi]));
    }
  TCanvas *cMeanPlus = new TCanvas("cMeanPlus", "", 1000, 800);
  cMeanPlus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiMeanPlus);
  hEtaPhiMeanPlus->Draw("colz");
  DrawCMSInternal((TPad *)gPad, false);
  cMeanPlus->SaveAs(Form("%s/HFEMaxPlusMeanMap.pdf", OutputFileName.c_str()));

  TCanvas *cStdPlus = new TCanvas("cStdPlus", "", 1000, 800);
  cStdPlus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiStdPlus);
  hEtaPhiStdPlus->Draw("colz");
  DrawCMSInternal((TPad *)gPad, false);
  cStdPlus->SaveAs(Form("%s/HFEMaxPlusStdMap.pdf", OutputFileName.c_str()));

  TCanvas *cNumPlus = new TCanvas("cNumPlus", "", 1000, 800);
  cNumPlus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiNumPlus);
  hEtaPhiNumPlus->Draw("colz");
  DrawCMSInternal((TPad *)gPad, false);
  cNumPlus->SaveAs(Form("%s/HFEMaxPlusNumEntriesMap.pdf", OutputFileName.c_str()));

  TCanvas *cMeanMinus = new TCanvas("cMeanMinus", "", 1000, 800);
  cMeanMinus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiMeanMinus);
  hEtaPhiMeanMinus->Draw("colz");
  DrawCMSInternal((TPad *)gPad, false);
  cMeanMinus->SaveAs(Form("%s/HFEMaxMinusMeanMap.pdf", OutputFileName.c_str()));

  TCanvas *cStdMinus = new TCanvas("cStdMinus", "", 1000, 800);
  cStdMinus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiStdMinus);
  hEtaPhiStdMinus->Draw("colz");
  DrawCMSInternal((TPad *)gPad, false);
  cStdMinus->SaveAs(Form("%s/HFEMaxMinusStdMap.pdf", OutputFileName.c_str()));

  TCanvas *cNumMinus = new TCanvas("cNumMinus", "", 1000, 800);
  cNumMinus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiNumMinus);
  hEtaPhiNumMinus->Draw("colz");
  DrawCMSInternal((TPad *)gPad, false);
  cNumMinus->SaveAs(Form("%s/HFEMaxMinusNumEntriesMap.pdf", OutputFileName.c_str()));

  ////////////////////////////////////////////////////////////
  // save
  ////////////////////////////////////////////////////////////

  cPlus->cd();
  DrawCMSInternal((TPad *)cPlus, false);
  cPlus->SaveAs(Form("%s/HFEMaxPlusMaps.pdf", OutputFileName.c_str()));

  cMinus->cd();
  DrawCMSInternal((TPad *)cMinus, false);
  cMinus->SaveAs(Form("%s/HFEMaxMinusMaps.pdf", OutputFileName.c_str()));

  cPlus->SaveAs(Form("%s/HFEMaxPlusMaps.png", OutputFileName.c_str()));

  cMinus->SaveAs(Form("%s/HFEMaxMinusMaps.png", OutputFileName.c_str()));

  return 0;
}
