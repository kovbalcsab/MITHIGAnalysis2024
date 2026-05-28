#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "CommandLine.h"

using namespace std;

namespace {
void NormalizeHistogram(TH1D *Histogram) {
  if (Histogram == nullptr)
    return;

  const double Integral = Histogram->Integral();
  if (Integral > 0.0)
    Histogram->Scale(1.0 / Integral);
}

void StyleHistogram(TH1D *Histogram, int Color, const string &XTitle) {
  if (Histogram == nullptr)
    return;

  Histogram->SetDirectory(nullptr);
  Histogram->SetStats(0);
  Histogram->SetLineColor(Color);
  Histogram->SetLineWidth(2);
  Histogram->SetMarkerColor(Color);
  Histogram->SetMarkerStyle(20);
  Histogram->SetMarkerSize(0.8);
  Histogram->SetTitle("");
  Histogram->GetXaxis()->SetTitle(XTitle.c_str());
  Histogram->GetYaxis()->SetTitle("Entries");
  Histogram->GetYaxis()->SetTitleOffset(1.3);
}

void StyleRatioHistogram(TH1D *Histogram, int Color, const string &XTitle) {
  if (Histogram == nullptr)
    return;

  Histogram->SetDirectory(nullptr);
  Histogram->SetStats(0);
  Histogram->SetLineColor(Color);
  Histogram->SetLineWidth(1);
  Histogram->SetMarkerColor(Color);
  Histogram->SetMarkerStyle(20);
  Histogram->SetMarkerSize(0.7);
  Histogram->SetTitle("");
  Histogram->GetXaxis()->SetTitle(XTitle.c_str());
  Histogram->GetYaxis()->SetTitle("Ratio");
  Histogram->GetYaxis()->SetTitleSize(0.10);
  Histogram->GetYaxis()->SetTitleOffset(0.45);
  Histogram->GetYaxis()->SetLabelSize(0.08);
  Histogram->GetYaxis()->SetNdivisions(505);
  Histogram->GetXaxis()->SetTitleSize(0.10);
  Histogram->GetXaxis()->SetTitleOffset(1.0);
  Histogram->GetXaxis()->SetLabelSize(0.08);
}

double GetDynamicXUpperEdge(const vector<TH1D *> &Histograms) {
  double MaximumX = 0.0;
  double AxisMaximum = 0.0;

  for (TH1D *Histogram : Histograms) {
    if (Histogram == nullptr)
      continue;

    AxisMaximum = max(AxisMaximum, Histogram->GetXaxis()->GetBinUpEdge(Histogram->GetNbinsX()));

    for (int i = Histogram->GetNbinsX(); i >= 1; --i) {
      if (Histogram->GetBinContent(i) == 0.0)
        continue;
      MaximumX = max(MaximumX, Histogram->GetXaxis()->GetBinUpEdge(i));
      break;
    }
  }

  if (MaximumX <= 0.0)
    MaximumX = AxisMaximum;

  return min(MaximumX * 1.1, AxisMaximum);
}

double GetPositiveMinimum(const vector<TH1D *> &Histograms) {
  double Minimum = -1.0;

  for (TH1D *Histogram : Histograms) {
    if (Histogram == nullptr)
      continue;

    for (int i = 1; i <= Histogram->GetNbinsX(); ++i) {
      const double Value = Histogram->GetBinContent(i);
      if (Value <= 0.0)
        continue;
      if (Minimum < 0.0 || Value < Minimum)
        Minimum = Value;
    }
  }

  return Minimum;
}

TH1D *LoadHistogram(TFile &InputFile, const string &Name, const string &CloneName) {
  TH1D *Histogram = dynamic_cast<TH1D *>(InputFile.Get(Name.c_str()));
  if (Histogram == nullptr)
    return nullptr;

  Histogram = dynamic_cast<TH1D *>(Histogram->Clone(CloneName.c_str()));
  if (Histogram != nullptr)
    Histogram->SetDirectory(nullptr);
  return Histogram;
}

TH1D *CloneHistogram(TH1D *Histogram, const string &CloneName) {
  if (Histogram == nullptr)
    return nullptr;

  TH1D *Clone = dynamic_cast<TH1D *>(Histogram->Clone(CloneName.c_str()));
  if (Clone != nullptr)
    Clone->SetDirectory(nullptr);
  return Clone;
}

TH1D *BuildRatioHistogram(TH1D *Numerator, TH1D *Denominator, const string &Name) {
  if (Numerator == nullptr || Denominator == nullptr)
    return nullptr;

  TH1D *Ratio = dynamic_cast<TH1D *>(Numerator->Clone(Name.c_str()));
  if (Ratio == nullptr)
    return nullptr;

  Ratio->SetDirectory(nullptr);
  Ratio->Sumw2();
  Ratio->Divide(Denominator);
  return Ratio;
}

pair<double, double> GetRatioRange(const vector<TH1D *> &Histograms, double XUpper) {
  vector<double> Values;

  for (TH1D *Histogram : Histograms) {
    if (Histogram == nullptr)
      continue;

    for (int i = 1; i <= Histogram->GetNbinsX(); ++i) {
      if (Histogram->GetXaxis()->GetBinCenter(i) > XUpper)
        continue;

      const double Value = Histogram->GetBinContent(i);
      if (std::isfinite(Value) == false || Value <= 0.0)
        continue;
      Values.push_back(Value);
    }
  }

  if (Values.empty() == true)
    return {0.8, 1.2};

  sort(Values.begin(), Values.end());
  const size_t LowIndex = min(Values.size() - 1, Values.size() / 20);
  const size_t HighIndex = min(Values.size() - 1, (Values.size() * 19) / 20);
  double Minimum = Values[LowIndex];
  double Maximum = Values[HighIndex];

  Minimum = min(Minimum * 0.9, 0.95);
  Maximum = max(Maximum * 1.1, 1.05);
  Maximum = min(Maximum, 3.0);

  if (Maximum - Minimum < 0.2) {
    Minimum = min(Minimum, 0.9);
    Maximum = max(Maximum, 1.1);
  }

  return {max(0.0, Minimum), Maximum};
}

void DrawPlusMinusComparison(TH1D *HPlusInput, TH1D *HMinusInput, const string &OutputFileName, const string &PlotLabel,
                             bool DoNormalize) {
  TH1D *HPlus = CloneHistogram(HPlusInput, string(HPlusInput->GetName()) + (DoNormalize ? "_NormWork" : "_RawWork"));
  TH1D *HMinus = CloneHistogram(HMinusInput, string(HMinusInput->GetName()) + (DoNormalize ? "_NormWork" : "_RawWork"));

  if (HPlus == nullptr || HMinus == nullptr)
    return;

  if (DoNormalize == true) {
    NormalizeHistogram(HPlus);
    NormalizeHistogram(HMinus);
  }

  StyleHistogram(HPlus, kRed + 1, "Leading HF energy (GeV)");
  StyleHistogram(HMinus, kBlue + 1, "Leading HF energy (GeV)");

  const double MeanPlus = HPlus->GetMean();
  const double MeanMinus = HMinus->GetMean();
  const double Maximum = max(HPlus->GetMaximum(), HMinus->GetMaximum());
  const double XUpper = GetDynamicXUpperEdge({HPlus, HMinus});
  const double PositiveMinimum = GetPositiveMinimum({HPlus, HMinus});
  const double Minimum = (DoNormalize == true) ? 1e-6 : max(0.5, PositiveMinimum * 0.5);
  TH1D *Ratio = BuildRatioHistogram(HPlus, HMinus, string(HPlus->GetName()) + "_Ratio");

  HPlus->SetMinimum(Minimum);
  HMinus->SetMinimum(Minimum);
  HPlus->SetMaximum(Maximum * 1.15);
  HPlus->GetXaxis()->SetRangeUser(0.0, XUpper);
  HMinus->GetXaxis()->SetRangeUser(0.0, XUpper);

  TCanvas Canvas("CanvasPlusMinus", "", 900, 820);
  TPad TopPad("TopPad", "", 0.0, 0.30, 1.0, 1.0);
  TPad BottomPad("BottomPad", "", 0.0, 0.0, 1.0, 0.30);
  TopPad.SetLeftMargin(0.12);
  TopPad.SetRightMargin(0.04);
  TopPad.SetBottomMargin(0.02);
  TopPad.SetTopMargin(0.10);
  TopPad.SetLogy();
  BottomPad.SetLeftMargin(0.12);
  BottomPad.SetRightMargin(0.04);
  BottomPad.SetBottomMargin(0.35);
  BottomPad.SetTopMargin(0.03);
  TopPad.Draw();
  BottomPad.Draw();

  TopPad.cd();
  HPlus->GetXaxis()->SetLabelSize(0);
  HPlus->GetXaxis()->SetTitleSize(0);
  HMinus->GetXaxis()->SetLabelSize(0);
  HMinus->GetXaxis()->SetTitleSize(0);
  HPlus->Draw("hist");
  HMinus->Draw("hist same");

  TLine PlusMeanLine(MeanPlus, Minimum, MeanPlus, Maximum * 1.15);
  PlusMeanLine.SetLineColor(kRed + 1);
  PlusMeanLine.SetLineStyle(2);
  PlusMeanLine.SetLineWidth(2);
  PlusMeanLine.Draw();

  TLine MinusMeanLine(MeanMinus, Minimum, MeanMinus, Maximum * 1.15);
  MinusMeanLine.SetLineColor(kBlue + 1);
  MinusMeanLine.SetLineStyle(2);
  MinusMeanLine.SetLineWidth(2);
  MinusMeanLine.Draw();

  TLegend *Legend = new TLegend(0.58, 0.76, 0.84, 0.90);
  Legend->SetBorderSize(0);
  Legend->SetLineColor(0);
  Legend->SetFillStyle(0);
  Legend->SetTextSize(0.032);
  Legend->AddEntry(HPlus, "HF+", "l");
  Legend->AddEntry(HMinus, "HF-", "l");
  Legend->Draw();

  TLatex Label;
  Label.SetNDC();
  Label.SetTextFont(42);
  Label.SetTextSize(0.035);
  Label.SetTextColor(kBlack);
  Label.DrawLatex(0.16, 0.965, Form("%s (%s)", PlotLabel.c_str(), DoNormalize ? "normalized" : "raw counts"));
  Label.SetTextColor(kRed + 1);
  Label.DrawLatex(0.58, 0.70, Form("HF+ mean = %.3f GeV", MeanPlus));
  Label.SetTextColor(kBlue + 1);
  Label.DrawLatex(0.58, 0.65, Form("HF- mean = %.3f GeV", MeanMinus));

  BottomPad.cd();
  StyleRatioHistogram(Ratio, kRed + 1, "Leading HF energy (GeV)");
  Ratio->GetXaxis()->SetRangeUser(0.0, XUpper);
  Ratio->SetMinimum(0.0);
  Ratio->SetMaximum(2.0);
  Ratio->Draw("E1");

  TLine UnityLine(0.0, 1.0, XUpper, 1.0);
  UnityLine.SetLineStyle(2);
  UnityLine.SetLineWidth(2);
  UnityLine.Draw();

  TLatex RatioLabel;
  RatioLabel.SetNDC();
  RatioLabel.SetTextFont(42);
  RatioLabel.SetTextSize(0.08);
  RatioLabel.DrawLatex(0.16, 0.82, "Denominator: HF-");

  Canvas.SaveAs(OutputFileName.c_str());
}

void DrawSideOverlap(const vector<TH1D *> &InputHistograms, const vector<string> &LegendLabels, const string &OutputFileName,
                     const string &PlotLabel, bool DoNormalize) {
  vector<TH1D *> ValidHistograms;
  for (size_t i = 0; i < InputHistograms.size(); ++i) {
    TH1D *Histogram = CloneHistogram(InputHistograms[i],
                                     string(InputHistograms[i]->GetName()) + (DoNormalize ? "_NormWork" : "_RawWork") +
                                         "_" + to_string(i));
    if (Histogram == nullptr)
      continue;
    if (DoNormalize == true)
      NormalizeHistogram(Histogram);
    ValidHistograms.push_back(Histogram);
  }

  if (ValidHistograms.empty() == true)
    return;

  const int Colors[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
  const double XUpper = GetDynamicXUpperEdge(ValidHistograms);
  const double PositiveMinimum = GetPositiveMinimum(ValidHistograms);
  double Maximum = 0.0;
  for (TH1D *Histogram : ValidHistograms)
    Maximum = max(Maximum, Histogram->GetMaximum());
  const double Minimum = (DoNormalize == true) ? 1e-6 : max(0.5, PositiveMinimum * 0.5);

  const size_t DenominatorIndex = ValidHistograms.size() - 1;
  vector<TH1D *> RatioHistograms;
  vector<size_t> RatioColorIndices;
  for (size_t i = 0; i < ValidHistograms.size(); ++i) {
    if (i == DenominatorIndex)
      continue;
    RatioHistograms.push_back(BuildRatioHistogram(ValidHistograms[i], ValidHistograms[DenominatorIndex],
                                                  string(ValidHistograms[i]->GetName()) + "_Ratio"));
    RatioColorIndices.push_back(i);
  }
  for (size_t i = 0; i < ValidHistograms.size(); ++i) {
    StyleHistogram(ValidHistograms[i], Colors[min<size_t>(i, 3)], "Leading HF energy (GeV)");
    ValidHistograms[i]->SetMinimum(Minimum);
    ValidHistograms[i]->SetMaximum(Maximum * 1.15);
    ValidHistograms[i]->GetXaxis()->SetRangeUser(0.0, XUpper);
  }

  TCanvas Canvas("CanvasSideOverlap", "", 900, 820);
  TPad TopPad("TopPadSide", "", 0.0, 0.30, 1.0, 1.0);
  TPad BottomPad("BottomPadSide", "", 0.0, 0.0, 1.0, 0.30);
  TopPad.SetLeftMargin(0.12);
  TopPad.SetRightMargin(0.04);
  TopPad.SetBottomMargin(0.02);
  TopPad.SetTopMargin(0.10);
  TopPad.SetLogy();
  BottomPad.SetLeftMargin(0.12);
  BottomPad.SetRightMargin(0.04);
  BottomPad.SetBottomMargin(0.35);
  BottomPad.SetTopMargin(0.03);
  TopPad.Draw();
  BottomPad.Draw();

  TopPad.cd();
  ValidHistograms[0]->GetXaxis()->SetLabelSize(0);
  ValidHistograms[0]->GetXaxis()->SetTitleSize(0);
  ValidHistograms[0]->Draw("hist");
  for (size_t i = 1; i < ValidHistograms.size(); ++i)
    ValidHistograms[i]->GetXaxis()->SetLabelSize(0);
  for (size_t i = 1; i < ValidHistograms.size(); ++i)
    ValidHistograms[i]->GetXaxis()->SetTitleSize(0);
  for (size_t i = 1; i < ValidHistograms.size(); ++i)
    ValidHistograms[i]->Draw("hist same");

  TLegend *Legend = new TLegend(0.42, 0.62, 0.84, 0.90);
  Legend->SetBorderSize(0);
  Legend->SetLineColor(0);
  Legend->SetFillStyle(0);
  Legend->SetTextSize(0.030);
  for (size_t i = 0; i < ValidHistograms.size() && i < LegendLabels.size(); ++i)
    Legend->AddEntry(ValidHistograms[i], LegendLabels[i].c_str(), "l");
  Legend->Draw();

  TLatex Label;
  Label.SetNDC();
  Label.SetTextFont(42);
  Label.SetTextSize(0.035);
  Label.DrawLatex(0.16, 0.965, Form("%s (%s)", PlotLabel.c_str(), DoNormalize ? "normalized" : "raw counts"));

  BottomPad.cd();
  bool FirstRatio = true;
  for (size_t i = 0; i < RatioHistograms.size(); ++i) {
    if (RatioHistograms[i] == nullptr)
      continue;

    StyleRatioHistogram(RatioHistograms[i], Colors[min<size_t>(RatioColorIndices[i], 3)], "Leading HF energy (GeV)");
    RatioHistograms[i]->GetXaxis()->SetRangeUser(0.0, XUpper);
    RatioHistograms[i]->SetMinimum(0.0);
    RatioHistograms[i]->SetMaximum(2.0);
    RatioHistograms[i]->Draw(FirstRatio ? "E1" : "E1 same");
    FirstRatio = false;
  }

  TLine UnityLine(0.0, 1.0, XUpper, 1.0);
  UnityLine.SetLineStyle(2);
  UnityLine.SetLineWidth(2);
  UnityLine.Draw();

  TLatex RatioLabel;
  RatioLabel.SetNDC();
  RatioLabel.SetTextFont(42);
  RatioLabel.SetTextSize(0.08);
  RatioLabel.DrawLatex(0.16, 0.82, Form("Denominator: %s", LegendLabels[DenominatorIndex].c_str()));

  Canvas.SaveAs(OutputFileName.c_str());
}
} // namespace

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  const string InputFileName = CL.Get("Input");
  const string OutputDirectory = CL.Get("Output");

  if (InputFileName.empty() == true || OutputDirectory.empty() == true) {
    cerr << "Usage: " << argv[0] << " --Input <input.root> --Output <output_directory>" << endl;
    return 1;
  }

  TFile InputFile(InputFileName.c_str(), "READ");
  if (InputFile.IsZombie() == true) {
    cerr << "Failed to open input file " << InputFileName << endl;
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(OutputDirectory.c_str(), true);

  struct VariableEntry {
    string Suffix;
    string Label;
  };
  const vector<VariableEntry> Variables = {
      {"Xn0nAll", "ZB + at least one side ADC <=19"},
      {"Xn0nOneSide", "ZB + side with ADC<=19 filled"},
      {"0n0n", "ZB+ both ADC<=19"},
      {"NoCut", "ZB events"},
  };

  vector<TH1D *> PlusHistograms;
  vector<TH1D *> MinusHistograms;
  vector<string> SideLegendLabels;

  for (const VariableEntry &Entry : Variables) {
    TH1D *HPlus = LoadHistogram(InputFile, "hHFPlusLeading_" + Entry.Suffix, "hHFPlusLeading_" + Entry.Suffix + "Clone");
    TH1D *HMinus = LoadHistogram(InputFile, "hHFMinusLeading_" + Entry.Suffix, "hHFMinusLeading_" + Entry.Suffix + "Clone");

    if (HPlus == nullptr || HMinus == nullptr) {
      cerr << "Failed to load hHFPlusLeading_" << Entry.Suffix << " or hHFMinusLeading_" << Entry.Suffix << endl;
      return 1;
    }

    DrawPlusMinusComparison(HPlus, HMinus, OutputDirectory + "/HF2Leading_" + Entry.Suffix + "_PlusMinus.pdf",
                            Entry.Label, true);
    DrawPlusMinusComparison(HPlus, HMinus, OutputDirectory + "/HF2Leading_" + Entry.Suffix + "_PlusMinus_Raw.pdf",
                            Entry.Label, false);

    PlusHistograms.push_back(HPlus);
    MinusHistograms.push_back(HMinus);
    SideLegendLabels.push_back(Entry.Label);
  }

  DrawSideOverlap(PlusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Plus_AllCategories.pdf",
                  "HF+ category overlap", true);
  DrawSideOverlap(PlusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Plus_AllCategories_Raw.pdf",
                  "HF+ category overlap", false);
  DrawSideOverlap(MinusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Minus_AllCategories.pdf",
                  "HF- category overlap", true);
  DrawSideOverlap(MinusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Minus_AllCategories_Raw.pdf",
                  "HF- category overlap", false);

  return 0;
}
