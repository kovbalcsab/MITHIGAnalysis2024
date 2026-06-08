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
  Histogram->GetXaxis()->SetTitleSize(0.060);
  Histogram->GetXaxis()->SetLabelSize(0.050);
  Histogram->GetYaxis()->SetTitleSize(0.060);
  Histogram->GetYaxis()->SetLabelSize(0.050);
  Histogram->GetYaxis()->SetTitleOffset(0.9);
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
  Histogram->GetYaxis()->SetTitleSize(0.12);
  Histogram->GetYaxis()->SetTitleOffset(0.45);
  Histogram->GetYaxis()->SetLabelSize(0.10);
  Histogram->GetYaxis()->SetNdivisions(505);
  Histogram->GetXaxis()->SetTitleSize(0.12);
  Histogram->GetXaxis()->SetTitleOffset(1.0);
  Histogram->GetXaxis()->SetLabelSize(0.10);
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

TH1D *BuildRatioHistogram(TH1D *Numerator, TH1D *Denominator, const string &Name, bool UseBinomial = false) {
  if (Numerator == nullptr || Denominator == nullptr)
    return nullptr;

  TH1D *Ratio = dynamic_cast<TH1D *>(Numerator->Clone(Name.c_str()));
  if (Ratio == nullptr)
    return nullptr;

  Ratio->SetDirectory(nullptr);
  Ratio->Sumw2();
  Ratio->Divide(Numerator, Denominator, 1.0, 1.0, UseBinomial == true ? "B" : "");
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

pair<double, double> GetRequestedRatioBounds(double RatioDistanceAroundUnity) {
  if (RatioDistanceAroundUnity > 0.0)
    return {1.0 - RatioDistanceAroundUnity, 1.0 + RatioDistanceAroundUnity};
  return {0.0, 2.0};
}

void DrawPlusMinusComparison(TH1D *HPlusInput, TH1D *HMinusInput, const string &OutputFileName, const string &PlotLabel,
                             bool DoNormalize, double RatioDistanceAroundUnity) {
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
  Legend->SetTextSize(0.045);
  Legend->AddEntry(HPlus, "HF+", "l");
  Legend->AddEntry(HMinus, "HF-", "l");
  Legend->Draw();

  TLatex Label;
  Label.SetNDC();
  Label.SetTextFont(42);
  Label.SetTextSize(0.040);
  Label.SetTextColor(kBlack);
  Label.DrawLatex(0.16, 0.965, Form("%s (%s)", PlotLabel.c_str(), DoNormalize ? "normalized" : "raw counts"));
  Label.SetTextColor(kRed + 1);
  Label.DrawLatex(0.58, 0.70, Form("HF+ mean = %.3f GeV", MeanPlus));
  Label.SetTextColor(kBlue + 1);
  Label.DrawLatex(0.58, 0.65, Form("HF- mean = %.3f GeV", MeanMinus));

  BottomPad.cd();
  StyleRatioHistogram(Ratio, kRed + 1, "Leading HF energy (GeV)");
  Ratio->GetXaxis()->SetRangeUser(0.0, XUpper);
  const pair<double, double> RatioBounds = GetRequestedRatioBounds(RatioDistanceAroundUnity);
  Ratio->SetMinimum(RatioBounds.first);
  Ratio->SetMaximum(RatioBounds.second);
  Ratio->Draw("E1P");

  TLine UnityLine(0.0, 1.0, XUpper, 1.0);
  UnityLine.SetLineStyle(2);
  UnityLine.SetLineWidth(2);
  UnityLine.Draw();

  TLatex RatioLabel;
  RatioLabel.SetNDC();
  RatioLabel.SetTextFont(42);
  RatioLabel.SetTextSize(0.09);
  RatioLabel.DrawLatex(0.16, 0.82, "Denominator: HF-");

  Canvas.SaveAs(OutputFileName.c_str());
}

void DrawDefaultTrigComparison(TH1D *HDefaultInput, TH1D *HTrigInput, const string &OutputFileName,
                               const string &PlotLabel, bool DoNormalize, double RatioDistanceAroundUnity) {
  TH1D *HDefault = CloneHistogram(HDefaultInput, string(HDefaultInput->GetName()) +
                                                     (DoNormalize ? "_OfflineNormWork" : "_OfflineRawWork"));
  TH1D *HTrig =
      CloneHistogram(HTrigInput, string(HTrigInput->GetName()) + (DoNormalize ? "_TrigNormWork" : "_TrigRawWork"));

  if (HDefault == nullptr || HTrig == nullptr)
    return;

  if (DoNormalize == true) {
    NormalizeHistogram(HDefault);
    NormalizeHistogram(HTrig);
  }

  StyleHistogram(HDefault, kBlack, "Leading HF energy (GeV)");
  StyleHistogram(HTrig, kRed + 1, "Leading HF energy (GeV)");

  const double MeanDefault = HDefault->GetMean();
  const double MeanTrig = HTrig->GetMean();
  const double Maximum = max(HDefault->GetMaximum(), HTrig->GetMaximum());
  const double XUpper = GetDynamicXUpperEdge({HDefault, HTrig});
  const double PositiveMinimum = GetPositiveMinimum({HDefault, HTrig});
  const double Minimum = (DoNormalize == true) ? 1e-6 : max(0.5, PositiveMinimum * 0.5);
  TH1D *Ratio = BuildRatioHistogram(HTrig, HDefault, string(HTrig->GetName()) + "_Ratio", true);

  HDefault->SetMinimum(Minimum);
  HTrig->SetMinimum(Minimum);
  HDefault->SetMaximum(Maximum * 1.15);
  HDefault->GetXaxis()->SetRangeUser(0.0, XUpper);
  HTrig->GetXaxis()->SetRangeUser(0.0, XUpper);

  TCanvas Canvas("CanvasDefaultTrig", "", 900, 820);
  TPad TopPad("TopPadDefaultTrig", "", 0.0, 0.30, 1.0, 1.0);
  TPad BottomPad("BottomPadDefaultTrig", "", 0.0, 0.0, 1.0, 0.30);
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
  HDefault->GetXaxis()->SetLabelSize(0);
  HDefault->GetXaxis()->SetTitleSize(0);
  HTrig->GetXaxis()->SetLabelSize(0);
  HTrig->GetXaxis()->SetTitleSize(0);
  HDefault->Draw("hist");
  HTrig->Draw("hist same");

  TLine DefaultMeanLine(MeanDefault, Minimum, MeanDefault, Maximum * 1.15);
  DefaultMeanLine.SetLineColor(kBlack);
  DefaultMeanLine.SetLineStyle(2);
  DefaultMeanLine.SetLineWidth(2);
  DefaultMeanLine.Draw();

  TLine TrigMeanLine(MeanTrig, Minimum, MeanTrig, Maximum * 1.15);
  TrigMeanLine.SetLineColor(kRed + 1);
  TrigMeanLine.SetLineStyle(2);
  TrigMeanLine.SetLineWidth(2);
  TrigMeanLine.Draw();

  TLegend *Legend = new TLegend(0.56, 0.76, 0.84, 0.90);
  Legend->SetBorderSize(0);
  Legend->SetLineColor(0);
  Legend->SetFillStyle(0);
  Legend->SetTextSize(0.045);
  Legend->AddEntry(HDefault, "Offline cuts", "l");
  Legend->AddEntry(HTrig, "Offline cuts + ADC <= 19", "l");
  Legend->Draw();

  TLatex Label;
  Label.SetNDC();
  Label.SetTextFont(42);
  Label.SetTextSize(0.040);
  Label.SetTextColor(kBlack);
  Label.DrawLatex(0.16, 0.965, Form("%s (%s)", PlotLabel.c_str(), DoNormalize ? "normalized" : "raw counts"));
  Label.SetTextColor(kBlack);
  Label.DrawLatex(0.56, 0.70, Form("Offline mean = %.3f GeV", MeanDefault));
  Label.SetTextColor(kRed + 1);
  Label.DrawLatex(0.56, 0.65, Form("Trig mean = %.3f GeV", MeanTrig));

  BottomPad.cd();
  StyleRatioHistogram(Ratio, kRed + 1, "Leading HF energy (GeV)");
  Ratio->GetXaxis()->SetRangeUser(0.0, XUpper);
  const pair<double, double> RatioBounds = GetRequestedRatioBounds(RatioDistanceAroundUnity);
  Ratio->SetMinimum(RatioBounds.first);
  Ratio->SetMaximum(RatioBounds.second);
  Ratio->Draw("E1P");

  TLine UnityLine(0.0, 1.0, XUpper, 1.0);
  UnityLine.SetLineStyle(2);
  UnityLine.SetLineWidth(2);
  UnityLine.Draw();

  TLatex RatioLabel;
  RatioLabel.SetNDC();
  RatioLabel.SetTextFont(42);
  RatioLabel.SetTextSize(0.09);
  RatioLabel.DrawLatex(0.16, 0.82, "Denominator: Offline cuts");

  Canvas.SaveAs(OutputFileName.c_str());
}

void DrawSideOverlap(const vector<TH1D *> &InputHistograms, const vector<string> &LegendLabels,
                     const string &OutputFileName, const string &PlotLabel, bool DoNormalize,
                     double RatioDistanceAroundUnity, const string &XTitle, double FixedXUpper = -1.0) {
  vector<TH1D *> ValidHistograms;
  for (size_t i = 0; i < InputHistograms.size(); ++i) {
    TH1D *Histogram =
        CloneHistogram(InputHistograms[i], string(InputHistograms[i]->GetName()) +
                                               (DoNormalize ? "_NormWork" : "_RawWork") + "_" + to_string(i));
    if (Histogram == nullptr)
      continue;
    if (DoNormalize == true)
      NormalizeHistogram(Histogram);
    ValidHistograms.push_back(Histogram);
  }

  if (ValidHistograms.empty() == true)
    return;

  const int Colors[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
  const double XUpper = (FixedXUpper > 0.0) ? FixedXUpper : GetDynamicXUpperEdge(ValidHistograms);
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
                                                  string(ValidHistograms[i]->GetName()) + "_Ratio", true));
    RatioColorIndices.push_back(i);
  }
  for (size_t i = 0; i < ValidHistograms.size(); ++i) {
    StyleHistogram(ValidHistograms[i], Colors[min<size_t>(i, 3)], XTitle);
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
  Legend->SetTextSize(0.040);
  for (size_t i = 0; i < ValidHistograms.size() && i < LegendLabels.size(); ++i)
    Legend->AddEntry(ValidHistograms[i], LegendLabels[i].c_str(), "l");
  Legend->Draw();

  TLatex Label;
  Label.SetNDC();
  Label.SetTextFont(42);
  Label.SetTextSize(0.040);
  Label.DrawLatex(0.16, 0.965, Form("%s (%s)", PlotLabel.c_str(), DoNormalize ? "normalized" : "raw counts"));

  BottomPad.cd();
  const pair<double, double> RatioBounds = GetRequestedRatioBounds(RatioDistanceAroundUnity);
  bool FirstRatio = true;
  for (size_t i = 0; i < RatioHistograms.size(); ++i) {
    if (RatioHistograms[i] == nullptr)
      continue;

    StyleRatioHistogram(RatioHistograms[i], Colors[min<size_t>(RatioColorIndices[i], 3)], XTitle);
    RatioHistograms[i]->GetXaxis()->SetRangeUser(0.0, XUpper);
    RatioHistograms[i]->SetMinimum(RatioBounds.first);
    RatioHistograms[i]->SetMaximum(RatioBounds.second);
    RatioHistograms[i]->Draw(FirstRatio ? "E1P" : "E1P same");
    FirstRatio = false;
  }

  TLine UnityLine(0.0, 1.0, XUpper, 1.0);
  UnityLine.SetLineStyle(2);
  UnityLine.SetLineWidth(2);
  UnityLine.Draw();

  TLatex RatioLabel;
  RatioLabel.SetNDC();
  RatioLabel.SetTextFont(42);
  RatioLabel.SetTextSize(0.09);
  RatioLabel.DrawLatex(0.16, 0.82, Form("Denominator: %s", LegendLabels[DenominatorIndex].c_str()));

  Canvas.SaveAs(OutputFileName.c_str());
}
} // namespace

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  const string InputFileName = CL.Get("Input");
  const string OutputDirectory = CL.Get("Output");
  const double RatioDistanceAroundUnity = CL.GetDouble("RatioDistanceAroundUnity", -1.0);

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
      {"Xn0nAll", "Trigger + at least one side ADC <=19"},
      {"Xn0nOneSide", "Trigger + side with ADC<=19 filled"},
      {"0n0n", "Trigger + both ADC<=19"},
      {"NoCut", "Trigger events"},
  };

  vector<TH1D *> PlusHistograms;
  vector<TH1D *> MinusHistograms;
  vector<string> SideLegendLabels;

  for (const VariableEntry &Entry : Variables) {
    TH1D *HPlus =
        LoadHistogram(InputFile, "hHFPlusLeading_" + Entry.Suffix, "hHFPlusLeading_" + Entry.Suffix + "Clone");
    TH1D *HMinus =
        LoadHistogram(InputFile, "hHFMinusLeading_" + Entry.Suffix, "hHFMinusLeading_" + Entry.Suffix + "Clone");

    if (HPlus == nullptr || HMinus == nullptr) {
      cerr << "Failed to load hHFPlusLeading_" << Entry.Suffix << " or hHFMinusLeading_" << Entry.Suffix << endl;
      return 1;
    }

    DrawPlusMinusComparison(HPlus, HMinus, OutputDirectory + "/HF2Leading_" + Entry.Suffix + "_PlusMinus.pdf",
                            Entry.Label, true, RatioDistanceAroundUnity);
    DrawPlusMinusComparison(HPlus, HMinus, OutputDirectory + "/HF2Leading_" + Entry.Suffix + "_PlusMinus_Raw.pdf",
                            Entry.Label, false, RatioDistanceAroundUnity);

    PlusHistograms.push_back(HPlus);
    MinusHistograms.push_back(HMinus);
    SideLegendLabels.push_back(Entry.Label);
  }

  DrawSideOverlap(PlusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Plus_AllCategories.pdf",
                  "HF+ category overlap", true, RatioDistanceAroundUnity, "Leading HF energy (GeV)");
  DrawSideOverlap(PlusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Plus_AllCategories_Raw.pdf",
                  "HF+ category overlap", false, RatioDistanceAroundUnity, "Leading HF energy (GeV)");
  DrawSideOverlap(MinusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Minus_AllCategories.pdf",
                  "HF- category overlap", true, RatioDistanceAroundUnity, "Leading HF energy (GeV)");
  DrawSideOverlap(MinusHistograms, SideLegendLabels, OutputDirectory + "/HF2Leading_Minus_AllCategories_Raw.pdf",
                  "HF- category overlap", false, RatioDistanceAroundUnity, "Leading HF energy (GeV)");

  TH1D *H0nXnPlusDefault =
      LoadHistogram(InputFile, "hHFPlusLeading_0nXn_ZB_OfflineSel", "hHFPlusLeading_0nXn_ZB_OfflineSelClone");
  TH1D *H0nXnPlusTrig =
      LoadHistogram(InputFile, "hHFPlusLeading_0nXn_ZB_OfflineSel_Trig", "hHFPlusLeading_0nXn_ZB_OfflineSel_TrigClone");
  TH1D *H0nXnMinusDefault =
      LoadHistogram(InputFile, "hHFMinusLeading_0nXn_ZB_OfflineSel", "hHFMinusLeading_0nXn_ZB_OfflineSelClone");
  TH1D *H0nXnMinusTrig = LoadHistogram(InputFile, "hHFMinusLeading_0nXn_ZB_OfflineSel_Trig",
                                       "hHFMinusLeading_0nXn_ZB_OfflineSel_TrigClone");
  TH1D *H0nAnMinusDefault =
      LoadHistogram(InputFile, "hHFMinusLeading_0nAn_ZB_OfflineSel", "hHFMinusLeading_0nAn_ZB_OfflineSelClone");
  TH1D *H0nAnMinusTrig = LoadHistogram(InputFile, "hHFMinusLeading_0nAn_ZB_OfflineSel_Trig",
                                       "hHFMinusLeading_0nAn_ZB_OfflineSel_TrigClone");
  TH1D *H0nAnPlusDefault =
      LoadHistogram(InputFile, "hHFPlusLeading_0nAn_ZB_OfflineSel", "hHFPlusLeading_0nAn_ZB_OfflineSelClone");
  TH1D *H0nAnPlusTrig =
      LoadHistogram(InputFile, "hHFPlusLeading_0nAn_ZB_OfflineSel_Trig", "hHFPlusLeading_0nAn_ZB_OfflineSel_TrigClone");
  TH1D *H0n0nMinusDefault =
      LoadHistogram(InputFile, "hHFMinusLeading_0n0n_ZB_OfflineSel", "hHFMinusLeading_0n0n_ZB_OfflineSelClone");
  TH1D *H0n0nMinusTrig = LoadHistogram(InputFile, "hHFMinusLeading_0n0n_ZB_OfflineSel_Trig",
                                       "hHFMinusLeading_0n0n_ZB_OfflineSel_TrigClone");
  TH1D *H0n0nPlusDefault =
      LoadHistogram(InputFile, "hHFPlusLeading_0n0n_ZB_OfflineSel", "hHFPlusLeading_0n0n_ZB_OfflineSelClone");
  TH1D *H0n0nPlusTrig =
      LoadHistogram(InputFile, "hHFPlusLeading_0n0n_ZB_OfflineSel_Trig", "hHFPlusLeading_0n0n_ZB_OfflineSel_TrigClone");
  TH1D *HMultNoCut = LoadHistogram(InputFile, "hMult_NoCut", "hMult_NoCutClone");
  TH1D *HMult0nAn = LoadHistogram(InputFile, "hMult_0nAn", "hMult_0nAnClone");
  TH1D *HMult0nXn = LoadHistogram(InputFile, "hMult_0nXn", "hMult_0nXnClone");
  TH1D *HMult0n0n = LoadHistogram(InputFile, "hMult_0n0n", "hMult_0n0nClone");

  if (H0nXnPlusDefault == nullptr || H0nXnPlusTrig == nullptr || H0nXnMinusDefault == nullptr ||
      H0nXnMinusTrig == nullptr || H0nAnMinusDefault == nullptr || H0nAnMinusTrig == nullptr ||
      H0nAnPlusDefault == nullptr || H0nAnPlusTrig == nullptr || H0n0nMinusDefault == nullptr ||
      H0n0nMinusTrig == nullptr || H0n0nPlusDefault == nullptr || H0n0nPlusTrig == nullptr || HMultNoCut == nullptr ||
      HMult0nAn == nullptr || HMult0nXn == nullptr || HMult0n0n == nullptr) {
    cerr << "Failed to load one or more 0nXn/0nAn/0n0n or multiplicity histograms" << endl;
    return 1;
  }

  DrawDefaultTrigComparison(H0nXnPlusDefault, H0nXnPlusTrig,
                            OutputDirectory + "/HF2Leading_0nXn_HFPlus_OfflineVsTrig.pdf", "0nXn HF+", true,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nXnPlusDefault, H0nXnPlusTrig,
                            OutputDirectory + "/HF2Leading_0nXn_HFPlus_OfflineVsTrig_Raw.pdf", "0nXn HF+", false,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nXnMinusDefault, H0nXnMinusTrig,
                            OutputDirectory + "/HF2Leading_0nXn_HFMinus_OfflineVsTrig.pdf", "0nXn HF-", true,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nXnMinusDefault, H0nXnMinusTrig,
                            OutputDirectory + "/HF2Leading_0nXn_HFMinus_OfflineVsTrig_Raw.pdf", "0nXn HF-", false,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nAnMinusDefault, H0nAnMinusTrig,
                            OutputDirectory + "/HF2Leading_0nAn_HFMinus_OfflineVsTrig.pdf", "0nAn HF-", true,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nAnMinusDefault, H0nAnMinusTrig,
                            OutputDirectory + "/HF2Leading_0nAn_HFMinus_OfflineVsTrig_Raw.pdf", "0nAn HF-", false,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nAnPlusDefault, H0nAnPlusTrig,
                            OutputDirectory + "/HF2Leading_0nAn_HFPlus_OfflineVsTrig.pdf", "0nAn HF+", true,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0nAnPlusDefault, H0nAnPlusTrig,
                            OutputDirectory + "/HF2Leading_0nAn_HFPlus_OfflineVsTrig_Raw.pdf", "0nAn HF+", false,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0n0nMinusDefault, H0n0nMinusTrig,
                            OutputDirectory + "/HF2Leading_0n0n_HFMinus_OfflineVsTrig.pdf", "0n0n HF-", true,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0n0nMinusDefault, H0n0nMinusTrig,
                            OutputDirectory + "/HF2Leading_0n0n_HFMinus_OfflineVsTrig_Raw.pdf", "0n0n HF-", false,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0n0nPlusDefault, H0n0nPlusTrig,
                            OutputDirectory + "/HF2Leading_0n0n_HFPlus_OfflineVsTrig.pdf", "0n0n HF+", true,
                            RatioDistanceAroundUnity);
  DrawDefaultTrigComparison(H0n0nPlusDefault, H0n0nPlusTrig,
                            OutputDirectory + "/HF2Leading_0n0n_HFPlus_OfflineVsTrig_Raw.pdf", "0n0n HF+", false,
                            RatioDistanceAroundUnity);

  const vector<TH1D *> MultiplicityHistograms = {HMult0nAn, HMult0nXn, HMult0n0n, HMultNoCut};
  const vector<string> MultiplicityLabels = {"0nAn", "0nXn", "0n0n", "No cut"};
  DrawSideOverlap(MultiplicityHistograms, MultiplicityLabels, OutputDirectory + "/HF2Multiplicity_AllCategories.pdf",
                  "Track multiplicity overlap", true, RatioDistanceAroundUnity, "Track multiplicity", 300.0);
  DrawSideOverlap(MultiplicityHistograms, MultiplicityLabels,
                  OutputDirectory + "/HF2Multiplicity_AllCategories_Raw.pdf", "Track multiplicity overlap", false,
                  RatioDistanceAroundUnity, "Track multiplicity", 300.0);

  return 0;
}
