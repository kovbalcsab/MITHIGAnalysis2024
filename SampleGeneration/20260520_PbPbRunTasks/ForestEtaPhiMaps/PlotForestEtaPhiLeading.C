#include <iostream>
#include <string>

#include "RtypesCore.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
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

void StyleHistogram(TH1D *Histogram, int Color) {
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
  Histogram->GetXaxis()->SetTitle("Leading HF energy (GeV)");
  Histogram->GetYaxis()->SetTitle("Entries");
  Histogram->GetYaxis()->SetTitleOffset(1.3);
}
} // namespace

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  const string InputFileName = CL.Get("Input");
  const string OutputFileName = CL.Get("Output");

  if (InputFileName.empty() == true || OutputFileName.empty() == true) {
    cerr << "Usage: " << argv[0] << " --Input <input.root> --Output <output.pdf>" << endl;
    return 1;
  }

  TFile InputFile(InputFileName.c_str(), "READ");
  if (InputFile.IsZombie() == true) {
    cerr << "Failed to open input file " << InputFileName << endl;
    return 1;
  }

  TH1D *HPlus = dynamic_cast<TH1D *>(InputFile.Get("hHFEMaxPlusLeading"));
  TH1D *HMinus = dynamic_cast<TH1D *>(InputFile.Get("hHFEMaxMinusLeading"));

  if (HPlus == nullptr || HMinus == nullptr) {
    cerr << "Failed to load hHFEMaxPlusLeading or hHFEMaxMinusLeading from " << InputFileName << endl;
    return 1;
  }

  HPlus = dynamic_cast<TH1D *>(HPlus->Clone("hHFEMaxPlusLeadingClone"));
  HMinus = dynamic_cast<TH1D *>(HMinus->Clone("hHFEMaxMinusLeadingClone"));
  HPlus->SetDirectory(nullptr);
  HMinus->SetDirectory(nullptr);

  NormalizeHistogram(HPlus);
  NormalizeHistogram(HMinus);

  TH1D *hPlusCum = HPlus ? dynamic_cast<TH1D *>(HPlus->GetCumulative(kTRUE)) : nullptr;
  TH1D *hMinusCum = HMinus ? dynamic_cast<TH1D *>(HMinus->GetCumulative(kTRUE)) : nullptr;
  double plusProb = 0;
  double minusProb = 0;
  double sameBin = 0;
  for (int iBin = 1; iBin <= HPlus->GetNbinsX() + 1; iBin++) {
    plusProb += HPlus->GetBinContent(iBin) * hMinusCum->GetBinContent(iBin - 1);
    sameBin += HPlus->GetBinContent(iBin) * HMinus->GetBinContent(iBin);
  }
  for (int iBin = 1; iBin <= HMinus->GetNbinsX() + 1; iBin++) {
    minusProb += HMinus->GetBinContent(iBin) * hPlusCum->GetBinContent(iBin - 1);
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const string OutputDirectory = gSystem->DirName(OutputFileName.c_str());
  if (OutputDirectory.empty() == false && OutputDirectory != ".")
    gSystem->mkdir(OutputDirectory.c_str(), true);

  StyleHistogram(HPlus, kRed + 1);
  StyleHistogram(HMinus, kBlue + 1);

  const double MeanPlus = HPlus->GetMean();
  const double MeanMinus = HMinus->GetMean();
  const double Maximum = max(HPlus->GetMaximum(), HMinus->GetMaximum());
  const double Minimum = 1e-6;
  HPlus->SetMinimum(1e-6);
  HPlus->SetMaximum(Maximum * 1.15);
  HPlus->GetXaxis()->SetRangeUser(0.0, 30.0);

  cout << "HF+ mean: " << MeanPlus << endl;
  cout << "HF- mean: " << MeanMinus << endl;
  cout << "HF+ probability of being more than HF-: " << plusProb << endl;
  cout << "HF- probability of being more than HF+: " << minusProb << endl;
  cout << "Probability of being in the same bin: " << sameBin << endl;
  cout << "Total probability: " << plusProb + minusProb + sameBin << endl;

  TCanvas Canvas("Canvas", "", 900, 700);
  Canvas.cd();
  Canvas.SetLeftMargin(0.12);
  Canvas.SetRightMargin(0.04);
  Canvas.SetBottomMargin(0.12);
  Canvas.SetTopMargin(0.06);
  Canvas.SetLogy();

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

  TLatex Label;
  Label.SetNDC();
  Label.SetTextFont(42);
  Label.SetTextSize(0.035);
  Label.SetTextColor(kRed + 1);
  Label.DrawLatex(0.68, 0.72, Form("HF+ mean = %.3f GeV", MeanPlus));
  Label.SetTextColor(kBlue + 1);
  Label.DrawLatex(0.68, 0.67, Form("HF- mean = %.3f GeV", MeanMinus));

  TLegend Legend(0.68, 0.78, 0.88, 0.90);
  Legend.SetBorderSize(0);
  Legend.SetFillStyle(0);
  Legend.AddEntry(HPlus, "HF+", "l");
  Legend.AddEntry(HMinus, "HF-", "l");
  Legend.Draw();

  Canvas.SaveAs(OutputFileName.c_str());

  delete HPlus;
  delete HMinus;

  return 0;
}
