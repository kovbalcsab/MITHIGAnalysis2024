#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "CommandLine.h"

using namespace std;

////////////////////////////////////////////////////////////
// CMS style (minimal clean version)
////////////////////////////////////////////////////////////
void SetCMSStyle() {
  gStyle->SetOptStat(0);

  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetFrameLineWidth(1);

  gStyle->SetEndErrorSize(0);
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
  pad->SetTopMargin(0.05);
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

void DrawLeadingOverlap(TH1D *hist1, TH1D *hist2, const string &outputFileName, int firstColor = kBlue,
                        int secondColor = kRed) {
  if (hist1 == nullptr)
    return;

  TCanvas canvas("cLeading", "", 900, 700);
  canvas.cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.06);

  DrawHistogram(hist1, true, true, firstColor, "hist");
  if (hist2 != nullptr)
    DrawHistogram(hist2, true, true, secondColor, "hist same");

  canvas.SaveAs(outputFileName.c_str());
}

double GetDynamicXUpperEdge(TH1D *hist1, TH1D *hist2) {
  double maximumX = 0.0;

  for (TH1D *hist : {hist1, hist2}) {
    if (hist == nullptr)
      continue;

    for (int i = hist->GetNbinsX(); i >= 1; --i) {
      if (hist->GetBinContent(i) == 0.0)
        continue;

      maximumX = max(maximumX, hist->GetXaxis()->GetBinUpEdge(i));
      break;
    }
  }

  if (maximumX <= 0.0) {
    if (hist1 != nullptr)
      maximumX = hist1->GetXaxis()->GetBinUpEdge(hist1->GetNbinsX());
    else if (hist2 != nullptr)
      maximumX = hist2->GetXaxis()->GetBinUpEdge(hist2->GetNbinsX());
  }

  return maximumX * 1.1;
}

void DrawHFnPFOverlap(TH1D *hist1, TH1D *hist2, const string &outputFileName, int firstColor = kBlue,
                      int secondColor = kRed) {
  if (hist1 == nullptr)
    return;

  NormalizeHistogram(hist1);
  NormalizeHistogram(hist2);
  const double xUpperEdge = GetDynamicXUpperEdge(hist1, hist2);

  hist1->SetLineColor(firstColor);
  hist1->SetLineWidth(2);
  hist1->SetMarkerColor(firstColor);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(0.8);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitle("HF nPF");
  hist1->GetYaxis()->SetTitle("Entries");
  hist1->GetYaxis()->SetTitleOffset(1.3);
  hist1->SetMinimum(1e-5);
  hist1->GetXaxis()->SetRangeUser(0.0, xUpperEdge);

  if (hist2 != nullptr) {
    hist2->SetLineColor(secondColor);
    hist2->SetLineWidth(2);
    hist2->SetMarkerColor(secondColor);
    hist2->SetMarkerStyle(20);
    hist2->SetMarkerSize(0.8);
    hist2->SetTitle("");
    hist2->GetXaxis()->SetTitle("HF nPF");
    hist2->GetYaxis()->SetTitle("Entries");
    hist2->GetYaxis()->SetTitleOffset(1.3);
    hist2->SetMinimum(1e-5);
    hist2->GetXaxis()->SetRangeUser(0.0, xUpperEdge);
  }

  TCanvas canvas("cHFnPF", "", 900, 700);
  canvas.cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.06);

  hist1->Draw("hist");
  if (hist2 != nullptr)
    hist2->Draw("hist same");

  canvas.SaveAs(outputFileName.c_str());
}

void Draw1DOverlap(TH1D *hist1, TH1D *hist2, const string &outputFileName, const string &xTitle, int firstColor = kBlue,
                   int secondColor = kRed) {
  if (hist1 == nullptr)
    return;

  NormalizeHistogram(hist1);
  NormalizeHistogram(hist2);
  const double xUpperEdge = GetDynamicXUpperEdge(hist1, hist2);

  hist1->SetLineColor(firstColor);
  hist1->SetLineWidth(2);
  hist1->SetMarkerColor(firstColor);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(0.8);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitle(xTitle.c_str());
  hist1->GetYaxis()->SetTitle("Entries");
  hist1->GetYaxis()->SetTitleOffset(1.3);
  hist1->SetMinimum(1e-5);
  hist1->GetXaxis()->SetRangeUser(0.0, xUpperEdge);

  if (hist2 != nullptr) {
    hist2->SetLineColor(secondColor);
    hist2->SetLineWidth(2);
    hist2->SetMarkerColor(secondColor);
    hist2->SetMarkerStyle(20);
    hist2->SetMarkerSize(0.8);
    hist2->SetTitle("");
    hist2->GetXaxis()->SetTitle(xTitle.c_str());
    hist2->GetYaxis()->SetTitle("Entries");
    hist2->GetYaxis()->SetTitleOffset(1.3);
    hist2->SetMinimum(1e-5);
    hist2->GetXaxis()->SetRangeUser(0.0, xUpperEdge);
  }

  TCanvas canvas("c1DOverlap", "", 900, 700);
  canvas.cd();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.06);

  hist1->Draw("hist");
  if (hist2 != nullptr)
    hist2->Draw("hist same");

  canvas.SaveAs(outputFileName.c_str());
}

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);

  string InputFileName = CL.Get("Input");
  string InputFileName2 = CL.Get("Input2", "");
  string OutputFileName = CL.Get("Output");
  bool doCumulative = CL.GetBool("DoCumulative", false);

  SetCMSStyle();

  TFile *InputFile = TFile::Open(InputFileName.c_str());
  TFile *InputFile2 = (InputFileName2.empty() == false) ? TFile::Open(InputFileName2.c_str()) : nullptr;
  const bool hasSecondInput = (InputFile2 != nullptr);

  TH1D *hPlusLeading = (TH1D *)InputFile->Get("hHFEMaxPlusLeading");
  TH1D *hMinusLeading = (TH1D *)InputFile->Get("hHFEMaxMinusLeading");
  TH1D *hPlusLeading2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFEMaxPlusLeading") : nullptr;
  TH1D *hMinusLeading2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFEMaxMinusLeading") : nullptr;
  TH1D *hHFnPF = (TH1D *)InputFile->Get("hHFnPF");
  TH1D *hHFnPF2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFnPF") : nullptr;
  TH1D *hHFnPF6p = (TH1D *)InputFile->Get("hHFnPF6p");
  TH1D *hHFnPF7p = (TH1D *)InputFile->Get("hHFnPF7p");
  TH1D *hHFnPF6m = (TH1D *)InputFile->Get("hHFnPF6m");
  TH1D *hHFnPF7m = (TH1D *)InputFile->Get("hHFnPF7m");
  TH1D *hHFnPF6p2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFnPF6p") : nullptr;
  TH1D *hHFnPF7p2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFnPF7p") : nullptr;
  TH1D *hHFnPF6m2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFnPF6m") : nullptr;
  TH1D *hHFnPF7m2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hHFnPF7m") : nullptr;
  TH1D *hAllEnergyp = (TH1D *)InputFile->Get("hAllEnergyp");
  TH1D *hAllEnergym = (TH1D *)InputFile->Get("hAllEnergym");
  TH1D *hAllEnergy6p = (TH1D *)InputFile->Get("hAllEnergy6p");
  TH1D *hAllEnergy7p = (TH1D *)InputFile->Get("hAllEnergy7p");
  TH1D *hAllEnergy6m = (TH1D *)InputFile->Get("hAllEnergy6m");
  TH1D *hAllEnergy7m = (TH1D *)InputFile->Get("hAllEnergy7m");
  TH1D *hAllEnergyp2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hAllEnergyp") : nullptr;
  TH1D *hAllEnergym2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hAllEnergym") : nullptr;
  TH1D *hAllEnergy6p2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hAllEnergy6p") : nullptr;
  TH1D *hAllEnergy7p2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hAllEnergy7p") : nullptr;
  TH1D *hAllEnergy6m2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hAllEnergy6m") : nullptr;
  TH1D *hAllEnergy7m2 = (hasSecondInput == true) ? (TH1D *)InputFile2->Get("hAllEnergy7m") : nullptr;

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
    NormalizeHistogram(hPlusLeading);
    NormalizeHistogram(hMinusLeading);
    NormalizeHistogram(hPlusLeading2);
    NormalizeHistogram(hMinusLeading2);

    if (doCumulative == true) {
      hPlusLeading = (TH1D *)hPlusLeading->GetCumulative(kFALSE);
      hMinusLeading = (TH1D *)hMinusLeading->GetCumulative(kFALSE);
      if (hPlusLeading2 != nullptr)
        hPlusLeading2 = (TH1D *)hPlusLeading2->GetCumulative(kFALSE);
      if (hMinusLeading2 != nullptr)
        hMinusLeading2 = (TH1D *)hMinusLeading2->GetCumulative(kFALSE);
    }

    DrawLeadingOverlap(hPlusLeading, hPlusLeading2, Form("%s/HFEMaxPlusLeadingOverlap.pdf", OutputFileName.c_str()));
    DrawLeadingOverlap(hMinusLeading, hMinusLeading2, Form("%s/HFEMaxMinusLeadingOverlap.pdf", OutputFileName.c_str()));
    if (hasEtaPhiHistograms == false)
      return 0;
  }

  if (hHFnPF != nullptr) {
    DrawHFnPFOverlap(hHFnPF, hHFnPF2, Form("%s/HFnPFOverlap.pdf", OutputFileName.c_str()));
  }
  vector<pair<pair<TH1D *, TH1D *>, pair<string, string>>> extra1D = {
      {{hHFnPF6p, hHFnPF6p2}, {"HFnPF6pOverlap.pdf", "HF nPF id6+"}},
      {{hHFnPF7p, hHFnPF7p2}, {"HFnPF7pOverlap.pdf", "HF nPF id7+"}},
      {{hHFnPF6m, hHFnPF6m2}, {"HFnPF6mOverlap.pdf", "HF nPF id6-"}},
      {{hHFnPF7m, hHFnPF7m2}, {"HFnPF7mOverlap.pdf", "HF nPF id7-"}},
      {{hAllEnergyp, hAllEnergyp2}, {"AllEnergypOverlap.pdf", "HF PF candidate energy + (GeV)"}},
      {{hAllEnergym, hAllEnergym2}, {"AllEnergymOverlap.pdf", "HF PF candidate energy - (GeV)"}},
      {{hAllEnergy6p, hAllEnergy6p2}, {"AllEnergy6pOverlap.pdf", "HF PF candidate energy id6+ (GeV)"}},
      {{hAllEnergy7p, hAllEnergy7p2}, {"AllEnergy7pOverlap.pdf", "HF PF candidate energy id7+ (GeV)"}},
      {{hAllEnergy6m, hAllEnergy6m2}, {"AllEnergy6mOverlap.pdf", "HF PF candidate energy id6- (GeV)"}},
      {{hAllEnergy7m, hAllEnergy7m2}, {"AllEnergy7mOverlap.pdf", "HF PF candidate energy id7- (GeV)"}}};
  for (const auto &entry : extra1D) {
    if (entry.first.first != nullptr)
      Draw1DOverlap(entry.first.first, entry.first.second,
                    Form("%s/%s", OutputFileName.c_str(), entry.second.first.c_str()), entry.second.second);
  }
  vector<pair<pair<TH1D *, TH1D *>, pair<string, string>>> sideOverlaps = {
      {{hHFnPF6p, hHFnPF6m}, {"HFnPF6SideOverlap.pdf", "HF nPF id6"}},
      {{hHFnPF7p, hHFnPF7m}, {"HFnPF7SideOverlap.pdf", "HF nPF id7"}},
      {{hAllEnergyp, hAllEnergym}, {"AllEnergySideOverlap.pdf", "HF PF candidate energy (GeV)"}},
      {{hAllEnergy6p, hAllEnergy6m}, {"AllEnergy6SideOverlap.pdf", "HF PF candidate energy id6 (GeV)"}},
      {{hAllEnergy7p, hAllEnergy7m}, {"AllEnergy7SideOverlap.pdf", "HF PF candidate energy id7 (GeV)"}}};
  for (const auto &entry : sideOverlaps) {
    if (entry.first.first != nullptr && entry.first.second != nullptr)
      Draw1DOverlap(entry.first.first, entry.first.second,
                    Form("%s/%s", OutputFileName.c_str(), entry.second.first.c_str()), entry.second.second);
  }

  ////////////////////////////////////////////////////////////
  // load histograms
  ////////////////////////////////////////////////////////////

  vector<vector<TH1D *>> hPlus(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hMinus(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hPlusCum(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hMinusCum(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hPlus2(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hMinus2(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hPlusCum2(nEta, vector<TH1D *>(nPhi, nullptr));
  vector<vector<TH1D *>> hMinusCum2(nEta, vector<TH1D *>(nPhi, nullptr));

  for (int iEta = 0; iEta < nEta; iEta++)
    for (int iPhi = 0; iPhi < nPhi; iPhi++) {
      const TString plusName = Form("hHFEMaxPlus_eta%.1f_%.1f_phi%.1f_%.1f", etaBorders[iEta], etaBorders[iEta + 1],
                                    phiBorders[iPhi], phiBorders[iPhi + 1]);
      const TString minusName = Form("hHFEMaxMinus_eta%.1f_%.1f_phi%.1f_%.1f", -etaBorders[iEta + 1], -etaBorders[iEta],
                                     phiBorders[iPhi], phiBorders[iPhi + 1]);

      hPlus[iEta][iPhi] = (TH1D *)InputFile->Get(plusName);

      if (hPlus[iEta][iPhi]->Integral() > 0)
        hPlus[iEta][iPhi]->Scale(1.0 / hPlus[iEta][iPhi]->Integral());
      if (doCumulative)
        hPlusCum[iEta][iPhi] = (TH1D *)hPlus[iEta][iPhi]->GetCumulative(kFALSE);

      hMinus[iEta][iPhi] = (TH1D *)InputFile->Get(minusName);
      if (hMinus[iEta][iPhi]->Integral() > 0)
        hMinus[iEta][iPhi]->Scale(1.0 / hMinus[iEta][iPhi]->Integral());
      if (doCumulative)
        hMinusCum[iEta][iPhi] = (TH1D *)hMinus[iEta][iPhi]->GetCumulative(kFALSE);

      if (hasSecondInput) {
        hPlus2[iEta][iPhi] = (TH1D *)InputFile2->Get(plusName);
        if (hPlus2[iEta][iPhi] != nullptr && hPlus2[iEta][iPhi]->Integral() > 0)
          hPlus2[iEta][iPhi]->Scale(1.0 / hPlus2[iEta][iPhi]->Integral());
        if (doCumulative && hPlus2[iEta][iPhi] != nullptr)
          hPlusCum2[iEta][iPhi] = (TH1D *)hPlus2[iEta][iPhi]->GetCumulative(kFALSE);

        hMinus2[iEta][iPhi] = (TH1D *)InputFile2->Get(minusName);
        if (hMinus2[iEta][iPhi] != nullptr && hMinus2[iEta][iPhi]->Integral() > 0)
          hMinus2[iEta][iPhi]->Scale(1.0 / hMinus2[iEta][iPhi]->Integral());
        if (doCumulative && hMinus2[iEta][iPhi] != nullptr)
          hMinusCum2[iEta][iPhi] = (TH1D *)hMinus2[iEta][iPhi]->GetCumulative(kFALSE);
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
  double topMargin = 0.02;

  int canvasW = padW * nEta / (1 - leftMargin - rightMargin);
  int canvasH = padH * nPhi / (1 - bottomMargin - topMargin);

  TCanvas *cPlus = new TCanvas("cPlus", "", canvasW, canvasH);
  TCanvas *cMinus = new TCanvas("cMinus", "", canvasW, canvasH);

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

      TH1D *histToDrawPlus = doCumulative ? hPlusCum[iEta][iPhi] : hPlus[iEta][iPhi];
      DrawHistogram(histToDrawPlus, showX, showY, kBlue, "hist");
      if (hasSecondInput) {
        TH1D *histToDrawPlus2 = doCumulative ? hPlusCum2[iEta][iPhi] : hPlus2[iEta][iPhi];
        DrawHistogram(histToDrawPlus2, showX, showY, kRed, "hist same");
      }
      label->DrawLatex(0.15, 0.85,
                       Form("#eta: [%.1f,%.1f], #phi: [%.1f,%.1f]", etaBorders[iEta], etaBorders[iEta + 1],
                            phiBorders[iPhi], phiBorders[iPhi + 1]));
      const double plusMean = calculateMean(hPlus[iEta][iPhi]);
      hEtaPhiMeanPlus->SetBinContent(iEta + 1, iPhi + 1, plusMean);
      hEtaPhiStdPlus->SetBinContent(iEta + 1, iPhi + 1, calculateStdDev(hPlus[iEta][iPhi], plusMean));
      hEtaPhiNumPlus->SetBinContent(iEta + 1, iPhi + 1, calculateIntegral(hPlus[iEta][iPhi]));

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

      TH1D *histToDrawMinus = doCumulative ? hMinusCum[iEta][iPhi] : hMinus[iEta][iPhi];
      DrawHistogram(histToDrawMinus, showX, showY, kBlue, "hist");
      if (hasSecondInput) {
        TH1D *histToDrawMinus2 = doCumulative ? hMinusCum2[iEta][iPhi] : hMinus2[iEta][iPhi];
        DrawHistogram(histToDrawMinus2, showX, showY, kRed, "hist same");
      }
      label->DrawLatex(0.15, 0.85,
                       Form("#eta: [%.1f,%.1f], #phi: [%.1f,%.1f]", -etaBorders[iEta + 1], -etaBorders[iEta],
                            phiBorders[iPhi], phiBorders[iPhi + 1]));
      const double minusMean = calculateMean(hMinus[iEta][iPhi]);
      hEtaPhiMeanMinus->SetBinContent(iEta + 1, iPhi + 1, minusMean);
      hEtaPhiStdMinus->SetBinContent(iEta + 1, iPhi + 1, calculateStdDev(hMinus[iEta][iPhi], minusMean));
      hEtaPhiNumMinus->SetBinContent(iEta + 1, iPhi + 1, calculateIntegral(hMinus[iEta][iPhi]));
    }
  TCanvas *cMeanPlus = new TCanvas("cMeanPlus", "", 1000, 800);
  cMeanPlus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiMeanPlus);
  hEtaPhiMeanPlus->Draw("colz");
  cMeanPlus->SaveAs(Form("%s/HFEMaxPlusMeanMap.pdf", OutputFileName.c_str()));

  TCanvas *cStdPlus = new TCanvas("cStdPlus", "", 1000, 800);
  cStdPlus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiStdPlus);
  hEtaPhiStdPlus->Draw("colz");
  cStdPlus->SaveAs(Form("%s/HFEMaxPlusStdMap.pdf", OutputFileName.c_str()));

  TCanvas *cNumPlus = new TCanvas("cNumPlus", "", 1000, 800);
  cNumPlus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiNumPlus);
  hEtaPhiNumPlus->Draw("colz");
  cNumPlus->SaveAs(Form("%s/HFEMaxPlusNumEntriesMap.pdf", OutputFileName.c_str()));

  TCanvas *cMeanMinus = new TCanvas("cMeanMinus", "", 1000, 800);
  cMeanMinus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiMeanMinus);
  hEtaPhiMeanMinus->Draw("colz");
  cMeanMinus->SaveAs(Form("%s/HFEMaxMinusMeanMap.pdf", OutputFileName.c_str()));

  TCanvas *cStdMinus = new TCanvas("cStdMinus", "", 1000, 800);
  cStdMinus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiStdMinus);
  hEtaPhiStdMinus->Draw("colz");
  cStdMinus->SaveAs(Form("%s/HFEMaxMinusStdMap.pdf", OutputFileName.c_str()));

  TCanvas *cNumMinus = new TCanvas("cNumMinus", "", 1000, 800);
  cNumMinus->cd();
  SetEtaPhiMapPadMargins(gPad);
  StyleEtaPhiMap(hEtaPhiNumMinus);
  hEtaPhiNumMinus->Draw("colz");
  cNumMinus->SaveAs(Form("%s/HFEMaxMinusNumEntriesMap.pdf", OutputFileName.c_str()));

  ////////////////////////////////////////////////////////////
  // save
  ////////////////////////////////////////////////////////////

  cPlus->SaveAs(Form("%s/HFEMaxPlusMaps.pdf", OutputFileName.c_str()));

  cMinus->SaveAs(Form("%s/HFEMaxMinusMaps.pdf", OutputFileName.c_str()));

  cPlus->SaveAs(Form("%s/HFEMaxPlusMaps.png", OutputFileName.c_str()));

  cMinus->SaveAs(Form("%s/HFEMaxMinusMaps.png", OutputFileName.c_str()));

  return 0;
}
