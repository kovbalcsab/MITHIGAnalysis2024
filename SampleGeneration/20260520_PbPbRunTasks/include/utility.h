#ifndef UTILITY_20260522_H
#define UTILITY_20260522_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

#include "CommandLine.h"
#include "hfEdges_20260520.h"

const std::vector<std::string> modeLabels = {"CountCrossCheck", "AvgEnergy", "StdDevEnergy"};
const std::vector<std::string> idLabels = {"All", "id6", "id7"};

std::vector<std::string> SplitString(const std::string &Input, char Delimiter) {
  std::vector<std::string> Result;
  std::stringstream Stream(Input);
  std::string Item;
  while (std::getline(Stream, Item, Delimiter))
    if (Item.empty() == false)
      Result.push_back(Item);
  return Result;
}

std::vector<std::string> ExpandInputFiles(const std::string &InputSpec) {
  std::vector<std::string> Result;
  if (InputSpec.size() >= 4 && InputSpec.substr(InputSpec.size() - 4) == ".txt") {
    std::ifstream InputFile(InputSpec);
    std::string Line;
    while (std::getline(InputFile, Line))
      if (Line.empty() == false)
        Result.push_back(Line);
    return Result;
  }
  return SplitString(InputSpec, ',');
}

int determineHistLoc(double Eta, double Phi, bool do18BinMerging) {
  if (std::abs(Eta) < HFPFBinEdges::EtaEdges36Positive()[HFPFBinEdges::EtaEdges36Positive().size() - 1]) {
    return 0;
  } else if (Phi < HFPFBinEdges::PhiEdges18Main()[0]) {
    return 2;
  } else if (Phi >= HFPFBinEdges::PhiEdges18Main()[HFPFBinEdges::PhiEdges18Main().size() - 1]) {
    if (do18BinMerging)
      return 2;
    return 3;
  } else {
    return 1;
  }
}

std::vector<TH2D *> makeEnergyDistributionHistograms(std::string Name_prefix, bool isPositive) {
  std::vector<TH2D *> result;
  std::vector<double> etaEdges36 = isPositive ? HFPFBinEdges::EtaEdges36Positive() : HFPFBinEdges::EtaEdges36Negative();
  std::vector<double> etaEdges18 = isPositive ? HFPFBinEdges::EtaEdges18Positive() : HFPFBinEdges::EtaEdges18Negative();
  int n36phi = HFPFBinEdges::PhiEdges36().size() - 1;
  int n18phiMain = HFPFBinEdges::PhiEdges18Main().size() - 1;

  for (int i = 0; i < (etaEdges36.size() - 1) * n36phi; i++) {
    result.push_back(
        new TH2D((Name_prefix + "_36_" + std::to_string(i)).c_str(),
                 (Name_prefix + " 36-phi ring " + std::to_string(i) + "eta: " + std::to_string(etaEdges36[i / n36phi]) +
                  "phi: " + std::to_string(HFPFBinEdges::PhiEdges36()[i % n36phi]))
                     .c_str(),
                 400, 0, 100, 2, 6, 8));
  }
  for (int i = 0; i < (etaEdges18.size() - 1) * n18phiMain; i++) {
    result.push_back(new TH2D((Name_prefix + "_18_" + std::to_string(i)).c_str(),
                              (Name_prefix + " 18-phi ring " + std::to_string(i) +
                               "eta: " + std::to_string(etaEdges18[i / n18phiMain]) +
                               "phi: " + std::to_string(HFPFBinEdges::PhiEdges18Main()[i % n18phiMain]))
                                  .c_str(),
                              400, 0, 100, 2, 6, 8));
  }
  for (int i = 0; i < (etaEdges18.size() - 1); i++) {
    result.push_back(
        new TH2D((Name_prefix + "_18_wrap_low_" + std::to_string(i)).c_str(),
                 (Name_prefix + " 18-phi wrap low ring " + std::to_string(i) + "eta: " + std::to_string(etaEdges18[i]) +
                  "phi: <" + std::to_string(HFPFBinEdges::PhiEdges18WrapLow()[0]))
                     .c_str(),
                 400, 0, 100, 2, 6, 8));
  }
  for (int i = 0; i < (etaEdges18.size() - 1); i++) {
    result.push_back(new TH2D(
        (Name_prefix + "_18_wrap_high_" + std::to_string(i)).c_str(),
        (Name_prefix + " 18-phi wrap high ring " + std::to_string(i) + "eta: " + std::to_string(etaEdges18[i]) +
         "phi: >" + std::to_string(HFPFBinEdges::PhiEdges18WrapHigh()[HFPFBinEdges::PhiEdges18WrapHigh().size() - 1]))
            .c_str(),
        400, 0, 100, 2, 6, 8));
  }
  return result;
}

int findEnergyHistogramIndex(double Eta, double Phi, bool do18BinMerging) {
  std::vector<double> etaEdges36 = Eta > 0 ? HFPFBinEdges::EtaEdges36Positive() : HFPFBinEdges::EtaEdges36Negative();
  if (std::abs(Eta) < HFPFBinEdges::EtaEdges36Positive()[HFPFBinEdges::EtaEdges36Positive().size() - 1]) {
    for (int i = 0; i < etaEdges36.size() - 1; i++) {
      if (Eta >= etaEdges36[i] && Eta < etaEdges36[i + 1]) {
        int n36phi = HFPFBinEdges::PhiEdges36().size() - 1;
        for (int j = 0; j < n36phi; j++) {
          if (Phi >= HFPFBinEdges::PhiEdges36()[j] && Phi < HFPFBinEdges::PhiEdges36()[j + 1]) {
            return i * n36phi + j;
          }
        }
      }
    }
  }
  std::vector<double> etaEdges18 = Eta > 0 ? HFPFBinEdges::EtaEdges18Positive() : HFPFBinEdges::EtaEdges18Negative();
  int eta18bin = -1;
  for (int i = 0; i < etaEdges18.size() - 1; i++) {
    if (Eta >= etaEdges18[i] && Eta < etaEdges18[i + 1]) {
      eta18bin = i;
    }
  }
  if (eta18bin < 0) {
    return -1;
  } // Should never reach here if the input is valid

  if (Phi < HFPFBinEdges::PhiEdges18Main()[0] ||
      (Phi >= HFPFBinEdges::PhiEdges18Main()[HFPFBinEdges::PhiEdges18Main().size() - 1] && do18BinMerging)) {
    return (etaEdges36.size() - 1) * (HFPFBinEdges::PhiEdges36().size() - 1) +
           (etaEdges18.size() - 1) * (HFPFBinEdges::PhiEdges18Main().size() - 1) + eta18bin; // Wrap low
  } else if (Phi >= HFPFBinEdges::PhiEdges18Main()[HFPFBinEdges::PhiEdges18Main().size() - 1]) {
    return (etaEdges36.size() - 1) * (HFPFBinEdges::PhiEdges36().size() - 1) +
           (etaEdges18.size() - 1) * (HFPFBinEdges::PhiEdges18Main().size() - 1) + (etaEdges18.size() - 1) +
           eta18bin; // Wrap high
  }
  for (int j = 0; j < HFPFBinEdges::PhiEdges18Main().size() - 1; j++) {
    if (Phi >= HFPFBinEdges::PhiEdges18Main()[j] && Phi < HFPFBinEdges::PhiEdges18Main()[j + 1]) {
      return (etaEdges36.size() - 1) * (HFPFBinEdges::PhiEdges36().size() - 1) +
             eta18bin * (HFPFBinEdges::PhiEdges18Main().size() - 1) + j;
    }
  }
  return -1; // Should never reach here if the input is valid
}

void performOperation(TH2D *hist, TH2D *energyhist, int mode, int id, int iBinX, int iBinY, long long nEvents) {
  if (mode == 0) {
    if (id) {
      hist->SetBinContent(iBinX, iBinY, ((float)energyhist->Integral(1, energyhist->GetNbinsX(), id, id)) / nEvents);
    } else {
      hist->SetBinContent(iBinX, iBinY, ((float)energyhist->GetEntries()) / nEvents);
    }
  } else if (mode == 1) {
    if (id) {
      TH1D *tmp_px = energyhist->ProjectionX((std::string(energyhist->GetName()) + "_px").c_str(), id, id);
      hist->SetBinContent(iBinX, iBinY, tmp_px->GetMean());
      delete tmp_px;
    } else {
      hist->SetBinContent(iBinX, iBinY, energyhist->GetMean(1));
    }
  } else if (mode == 2) {
    if (id) {
      TH1D *tmp_px = energyhist->ProjectionX((std::string(energyhist->GetName()) + "_px").c_str(), id, id);
      hist->SetBinContent(iBinX, iBinY, tmp_px->GetStdDev());
      delete tmp_px;
    } else {
      hist->SetBinContent(iBinX, iBinY, energyhist->GetStdDev(1));
    }
  }
}

// mode 0 == counts, 1== avg energy, 2==energy std
void fillAvgHistograms(std::vector<TH2D *> &energyHists, TH2D *h36, TH2D *h18Main, TH2D *h18Low, TH2D *h18High,
                       int mode = 0, int id = 0, long long nEvents = 1) {
  for (int iBinX = 1; iBinX <= h36->GetNbinsX(); iBinX++) {
    for (int iBinY = 1; iBinY <= h36->GetNbinsY(); iBinY++) {
      performOperation(h36, energyHists[(iBinX - 1) * h36->GetNbinsY() + iBinY - 1], mode, id, iBinX, iBinY, nEvents);
    }
  }
  for (int iBinX = 1; iBinX <= h18Main->GetNbinsX(); iBinX++) {
    for (int iBinY = 1; iBinY <= h18Main->GetNbinsY(); iBinY++) {
      performOperation(
          h18Main, energyHists[(h36->GetNbinsX() * h36->GetNbinsY()) + (iBinX - 1) * h18Main->GetNbinsY() + iBinY - 1],
          mode, id, iBinX, iBinY, nEvents);
    }
  }
  for (int iBinX = 1; iBinX <= h18Low->GetNbinsX(); iBinX++) {
    performOperation(h18Low,
                     energyHists[(h36->GetNbinsX() * h36->GetNbinsY()) + (h18Main->GetNbinsX() * h18Main->GetNbinsY()) +
                                 (iBinX - 1)],
                     mode, id, iBinX, 1, nEvents);
  }
  for (int iBinX = 1; iBinX <= h18High->GetNbinsX(); iBinX++) {
    performOperation(h18High,
                     energyHists[(h36->GetNbinsX() * h36->GetNbinsY()) + (h18Main->GetNbinsX() * h18Main->GetNbinsY()) +
                                 (h18Low->GetNbinsX()) + (iBinX - 1)],
                     mode, id, iBinX, 1, nEvents);
  }
}
std::vector<TH2D *> loadEnergyDistributionHistograms(TFile &input, std::string Name_prefix, bool isPositive) {
  std::vector<TH2D *> result;
  std::vector<double> etaEdges36 = isPositive ? HFPFBinEdges::EtaEdges36Positive() : HFPFBinEdges::EtaEdges36Negative();
  std::vector<double> etaEdges18 = isPositive ? HFPFBinEdges::EtaEdges18Positive() : HFPFBinEdges::EtaEdges18Negative();
  int n36phi = HFPFBinEdges::PhiEdges36().size() - 1;
  int n18phiMain = HFPFBinEdges::PhiEdges18Main().size() - 1;

  for (int i = 0; i < (etaEdges36.size() - 1) * n36phi; i++) {
    result.push_back((TH2D *)input.Get((Name_prefix + "_36_" + std::to_string(i)).c_str()));
  }
  for (int i = 0; i < (etaEdges18.size() - 1) * n18phiMain; i++) {
    result.push_back((TH2D *)input.Get((Name_prefix + "_18_" + std::to_string(i)).c_str()));
  }
  for (int i = 0; i < (etaEdges18.size() - 1); i++) {
    result.push_back((TH2D *)input.Get((Name_prefix + "_18_wrap_low_" + std::to_string(i)).c_str()));
  }
  for (int i = 0; i < (etaEdges18.size() - 1); i++) {
    result.push_back((TH2D *)input.Get((Name_prefix + "_18_wrap_high_" + std::to_string(i)).c_str()));
  }
  return result;
}

std::vector<TH2D *> loadMetaHistograms(TFile &Input) {
  std::vector<TH2D *> result;

  for (int iMode = 0; iMode < 3; iMode++) {
    for (int iID = 0; iID < 3; iID++) {
      result.push_back((TH2D *)Input.Get(("h36Negative_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h18NegativeMain_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h18NegativeWrapLow_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h18NegativeWrapHigh_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h36Positive_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h18PositiveMain_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h18PositiveWrapLow_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
      result.push_back((TH2D *)Input.Get(("h18PositiveWrapHigh_" + modeLabels[iMode] + "_" + idLabels[iID]).c_str()));
    }
  }
  return result;
}

TH2D *diffHistograms(TH2D *hist1, TH2D *hist2) {
  if (hist1->GetNbinsX() != hist2->GetNbinsX() || hist1->GetNbinsY() != hist2->GetNbinsY()) {
    std::cerr << "Error: Histograms have different binning!" << std::endl;
    return nullptr;
  }
  TH2D *result = (TH2D *)hist1->Clone((std::string(hist1->GetName()) + "_diff").c_str());
  result->SetTitle((std::string(hist1->GetTitle()) + " - " + hist2->GetTitle()).c_str());
  result->Add(hist2, -1);
  result->Divide(hist2);
  return result;
}

#endif
