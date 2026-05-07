#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CommandLine.h"
#include "InfoManager.h"
#include "RootIOUtils.h"

// ====================================================================
//  Helper: Extract errors from RooUnfold covariance matrix
//  Covariance matrix diagonal elements are variances; errors = sqrt(diag)
// ====================================================================
static void applyUnfoldingErrors(TH1D *hist, const TMatrixD &covMatrix) {
  if (!hist)
    return;
  int nBins = hist->GetNbinsX();

  // Verify covariance matrix size matches histogram
  if (covMatrix.GetNrows() != nBins || covMatrix.GetNcols() != nBins) {
    std::cerr << "Warning: covariance matrix size (" << covMatrix.GetNrows() << "x" << covMatrix.GetNcols()
              << ") does not match histogram bins (" << nBins << ")" << std::endl;
    return;
  }

  // Extract errors from covariance matrix diagonal
  for (int i = 0; i < nBins; ++i) {
    double variance = covMatrix(i, i);
    double error = (variance > 0) ? std::sqrt(variance) : 0.0;
    hist->SetBinError(i + 1, error);
  }
}

// ====================================================================
//  Helper: Propagate errors from unfolded to refolded using response matrix
//
//  Theory: If A = M·x where x is unfolded and M is the response matrix,
//  then Cov(A) = M · Cov(x) · M^T (linear error propagation)
//  Refolded variances: σ²(A_i) = [M·Cov(x)·M^T]_{ii}
// ====================================================================
static void propagateUnfoldingErrorsToRefolded(TH1D *refoldedHist, const TH1D *unfoldedHist,
                                               const TMatrixD &unfoldCovMatrix, const RooUnfoldResponse &response) {
  if (!refoldedHist || !unfoldedHist)
    return;

  int nBins = unfoldedHist->GetNbinsX();
  int nMeasBins = refoldedHist->GetNbinsX();

  // Get the response matrix as TMatrixD (mapping from truth to measured)
  // Note: response.Hresponse() is a TH2 histogram
  const TH2 *responseH2 = response.Hresponse();
  if (!responseH2) {
    std::cerr << "Warning: cannot access response matrix for error propagation" << std::endl;
    return;
  }

  // Build response matrix M from the 2D histogram
  // M[measured,truth] gives the transfer probabilities
  TMatrixD M(nMeasBins, nBins);
  for (int i = 0; i < nMeasBins; ++i) {
    for (int j = 0; j < nBins; ++j) {
      // responseH2 is indexed with TH2, so GetBinContent(i+1, j+1) for 1-based indexing
      M(i, j) = responseH2->GetBinContent(i + 1, j + 1);
    }
  }

  // Compute refolded covariance: Cov_refolded = M · Cov_unfolded · M^T
  TMatrixD Mt(TMatrixD::kTransposed, M);
  TMatrixD covRefolded = M * (unfoldCovMatrix * Mt);

  // Extract errors from refolded covariance diagonal
  for (int i = 0; i < nMeasBins; ++i) {
    double variance = covRefolded(i, i);
    double error = (variance > 0) ? std::sqrt(variance) : 0.0;
    refoldedHist->SetBinError(i + 1, error);
  }
}

// ====================================================================
//  Helper: Save covariance matrix to ROOT file as a TMatrixD
// ====================================================================
static void saveCovarianceMatrix(const TMatrixD &covMatrix, const std::string &name) { covMatrix.Write(name.c_str()); }

static bool validateBinEdges(const std::vector<double> &edges, const std::string &label) {
  if (edges.empty())
    return true;
  if (edges.size() < 2) {
    std::cerr << label << " must contain at least 2 edge values." << std::endl;
    return false;
  }
  for (size_t i = 1; i < edges.size(); ++i) {
    if (!(edges[i] > edges[i - 1])) {
      std::cerr << label << " must be strictly increasing. Bad pair at index " << (i - 1) << " and " << i << "."
                << std::endl;
      return false;
    }
  }
  return true;
}

static std::vector<double> buildUniformEdges(double xMin, double xMax, int nBins) {
  std::vector<double> edges;
  edges.reserve(nBins + 1);
  const double width = (xMax - xMin) / nBins;
  for (int i = 0; i <= nBins; ++i)
    edges.push_back(xMin + i * width);
  edges.back() = xMax;
  return edges;
}

static long long gcdLongLong(long long a, long long b) {
  if (a < 0)
    a = -a;
  if (b < 0)
    b = -b;
  while (b != 0) {
    const long long r = a % b;
    a = b;
    b = r;
  }
  return a;
}

static bool appendBinWidthUnits(const std::vector<double> &edges, long long scale, const std::string &label,
                                std::vector<long long> &units) {
  if (edges.size() < 2) {
    std::cerr << "Cannot build fine response grid from " << label << ": insufficient bin edges." << std::endl;
    return false;
  }

  for (size_t i = 1; i < edges.size(); ++i) {
    const double width = edges[i] - edges[i - 1];
    if (width <= 0.0) {
      std::cerr << "Cannot build fine response grid from " << label << ": non-positive width at bin " << (i - 1)
                << "." << std::endl;
      return false;
    }
    const long long unit = static_cast<long long>(std::llround(width * scale));
    if (unit <= 0) {
      std::cerr << "Cannot build fine response grid from " << label
                << ": width too small for numeric resolution at bin " << (i - 1) << "." << std::endl;
      return false;
    }
    units.push_back(unit);
  }
  return true;
}

static bool computeCompatibleFineBinWidth(const std::vector<double> &measuredEdges, const std::vector<double> &unfoldedEdges,
                                          double &fineBinWidth) {
  const long long scale = 1000000LL;
  std::vector<long long> units;
  if (!appendBinWidthUnits(measuredEdges, scale, "MeasuredBinEdges", units))
    return false;
  if (!appendBinWidthUnits(unfoldedEdges, scale, "UnfoldedBinEdges", units))
    return false;
  if (units.empty())
    return false;

  long long commonUnit = units[0];
  for (size_t i = 1; i < units.size(); ++i)
    commonUnit = gcdLongLong(commonUnit, units[i]);
  if (commonUnit <= 0) {
    std::cerr << "Failed to derive positive compatible fine bin width." << std::endl;
    return false;
  }

  fineBinWidth = static_cast<double>(commonUnit) / scale;
  if (fineBinWidth <= 0.0) {
    std::cerr << "Failed to derive positive compatible fine bin width." << std::endl;
    return false;
  }
  return true;
}

static int computeFineNBins(double xMin, double xMax, double fineBinWidth, const std::string &label) {
  const double span = xMax - xMin;
  if (span <= 0.0 || fineBinWidth <= 0.0)
    return -1;

  const double rawBins = span / fineBinWidth;
  const int nBins = static_cast<int>(std::lround(rawBins));
  const double tolerance = 1e-6 * std::max(1.0, std::fabs(rawBins));
  if (nBins < 1 || std::fabs(rawBins - nBins) > tolerance) {
    std::cerr << "Incompatible " << label << " range with derived fine bin width. "
              << "Range span = " << span << ", fine width = " << fineBinWidth << ", raw bins = " << rawBins
              << std::endl;
    return -1;
  }
  return nBins;
}
                               
static TH1D *loadTreeHistogramWithOptionalEdges(const std::string &fileName, const std::string &treeName,
                                                const std::string &branchName, const std::string &histName,
                                                const std::string &histTitle, int nBins, double xMin, double xMax,
                                                const std::vector<double> &customEdges, int dataQuarter,
                                                const std::string &label) {
  if (customEdges.empty()) {
    return RootIOUtils::LoadTreeBranchHistogramOrNull(fileName, treeName, branchName, histName, histTitle, nBins, xMin,
                                                      xMax, dataQuarter, true, true, label);
  }

  TFile *file = RootIOUtils::OpenFileOrNull(fileName, "READ", label);
  if (file == nullptr)
    return nullptr;

  TTree *tree = RootIOUtils::GetTreeOrNull(file, treeName, label + ": " + fileName);
  if (tree == nullptr) {
    RootIOUtils::CloseAndDeleteFile(file);
    return nullptr;
  }

  TBranch *branch = RootIOUtils::RequireBranchOrNull(tree, branchName, label + ": " + treeName);
  if (branch == nullptr) {
    RootIOUtils::CloseAndDeleteFile(file);
    return nullptr;
  }

  TLeaf *leaf = branch->GetLeaf(branchName.c_str());
  if (leaf == nullptr)
    leaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->First());
  const std::string leafType = (leaf != nullptr && leaf->GetTypeName() != nullptr) ? leaf->GetTypeName() : "";
  const bool isFloat = (leafType == "Float_t" || leafType == "float");

  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  TH1D *hist = new TH1D(histName.c_str(), histTitle.c_str(), customEdges.size() - 1, customEdges.data());
  TH1::AddDirectory(oldAddDirectory);
  hist->Sumw2();

  double valueDouble = 0.0;
  float valueFloat = 0.0f;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(branchName.c_str(), 1);
  if (isFloat)
    tree->SetBranchAddress(branchName.c_str(), &valueFloat);
  else
    tree->SetBranchAddress(branchName.c_str(), &valueDouble);

  const Long64_t nEntries = tree->GetEntries();
  Long64_t nSelected = 0;
  for (Long64_t i = 0; i < nEntries; ++i) {
    if (dataQuarter >= 0 && (i % 4) != dataQuarter)
      continue;
    tree->GetEntry(i);
    hist->Fill(isFloat ? static_cast<double>(valueFloat) : valueDouble);
    nSelected++;
  }

  tree->SetBranchStatus("*", 1);
  tree->ResetBranchAddresses();
  RootIOUtils::CloseAndDeleteFile(file);

  if (dataQuarter >= 0)
    std::cout << "  Loaded " << nSelected << "/" << nEntries << " entries from " << fileName << "/" << treeName << "/"
              << branchName << " (quarter " << dataQuarter << "/4)" << std::endl;
  else
    std::cout << "  Loaded " << nSelected << " entries from " << fileName << "/" << treeName << "/" << branchName
              << std::endl;

  return hist;
}

int main(int argc, char **argv) {
  CommandLine CL(argc, argv);

  std::string noiseFileName = CL.Get("NoiseFile", "noise.root");
  std::string DataFileName = CL.Get("DataFile", "data.root");
  std::string outputFileName = CL.Get("OutputFileName", "fit_result.root");
  std::string varNoiseName = CL.Get("VarNoiseName", "HFEMaxPlus_forest");
  std::string varDataName = CL.Get("VarDataName", "HFEMaxPlus_forest");
  int iterations = CL.GetInt("Iterations", 4);
  int dataQuarter = CL.GetInt("DataQuarter", 0);
  double xMin = CL.GetDouble("XMin", 0.0);
  double xMax = CL.GetDouble("XMax", 100.0);
  int binsPerGeV = CL.GetInt("BinsPerGeV", 2);
  std::vector<double> measuredBinEdges = CL.GetDoubleVector("MeasuredBinEdges", std::vector<double>());
  std::vector<double> unfoldedBinEdges = CL.GetDoubleVector("UnfoldedBinEdges", std::vector<double>());
  std::vector<double> legacyBinEdges = CL.GetDoubleVector("BinEdges", std::vector<double>());

  if (!legacyBinEdges.empty()) {
    if (measuredBinEdges.empty())
      measuredBinEdges = legacyBinEdges;
    if (unfoldedBinEdges.empty())
      unfoldedBinEdges = legacyBinEdges;
  }

  const bool useCustomMeasuredEdges = !measuredBinEdges.empty();
  const bool useCustomUnfoldedEdges = !unfoldedBinEdges.empty();
  if (!validateBinEdges(measuredBinEdges, "MeasuredBinEdges"))
    return -1;
  if (!validateBinEdges(unfoldedBinEdges, "UnfoldedBinEdges"))
    return -1;

  const bool needFixedBinning = (!useCustomMeasuredEdges || !useCustomUnfoldedEdges);
  if (needFixedBinning && xMax <= xMin) {
    std::cerr << "Invalid range: XMax must be greater than XMin." << std::endl;
    return -1;
  }
  if (needFixedBinning && binsPerGeV <= 0) {
    std::cerr << "Invalid binning: BinsPerGeV must be > 0." << std::endl;
    return -1;
  }
  if (dataQuarter < -1 || dataQuarter > 3) {
    std::cerr << "Invalid DataQuarter: must be in [-1,3]." << std::endl;
    return -1;
  }

  int requestedNBins = 0;
  if (needFixedBinning) {
    requestedNBins = static_cast<int>(std::lround((xMax - xMin) * binsPerGeV));
    if (requestedNBins < 1) {
      std::cerr << "Invalid binning: computed number of bins is < 1." << std::endl;
      return -1;
    }
  }

  const int measuredNBins = useCustomMeasuredEdges ? static_cast<int>(measuredBinEdges.size() - 1) : requestedNBins;
  const double measuredMin = useCustomMeasuredEdges ? measuredBinEdges.front() : xMin;
  const double measuredMax = useCustomMeasuredEdges ? measuredBinEdges.back() : xMax;
  const int unfoldedNBins = useCustomUnfoldedEdges ? static_cast<int>(unfoldedBinEdges.size() - 1) : requestedNBins;
  const double unfoldedMin = useCustomUnfoldedEdges ? unfoldedBinEdges.front() : xMin;
  const double unfoldedMax = useCustomUnfoldedEdges ? unfoldedBinEdges.back() : xMax;
  const std::vector<double> measuredEdges = useCustomMeasuredEdges ? measuredBinEdges : buildUniformEdges(measuredMin, measuredMax, measuredNBins);
  const std::vector<double> unfoldedEdges = useCustomUnfoldedEdges ? unfoldedBinEdges : buildUniformEdges(unfoldedMin, unfoldedMax, unfoldedNBins);

  double responseFineBinWidth = 0.0;
  if (!computeCompatibleFineBinWidth(measuredEdges, unfoldedEdges, responseFineBinWidth))
    return -1;

  const int fineMeasuredNBins = computeFineNBins(measuredMin, measuredMax, responseFineBinWidth, "measured");
  const int fineUnfoldedNBins = computeFineNBins(unfoldedMin, unfoldedMax, responseFineBinWidth, "unfolded");
  if (fineMeasuredNBins < 1 || fineUnfoldedNBins < 1)
    return -1;

  const int maxFineBins = 200000;
  if (fineMeasuredNBins > maxFineBins || fineUnfoldedNBins > maxFineBins) {
    std::cerr << "Derived fine binning is too granular for this macro. "
              << "Measured fine bins = " << fineMeasuredNBins << ", unfolded fine bins = " << fineUnfoldedNBins
              << "." << std::endl;
    return -1;
  }

  // The noise histogram is used to build the response matrix, so we load it with the full dataset (dataQuarter = -1) to
  // get the best possible statistics for the response. The signal+noise histogram is loaded with the specified quarter
  // to simulate a realistic measurement scenario.
  TH1D *noiseHist =
      loadTreeHistogramWithOptionalEdges(noiseFileName, "OutputTree", varNoiseName, "hNoise",
                                         "hNoise;" + varNoiseName + ";Events", measuredNBins, measuredMin, measuredMax,
                                         measuredBinEdges, -1, "noise file");
  std::unique_ptr<TH1D> noiseFineHist(loadTreeHistogramWithOptionalEdges(
      noiseFileName, "OutputTree", varNoiseName, "hNoiseFine", "hNoiseFine;" + varNoiseName + ";Events",
      fineMeasuredNBins, measuredMin, measuredMax, std::vector<double>(), -1, "noise file (fine)"));
  TH1D *signalPlusNoiseHist = loadTreeHistogramWithOptionalEdges(
      DataFileName, "OutputTree", varDataName, "hData", "hData;" + varDataName + ";Events", measuredNBins, measuredMin,
      measuredMax, measuredBinEdges, dataQuarter, "data file");

  if (!noiseHist || !noiseFineHist || !signalPlusNoiseHist) {
    std::cerr << "Error retrieving histograms from file." << std::endl;
    return -1;
  }

  // Normalize noise distribution (important!)
  if (noiseHist->Integral() > 0)
    noiseHist->Scale(1.0 / noiseHist->Integral());
  else {
    std::cerr << "Noise histogram has zero integral." << std::endl;
    return -1;
  }
  if (noiseFineHist->Integral() > 0)
    noiseFineHist->Scale(1.0 / noiseFineHist->Integral());
  else {
    std::cerr << "Fine noise histogram has zero integral." << std::endl;
    return -1;
  }

  std::unique_ptr<TH1D> truthTemplate;
  truthTemplate.reset(new TH1D("hTruthTemplate", "hTruthTemplate", unfoldedNBins, unfoldedEdges.data()));
  truthTemplate->SetDirectory(nullptr);
  std::unique_ptr<TH1D> truthFineTemplate(
      new TH1D("hTruthFineTemplate", "hTruthFineTemplate", fineUnfoldedNBins, unfoldedMin, unfoldedMax));
  truthFineTemplate->SetDirectory(nullptr);

  // Create response matrix
  RooUnfoldResponse response(signalPlusNoiseHist, truthTemplate.get());

  const TAxis *truthAxis = truthTemplate->GetXaxis();
  const TAxis *truthFineAxis = truthFineTemplate->GetXaxis();
  const TAxis *noiseFineAxis = noiseFineHist->GetXaxis();
  const TAxis *measuredAxis = signalPlusNoiseHist->GetXaxis();
  const int measuredNBinsAxis = measuredAxis->GetNbins();
  std::vector<int> truthFineCounts(truthTemplate->GetNbinsX() + 2, 0);

  for (int jFine = 1; jFine <= truthFineTemplate->GetNbinsX(); ++jFine) {
    const double xTrue = truthFineAxis->GetBinCenter(jFine);
    const int truthBin = truthAxis->FindFixBin(xTrue);
    if (truthBin >= 1 && truthBin <= truthTemplate->GetNbinsX())
      truthFineCounts[truthBin]++;
  }

  std::cout << "Building response matrix with fine tiling width " << responseFineBinWidth
            << " (measured fine bins = " << fineMeasuredNBins
            << ", unfolded fine bins = " << fineUnfoldedNBins << ")." << std::endl;
  for (int jFine = 1; jFine <= truthFineTemplate->GetNbinsX(); ++jFine) {
    const double xTrue = truthFineAxis->GetBinCenter(jFine);
    const int truthBin = truthAxis->FindFixBin(xTrue);
    if (truthBin < 1 || truthBin > truthTemplate->GetNbinsX())
      continue;

    const int nTruthFine = truthFineCounts[truthBin];
    if (nTruthFine < 1)
      continue;

    for (int bFine = 1; bFine <= noiseFineHist->GetNbinsX(); ++bFine) {
      const double noiseWeight = noiseFineHist->GetBinContent(bFine);
      if (noiseWeight <= 0)
        continue;

      const double noise = noiseFineAxis->GetBinCenter(bFine);
      const double xMeasured = xTrue + noise;
      const double weight = noiseWeight / nTruthFine;
      const int measuredBin = measuredAxis->FindFixBin(xMeasured);

      if (measuredBin < 1 || measuredBin > measuredNBinsAxis)
        response.Miss(xTrue, weight);
      else
        response.Fill(xMeasured, xTrue, weight);
    }
  }

  if (iterations < 1) {
    std::cerr << "Iterations must be >= 1" << std::endl;
    return -1;
  }

  TFile *outputFile = RootIOUtils::OpenFileOrNull(outputFileName, "RECREATE", "output file");
  if (outputFile == nullptr) {
    return -1;
  }
  TTimeStamp *currentTime = new TTimeStamp();
  GeneralInfoManager man(outputFile, "InfoDir", false);
  man.AddSourceFile(noiseFileName, currentTime);
  man.AddSourceFile(DataFileName, currentTime);
  man.AddCutParameter("Iterations", iterations, currentTime);
  man.AddCutParameter("DataQuarter", dataQuarter, currentTime);
  man.AddCutParameter("XMin", xMin, currentTime);
  man.AddCutParameter("XMax", xMax, currentTime);
  man.AddCutParameter("BinsPerGeV", binsPerGeV, currentTime);
  man.AddCutParameter("UseMeasuredBinEdges", useCustomMeasuredEdges ? 1.0 : 0.0, currentTime);
  man.AddCutParameter("UseUnfoldedBinEdges", useCustomUnfoldedEdges ? 1.0 : 0.0, currentTime);
  man.AddCutParameter("ResponseFineBinWidth", responseFineBinWidth, currentTime);
  man.AddCutParameter("ResponseFineMeasuredNBins", fineMeasuredNBins, currentTime);
  man.AddCutParameter("ResponseFineUnfoldedNBins", fineUnfoldedNBins, currentTime);

  signalPlusNoiseHist->SetName("hMeasured");
  noiseHist->SetName("hNoisePDF");
  noiseFineHist->SetName("hNoisePDFFine");
  signalPlusNoiseHist->Write();
  noiseHist->Write();
  noiseFineHist->Write();

  TH2 *responseMatrix = dynamic_cast<TH2 *>(response.Hresponse()->Clone("hResponseMatrix"));
  if (responseMatrix) {
    responseMatrix->SetDirectory(outputFile);
    responseMatrix->Write();
  }

  // Run unfolding for all iteration counts from 1..Iterations and save both unfolded and refolded spectra
  // with native RooUnfold error propagation.
  TH1D *lastUnfoldedSignal = nullptr;
  TH1 *lastRefoldedSignal = nullptr;
  for (int iter = 1; iter <= iterations; iter++) {
    RooUnfoldBayes unfold(&response, signalPlusNoiseHist, iter);
    TH1D *unfoldedTmp = dynamic_cast<TH1D *>(unfold.Hreco());
    if (!unfoldedTmp) {
      std::cerr << "Failed to build unfolded histogram for iteration " << iter << std::endl;
      continue;
    }

    const std::string unfoldedName = "hUnfolded_iter" + std::to_string(iter);
    TH1D *unfoldedSignal = dynamic_cast<TH1D *>(unfoldedTmp->Clone(unfoldedName.c_str()));
    if (!unfoldedSignal) {
      std::cerr << "Failed to clone unfolded histogram for iteration " << iter << std::endl;
      continue;
    }

    // Extract covariance matrix from RooUnfold (public API) and apply errors.
    // Ereco() returns the unfolded covariance matrix in RooUnfold.
    TMatrixD unfoldCovMatrix = unfold.Ereco();
    applyUnfoldingErrors(unfoldedSignal, unfoldCovMatrix);

    printf("  Iteration %d: unfolded errors set from covariance matrix (%dx%d)\n", iter, unfoldCovMatrix.GetNrows(),
           unfoldCovMatrix.GetNcols());

    unfoldedSignal->SetDirectory(outputFile);
    unfoldedSignal->Write();

    // Save covariance matrix for later analysis if needed
    const std::string covName = "cov_unfolded_iter" + std::to_string(iter);
    outputFile->cd();
    saveCovarianceMatrix(unfoldCovMatrix, covName);

    // Apply response matrix to get refolded histogram and propagate errors
    const std::string refoldedName = "hRefolded_iter" + std::to_string(iter);
    TH1 *refoldedSignal = response.ApplyToTruth(unfoldedSignal, refoldedName.c_str());
    if (!refoldedSignal) {
      std::cerr << "Failed to build refolded histogram for iteration " << iter << std::endl;
      continue;
    }

    // Propagate errors from unfolded to refolded using response matrix
    TH1D *refoldedFor1D = dynamic_cast<TH1D *>(refoldedSignal);
    if (refoldedFor1D) {
      propagateUnfoldingErrorsToRefolded(refoldedFor1D, unfoldedSignal, unfoldCovMatrix, response);
      printf("  Iteration %d: refolded errors propagated from unfolded via response matrix\n", iter);
    }

    refoldedSignal->SetDirectory(outputFile);
    refoldedSignal->Write();

    lastUnfoldedSignal = unfoldedSignal;
    lastRefoldedSignal = refoldedSignal;
  }

  if (lastUnfoldedSignal && lastRefoldedSignal) {
    // Store one quick-check canvas comparing measured, unfolded and refolded for the highest successful iteration.
    TCanvas *c1 = new TCanvas("c1", "Unfolding Result", 800, 800);
    c1->SetLogy();

    signalPlusNoiseHist->SetMarkerStyle(20);
    signalPlusNoiseHist->SetMarkerColor(kBlue + 1);
    signalPlusNoiseHist->SetLineColor(kBlue + 1);

    lastUnfoldedSignal->SetMarkerStyle(20);
    lastUnfoldedSignal->SetMarkerColor(kRed + 1);
    lastUnfoldedSignal->SetLineColor(kRed + 1);

    lastRefoldedSignal->SetMarkerStyle(24);
    lastRefoldedSignal->SetMarkerColor(kGreen + 2);
    lastRefoldedSignal->SetLineColor(kGreen + 2);

    signalPlusNoiseHist->Draw("E");
    lastUnfoldedSignal->Draw("E SAME");
    lastRefoldedSignal->Draw("E SAME");

    TLegend *legend = new TLegend(0.55, 0.68, 0.88, 0.9);
    legend->AddEntry(signalPlusNoiseHist, "Measured", "lep");
    legend->AddEntry(lastUnfoldedSignal, "Unfolded (max iter)", "lep");
    legend->AddEntry(lastRefoldedSignal, "Refolded (max iter)", "lep");
    legend->Draw();

    outputFile->cd();
    c1->Write();
  }

  man.SaveToFile();
  RootIOUtils::CloseAndDeleteFile(outputFile);
  delete currentTime;

  return 0;
}
