#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMatrixD.h"

#include <iostream>
#include <cmath>
#include <limits>
#include <string>

#include "CommandLine.h"
#include "InfoManager.h"
#include "RootIOUtils.h"

static void applyUnfoldingErrors(TH1D *hist, const TMatrixD &covMatrix) {
  if (!hist)
    return;

  const int nBins = hist->GetNbinsX();
  if (covMatrix.GetNrows() != nBins || covMatrix.GetNcols() != nBins) {
    std::cerr << "Warning: covariance matrix size (" << covMatrix.GetNrows()
              << "x" << covMatrix.GetNcols() << ") does not match histogram bins ("
              << nBins << ")" << std::endl;
    return;
  }

  for (int i = 0; i < nBins; ++i) {
    const double variance = covMatrix(i, i);
    const double error = (variance > 0.0) ? std::sqrt(variance) : 0.0;
    hist->SetBinError(i + 1, error);
  }
}

static void propagateUnfoldingErrorsToRefolded(TH1D *refoldedHist,
                                               const TH1D *unfoldedHist,
                                               const TMatrixD &unfoldCovMatrix,
                                               const RooUnfoldResponse &response) {
  if (!refoldedHist || !unfoldedHist)
    return;

  const int nTruthBins = unfoldedHist->GetNbinsX();
  const int nMeasuredBins = refoldedHist->GetNbinsX();

  const TH2 *responseH2 = response.Hresponse();
  if (!responseH2) {
    std::cerr << "Warning: cannot access response matrix for refolded error propagation"
              << std::endl;
    return;
  }

  TMatrixD responseMatrix(nMeasuredBins, nTruthBins);
  for (int i = 0; i < nMeasuredBins; ++i)
    for (int j = 0; j < nTruthBins; ++j)
      responseMatrix(i, j) = responseH2->GetBinContent(i + 1, j + 1);

  TMatrixD responseMatrixT(TMatrixD::kTransposed, responseMatrix);
  TMatrixD covRefolded = responseMatrix * (unfoldCovMatrix * responseMatrixT);
  for (int i = 0; i < nMeasuredBins; ++i) {
    const double variance = covRefolded(i, i);
    const double error = (variance > 0.0) ? std::sqrt(variance) : 0.0;
    refoldedHist->SetBinError(i + 1, error);
  }
}

static void saveCovarianceMatrix(const TMatrixD &covMatrix, const std::string &name) {
  covMatrix.Write(name.c_str());
}

static double computeChi2(const TH1 *measured, const TH1 *refolded, int &ndf) {
  ndf = 0;
  double chi2 = 0.0;

  if (!measured || !refolded || measured->GetNbinsX() != refolded->GetNbinsX())
    return 0.0;

  for (int i = 1; i <= measured->GetNbinsX(); ++i) {
    const double m = measured->GetBinContent(i);
    const double r = refolded->GetBinContent(i);
    const double em = measured->GetBinError(i);
    const double er = refolded->GetBinError(i);
    const double variance = em * em + er * er;
    if (variance <= 0)
      continue;
    const double diff = m - r;
    chi2 += diff * diff / variance;
    ndf++;
  }

  return chi2;
}

static void covarianceSummary(const TMatrixD &cov, double &trace, double &diagMin, double &diagMax) {
  trace = 0.0;
  diagMin = std::numeric_limits<double>::infinity();
  diagMax = -std::numeric_limits<double>::infinity();

  const int nDiag = (cov.GetNrows() < cov.GetNcols()) ? cov.GetNrows() : cov.GetNcols();
  for (int i = 0; i < nDiag; ++i) {
    const double v = cov(i, i);
    trace += v;
    if (v < diagMin)
      diagMin = v;
    if (v > diagMax)
      diagMax = v;
  }

  if (nDiag == 0) {
    diagMin = 0.0;
    diagMax = 0.0;
  }
}


int main(int argc, char** argv)
{
    CommandLine CL(argc, argv);

    std::string noiseFileName = CL.Get("NoiseFile", "noise.root");
    std::string DataFileName = CL.Get("DataFile", "data.root");
    std::string outputFileName = CL.Get("OutputFileName", "fit_result.root");
    std::string varNoiseName = CL.Get("VarNoiseName", "HFEMaxPlus_forest");
    std::string varDataName = CL.Get("VarDataName", "HFEMaxPlus_forest");
    int iterations            = CL.GetInt("Iterations", 4);
    int dataQuarter           = CL.GetInt("DataQuarter", 0);
    int kTermMin              = CL.GetInt("KTermMin", 1);
    int kTermMaxInput         = CL.GetInt("KTermMax", -1);
    double xMin                = CL.GetDouble("XMin", 0.0);
    double xMax                = CL.GetDouble("XMax", 100.0);
    int binsPerGeV            = CL.GetInt("BinsPerGeV", 2);

    if(xMax <= xMin)
    {
      std::cerr << "Invalid range: XMax must be greater than XMin." << std::endl;
      return -1;
    }
    if(binsPerGeV <= 0)
    {
      std::cerr << "Invalid binning: BinsPerGeV must be > 0." << std::endl;
      return -1;
    }
    if(dataQuarter < -1 || dataQuarter > 3)
    {
      std::cerr << "Invalid DataQuarter: must be in [-1,3]." << std::endl;
      return -1;
    }

    const int requestedNBins = static_cast<int>(std::lround((xMax - xMin) * binsPerGeV));
    if(requestedNBins < 1)
    {
      std::cerr << "Invalid binning: computed number of bins is < 1." << std::endl;
      return -1;
    }

    // The noise histogram is used to build the response matrix, so we load it with the full dataset (dataQuarter = -1) to get the best possible statistics for the response. The signal+noise histogram is loaded with the specified quarter to simulate a realistic measurement scenario.
    TH1D* noiseHist = RootIOUtils::LoadTreeBranchHistogramOrNull(
      noiseFileName, "OutputTree", varNoiseName, "hNoise", "hNoise;" + varNoiseName + ";Events",
      requestedNBins, xMin, xMax, -1, true, true, "noise file");
    TH1D* signalPlusNoiseHist = RootIOUtils::LoadTreeBranchHistogramOrNull(
      DataFileName, "OutputTree", varDataName, "hData", "hData;" + varDataName + ";Events",
      requestedNBins, xMin, xMax, dataQuarter, true, true, "data file");

    if(!noiseHist || !signalPlusNoiseHist)
    {
        std::cerr << "Error retrieving histograms from file." << std::endl;
        return -1;
    }

    // Normalize noise distribution (important!)
    if(noiseHist->Integral() > 0)
        noiseHist->Scale(1.0 / noiseHist->Integral());

    int nBins = signalPlusNoiseHist->GetNbinsX();
    double xmin = signalPlusNoiseHist->GetXaxis()->GetXmin();
    double xmax = signalPlusNoiseHist->GetXaxis()->GetXmax();

    // Create response matrix
    RooUnfoldResponse response(nBins, xmin, xmax);

    const double noiseNorm = noiseHist->Integral();
    if(noiseNorm <= 0)
    {
      std::cerr << "Noise histogram has zero integral after loading/normalization." << std::endl;
      return -1;
    }

    for(int j = 1; j <= nBins; j++)
    {
      const double xTrue = signalPlusNoiseHist->GetXaxis()->GetBinCenter(j);

      for(int b = 1; b <= noiseHist->GetNbinsX(); b++)
        {
        const double binContent = noiseHist->GetBinContent(b);
        if(binContent <= 0)
          continue;
          
        const double noise = noiseHist->GetXaxis()->GetBinCenter(b);
        const double xMeasured = xTrue + noise;

        if(xMeasured < xmin || xMeasured >= xmax)
          response.Miss(xTrue, binContent);
        else
          response.Fill(xMeasured, xTrue, binContent);
        }
    }

    if(iterations < 1)
    {
      std::cerr << "Iterations must be >= 1" << std::endl;
      return -1;
    }

    int kTermMax = (kTermMaxInput > 0) ? kTermMaxInput : nBins;
    if(kTermMin < 1)
    {
      std::cerr << "KTermMin must be >= 1 for SVD regularization." << std::endl;
      return -1;
    }
    if(kTermMax < kTermMin)
    {
      std::cerr << "KTermMax must be >= KTermMin." << std::endl;
      return -1;
    }
    if(kTermMax > nBins)
      kTermMax = nBins;

    TFile *outputFile = RootIOUtils::OpenFileOrNull(outputFileName, "RECREATE", "output file");
    if(outputFile == nullptr)
    {
      return -1;
    }
    TTimeStamp *currentTime = new TTimeStamp();
    GeneralInfoManager man(outputFile, "InfoDir", false);
    man.AddSourceFile(noiseFileName, currentTime);
    man.AddSourceFile(DataFileName, currentTime);
    man.AddCutParameter("Iterations", iterations, currentTime);
    man.AddCutParameter("DataQuarter", dataQuarter, currentTime);
    man.AddCutParameter("KTermMin", kTermMin, currentTime);
    man.AddCutParameter("KTermMax", kTermMaxInput, currentTime);
    man.AddCutParameter("XMin", xMin, currentTime);
    man.AddCutParameter("XMax", xMax, currentTime);
    man.AddCutParameter("BinsPerGeV", binsPerGeV, currentTime);

    signalPlusNoiseHist->SetName("hMeasured");
    noiseHist->SetName("hNoisePDF");
    signalPlusNoiseHist->Write();
    noiseHist->Write();

    TTree *regTree = new TTree("RegularizationTree", "SVD regularization parameter scan");
    int regIteration = 0;
    int regKTerm = 0;
    int regNDF = 0;
    double regChi2 = 0.0;
    double regChi2NDF = 0.0;
    double regCovTrace = 0.0;
    double regCovDiagMin = 0.0;
    double regCovDiagMax = 0.0;
    regTree->Branch("Iteration", &regIteration, "Iteration/I");
    regTree->Branch("KTerm", &regKTerm, "KTerm/I");
    regTree->Branch("NDF", &regNDF, "NDF/I");
    regTree->Branch("Chi2", &regChi2, "Chi2/D");
    regTree->Branch("Chi2NDF", &regChi2NDF, "Chi2NDF/D");
    regTree->Branch("CovTrace", &regCovTrace, "CovTrace/D");
    regTree->Branch("CovDiagMin", &regCovDiagMin, "CovDiagMin/D");
    regTree->Branch("CovDiagMax", &regCovDiagMax, "CovDiagMax/D");

    TH2 *responseMatrix = dynamic_cast<TH2 *>(response.Hresponse()->Clone("hResponseMatrix"));
    if(responseMatrix)
    {
      responseMatrix->SetDirectory(outputFile);
      responseMatrix->Write();
    }

    // Run SVD unfolding for a k-term scan and save both unfolded and refolded spectra.
    TH1D *lastUnfoldedSignal = nullptr;
    TH1 *lastRefoldedSignal = nullptr;
    for(int iter = 1; iter <= iterations; iter++)
    {
      const double frac = (iterations > 1) ? static_cast<double>(iter - 1) / (iterations - 1) : 0.0;
      int kTerm = kTermMin + static_cast<int>(std::lround(frac * (kTermMax - kTermMin)));
      if(kTerm < 1)
        kTerm = 1;
      if(kTerm > nBins)
        kTerm = nBins;

      RooUnfoldSvd unfold(&response, signalPlusNoiseHist, kTerm);
      TH1D *unfoldedTmp = dynamic_cast<TH1D *>(unfold.Hreco());
      if(!unfoldedTmp)
      {
        std::cerr << "Failed to build unfolded histogram for iteration " << iter << std::endl;
        continue;
      }

      const std::string unfoldedName = "hUnfolded_iter" + std::to_string(iter);
      TH1D *unfoldedSignal = dynamic_cast<TH1D *>(unfoldedTmp->Clone(unfoldedName.c_str()));
      if(!unfoldedSignal)
      {
        std::cerr << "Failed to clone unfolded histogram for iteration " << iter << std::endl;
        continue;
      }

      TMatrixD unfoldCovMatrix = unfold.Ereco();
      applyUnfoldingErrors(unfoldedSignal, unfoldCovMatrix);

      unfoldedSignal->SetDirectory(outputFile);
      unfoldedSignal->Write();

      const std::string covName = "cov_unfolded_iter" + std::to_string(iter);
      outputFile->cd();
      saveCovarianceMatrix(unfoldCovMatrix, covName);

      const std::string refoldedName = "hRefolded_iter" + std::to_string(iter);
      TH1 *refoldedSignal = response.ApplyToTruth(unfoldedSignal, refoldedName.c_str());
      if(!refoldedSignal)
      {
        std::cerr << "Failed to build refolded histogram for iteration " << iter << std::endl;
        continue;
      }

      TH1D *refoldedFor1D = dynamic_cast<TH1D *>(refoldedSignal);
      if(refoldedFor1D)
        propagateUnfoldingErrorsToRefolded(refoldedFor1D, unfoldedSignal, unfoldCovMatrix, response);

      refoldedSignal->SetDirectory(outputFile);
      refoldedSignal->Write();

      covarianceSummary(unfoldCovMatrix, regCovTrace, regCovDiagMin, regCovDiagMax);
      regChi2 = computeChi2(signalPlusNoiseHist, refoldedSignal, regNDF);
      regChi2NDF = (regNDF > 0) ? regChi2 / regNDF : 0.0;
      regIteration = iter;
      regKTerm = kTerm;
      regTree->Fill();

      std::cout << "Iteration " << iter
                << " uses k-term = " << kTerm
                << ", diagnostic chi2/NDF = " << regChi2NDF
                << " (" << regChi2 << "/" << regNDF << ")"
                << ", cov-trace = " << regCovTrace
                << std::endl;

      lastUnfoldedSignal = unfoldedSignal;
      lastRefoldedSignal = refoldedSignal;
    }

    regTree->Write();

    if(lastUnfoldedSignal && lastRefoldedSignal)
    {
      // Store one quick-check canvas comparing measured, unfolded and refolded for the highest successful iteration.
      TCanvas *c1 = new TCanvas("c1","Unfolding Result",800,800);
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

      TLegend *legend = new TLegend(0.55,0.68,0.88,0.9);
      legend->AddEntry(signalPlusNoiseHist,"Measured","lep");
      legend->AddEntry(lastUnfoldedSignal,"Unfolded (max iter)","lep");
      legend->AddEntry(lastRefoldedSignal,"Refolded (max iter)","lep");
      legend->Draw();

      outputFile->cd();
      c1->Write();
    }

    man.SaveToFile();
    RootIOUtils::CloseAndDeleteFile(outputFile);
    delete currentTime;

    return 0;
}
