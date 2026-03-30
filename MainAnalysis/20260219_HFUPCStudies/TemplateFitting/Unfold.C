#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <iostream>
#include <cmath>
#include <string>

#include "CommandLine.h"
#include "InfoManager.h"
#include "RootIOUtils.h"

// ====================================================================
//  Helper: Extract errors from RooUnfold covariance matrix
//  Covariance matrix diagonal elements are variances; errors = sqrt(diag)
// ====================================================================
static void applyUnfoldingErrors(TH1D *hist, const TMatrixD &covMatrix) {
  if (!hist) return;
  int nBins = hist->GetNbinsX();
  
  // Verify covariance matrix size matches histogram
  if (covMatrix.GetNrows() != nBins || covMatrix.GetNcols() != nBins) {
    std::cerr << "Warning: covariance matrix size (" << covMatrix.GetNrows() 
              << "x" << covMatrix.GetNcols() << ") does not match histogram bins (" 
              << nBins << ")" << std::endl;
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
static void propagateUnfoldingErrorsToRefolded(TH1D *refoldedHist, 
                                               const TH1D *unfoldedHist,
                                               const TMatrixD &unfoldCovMatrix,
                                               const RooUnfoldResponse &response) {
  if (!refoldedHist || !unfoldedHist) return;
  
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
static void saveCovarianceMatrix(const TMatrixD &covMatrix, const std::string &name) {
  covMatrix.Write(name.c_str());
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
    man.AddCutParameter("XMin", xMin, currentTime);
    man.AddCutParameter("XMax", xMax, currentTime);
    man.AddCutParameter("BinsPerGeV", binsPerGeV, currentTime);

    signalPlusNoiseHist->SetName("hMeasured");
    noiseHist->SetName("hNoisePDF");
    signalPlusNoiseHist->Write();
    noiseHist->Write();

    TH2 *responseMatrix = dynamic_cast<TH2 *>(response.Hresponse()->Clone("hResponseMatrix"));
    if(responseMatrix)
    {
      responseMatrix->SetDirectory(outputFile);
      responseMatrix->Write();
    }

    // Run unfolding for all iteration counts from 1..Iterations and save both unfolded and refolded spectra
    // with native RooUnfold error propagation.
    TH1D *lastUnfoldedSignal = nullptr;
    TH1 *lastRefoldedSignal = nullptr;
    for(int iter = 1; iter <= iterations; iter++)
    {
      RooUnfoldBayes unfold(&response, signalPlusNoiseHist, iter);
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

      // Extract covariance matrix from RooUnfold (public API) and apply errors.
      // Ereco() returns the unfolded covariance matrix in RooUnfold.
      TMatrixD unfoldCovMatrix = unfold.Ereco();
      applyUnfoldingErrors(unfoldedSignal, unfoldCovMatrix);
      
      printf("  Iteration %d: unfolded errors set from covariance matrix (%dx%d)\n",
             iter, unfoldCovMatrix.GetNrows(), unfoldCovMatrix.GetNcols());

      unfoldedSignal->SetDirectory(outputFile);
      unfoldedSignal->Write();

      // Save covariance matrix for later analysis if needed
      const std::string covName = "cov_unfolded_iter" + std::to_string(iter);
      outputFile->cd();
      saveCovarianceMatrix(unfoldCovMatrix, covName);

      // Apply response matrix to get refolded histogram and propagate errors
      const std::string refoldedName = "hRefolded_iter" + std::to_string(iter);
      TH1 *refoldedSignal = response.ApplyToTruth(unfoldedSignal, refoldedName.c_str());
      if(!refoldedSignal)
      {
        std::cerr << "Failed to build refolded histogram for iteration " << iter << std::endl;
        continue;
      }

      // Propagate errors from unfolded to refolded using response matrix
      TH1D *refoldedFor1D = dynamic_cast<TH1D *>(refoldedSignal);
      if(refoldedFor1D)
      {
        propagateUnfoldingErrorsToRefolded(refoldedFor1D, unfoldedSignal, unfoldCovMatrix, response);
        printf("  Iteration %d: refolded errors propagated from unfolded via response matrix\n", iter);
      }

      refoldedSignal->SetDirectory(outputFile);
      refoldedSignal->Write();

      lastUnfoldedSignal = unfoldedSignal;
      lastRefoldedSignal = refoldedSignal;
    }

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
