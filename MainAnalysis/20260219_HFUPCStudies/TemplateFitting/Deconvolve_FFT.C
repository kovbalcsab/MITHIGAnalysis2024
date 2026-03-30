#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TVirtualFFT.h"

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "CommandLine.h"
#include "InfoManager.h"
#include "RootIOUtils.h"

// ====================================================================
//  Wiener FFT Deconvolution with damping factor iteration
//  
//  Theory:
//    - Observed: D(x) = S(x) ⊗ K(x)  (signal convolved with kernel)
//    - FFT domain: D̂ = Ŝ · K̂
//    - Wiener filter: Ŝ = (K̂* · D̂) / (|K̂|² + α)
//    - α is the damping/regularization factor (higher = more smoothing)
//
//  Implementation details (from Deconvolve.C):
//    1. Zero-pad to 2N for linear convolution (not circular)
//    2. Kernel center kc = round(-xMin/binWidth) maps physical zero-lag to bin index
//    3. Kernel shift: kShift[i] = k[(i+kc)%N] moves zero-lag to bin 0 before FFT
//    4. Real-to-complex FFT with proper Hermitian symmetry
// ====================================================================
static TH1D *wienerDeconvolveFFT(TH1D *hData, TH1D *hKernel, double dampingFactor, const char *outName) {
  int N = hData->GetNbinsX();
  int M = 2 * N;   // Zero-pad to 2N for linear (non-circular) convolution

  // Calculate kernel center (bin index where physical zero-displacement occurs)
  // Model: obs(x) = integral signal(x') * kernel(x-x') dx'
  // Kernel's zero-lag at x=0 maps to bin round(-xMin/binWidth)
  double binWidth = hData->GetXaxis()->GetBinWidth(1);
  double xMinVal  = hData->GetXaxis()->GetXmin();
  int kc = static_cast<int>(std::round(-xMinVal / binWidth));
  kc = std::max(0, std::min(kc, N - 1));

  printf("    Wiener deconv: α=%.3e, N=%d, M=%d, kc=%d\n", dampingFactor, N, M, kc);

  // Load data and kernel into zero-padded arrays
  std::vector<double> d(M, 0), k(M, 0);
  for (int i = 0; i < N; ++i) {
    d[i] = hData->GetBinContent(i + 1);
    k[i] = hKernel->GetBinContent(i + 1);
  }

  // Shift kernel so zero-lag bin (kc) lands at index 0 before FFT
  // This is the correct ifftshift for kernel with zero-displacement at bin kc
  std::vector<double> kShift(M, 0);
  for (int i = 0; i < N; ++i)
    kShift[i] = k[(i + kc) % N];

  // FFT helper: real-to-complex forward transform with Hermitian symmetry
  auto fwd = [&](const std::vector<double> &v, std::vector<double> &re, std::vector<double> &im) {
    int n = v.size();
    re.assign(n, 0);
    im.assign(n, 0);
    TVirtualFFT *f = TVirtualFFT::FFT(1, &n, "R2C EX K");
    f->SetPoints(const_cast<double *>(v.data()));
    f->Transform();
    for (int k2 = 0; k2 < n / 2 + 1; ++k2) {
      Double_t r, ii;
      f->GetPointComplex(k2, r, ii);
      re[k2] = r;
      im[k2] = ii;
      // Mirror for Hermitian symmetry: conj(F[N-k]) = F[k]
      if (k2 > 0 && k2 < n - k2) {
        re[n - k2] = r;
        im[n - k2] = -ii;
      }
    }
    delete f;
  };

  // Transform both data and shifted kernel to frequency domain
  std::vector<double> dRe, dIm, kRe, kIm;
  fwd(d, dRe, dIm);
  fwd(kShift, kRe, kIm);

  // Wiener filter in frequency domain
  // Ŝ = (K̂* · D̂) / (|K̂|² + α)
  std::vector<double> sRe(M, 0), sIm(M, 0);
  for (int k2 = 0; k2 < M / 2 + 1; ++k2) {
    double kp2 = kRe[k2] * kRe[k2] + kIm[k2] * kIm[k2];
    
    // Avoid division by zero: if kernel is too small, zero out this frequency component
    if (kp2 < 1e-30) {
      sRe[k2] = sIm[k2] = 0;
      continue;
    }

    // Numerator: K̂* · D̂ = (kRe + i·kIm)* · (dRe + i·dIm)
    //                    = (kRe - i·kIm) · (dRe + i·dIm)
    //                    = kRe·dRe + kIm·dIm + i(kRe·dIm - kIm·dRe)
    double numRe = kRe[k2] * dRe[k2] + kIm[k2] * dIm[k2];
    double numIm = kRe[k2] * dIm[k2] - kIm[k2] * dRe[k2];

    // Denominator: |K̂|² + α
    double denom = kp2 + dampingFactor;

    sRe[k2] = numRe / denom;
    sIm[k2] = numIm / denom;

    // Maintain Hermitian symmetry for inverse FFT
    if (k2 > 0 && k2 < M - k2) {
      sRe[M - k2] = sRe[k2];
      sIm[M - k2] = -sIm[k2];
    }
  }

  // Inverse FFT: complex-to-real
  TVirtualFFT *ifft = TVirtualFFT::FFT(1, &M, "C2R EX K");
  ifft->SetPointsComplex(sRe.data(), sIm.data());
  ifft->Transform();

  // Extract result: keep only first N bins (linear part), normalize by FFT size
  TH1D *hOut = static_cast<TH1D *>(hData->Clone(outName));
  hOut->SetTitle(Form("Wiener deconv (α=%.3e);%s;Amplitude", dampingFactor, hData->GetXaxis()->GetTitle()));
  hOut->Reset();
  double norm = 1.0 / M;
  for (int i = 0; i < N; ++i) {
    hOut->SetBinContent(i + 1, std::max(0.0, ifft->GetPointReal(i) * norm));
  }
  delete ifft;

  // Normalize to unit integral
  double intg = hOut->Integral();
  if (intg > 0)
    hOut->Scale(1.0 / intg);

  return hOut;
}

// Build reconvolved spectrum from deconvolved signal and kernel.
// This mirrors the validation macro's discrete convolution convention so
// hReconvolved_iterX is directly consumable by validation without fallback.
static TH1D *reconvolveFromDeconvolved(const TH1D *hDeconv, const TH1D *hKernel, const char *outName) {
  if (hDeconv == nullptr || hKernel == nullptr)
    return nullptr;
  if (hDeconv->GetNbinsX() != hKernel->GetNbinsX())
    return nullptr;

  const int N = hDeconv->GetNbinsX();
  std::vector<double> s(N, 0.0), k(N, 0.0);
  for (int i = 0; i < N; ++i) {
    s[i] = std::max(0.0, hDeconv->GetBinContent(i + 1));
    k[i] = std::max(0.0, hKernel->GetBinContent(i + 1));
  }

  double ksum = 0.0;
  for (double v : k)
    ksum += v;
  if (ksum <= 0.0)
    return nullptr;
  for (double &v : k)
    v /= ksum;

  const double binWidth = hDeconv->GetXaxis()->GetBinWidth(1);
  const double xMinVal = hDeconv->GetXaxis()->GetXmin();
  int kc = static_cast<int>(std::round(-xMinVal / binWidth));
  kc = std::max(0, std::min(kc, N - 1));

  TH1D *out = dynamic_cast<TH1D *>(hDeconv->Clone(outName));
  if (out == nullptr)
    return nullptr;
  out->Reset();

  for (int i = 0; i < N; ++i) {
    double sum = 0.0;
    for (int j = 0; j < N; ++j) {
      const int ki = i - j + kc;
      if (ki >= 0 && ki < N)
        sum += s[j] * k[ki];
    }
    out->SetBinContent(i + 1, std::max(0.0, sum));
    out->SetBinError(i + 1, 0.0);
  }

  const double integral = out->Integral();
  if (integral > 0.0)
    out->Scale(1.0 / integral);

  return out;
}

// ====================================================================
//  Main program
// ====================================================================
int main(int argc, char **argv) {
  gStyle->SetOptStat(0);
  CommandLine CL(argc, argv);

  // Parse command-line arguments (mirroring Unfold.C input structure)
  std::string kernelFileName = CL.Get("KernelFile", "kernel.root");
  std::string dataFileName = CL.Get("DataFile", "data.root");
  std::string outputFileName = CL.Get("OutputFileName", "deconv_result.root");
  std::string varKernelName = CL.Get("VarKernelName", "HFEMaxPlus_forest");
  std::string varDataName = CL.Get("VarDataName", "HFEMaxPlus_forest");
  int iterations = CL.GetInt("Iterations", 5);  // Number of damping factor steps
  int dataQuarter = CL.GetInt("DataQuarter", 0);
  double xMin = CL.GetDouble("XMin", 0.0);
  double xMax = CL.GetDouble("XMax", 100.0);
  int binsPerGeV = CL.GetInt("BinsPerGeV", 2);
  double dampingMin = CL.GetDouble("DampingMin", 1e-3);
  double dampingMax = CL.GetDouble("DampingMax", 1e1);
  bool logDamping = CL.GetBool("LogDamping", true);  // Log-scale damping scan

  // Input validation
  if (xMax <= xMin) {
    std::cerr << "Invalid range: XMax must be greater than XMin." << std::endl;
    return -1;
  }
  if (binsPerGeV <= 0) {
    std::cerr << "Invalid binning: BinsPerGeV must be > 0." << std::endl;
    return -1;
  }
  if (iterations < 1) {
    std::cerr << "Iterations must be >= 1." << std::endl;
    return -1;
  }
  if(dataQuarter < -1 || dataQuarter > 3)
  {
    std::cerr << "Invalid DataQuarter: must be in [-1,3]." << std::endl;
    return -1;
  }
  if (dampingMin <= 0 || dampingMax <= 0 || dampingMin > dampingMax) {
    std::cerr << "Invalid damping range: both must be > 0 and DampingMin <= DampingMax." << std::endl;
    return -1;
  }

  const int nBins = static_cast<int>(std::lround((xMax - xMin) * binsPerGeV));
  if (nBins < 1) {
    std::cerr << "Computed number of bins is < 1." << std::endl;
    return -1;
  }

  printf("Config: xMin=%.1f  xMax=%.1f  nBins=%d  BinsPerGeV=%d\n", xMin, xMax, nBins, binsPerGeV);
  printf("Damping scan: α ∈ [%.3e, %.3e], %d iterations, %s scale\n", 
         dampingMin, dampingMax, iterations, logDamping ? "log" : "linear");

  // Load histograms from trees
  printf("\nReading kernel/noise from %s [%s/%s]...\n", kernelFileName.c_str(), "OutputTree", varKernelName.c_str());
  // The kernel histogram is used to build the response matrix, so we load it with the full dataset (dataQuarter = -1) to get the best possible statistics for the response. The signal+noise histogram is loaded with the specified quarter to simulate a realistic measurement scenario.
  TH1D *hKernel = RootIOUtils::LoadTreeBranchHistogramOrNull(
    kernelFileName, "OutputTree", varKernelName, "hKernel", "hKernel;" + varKernelName + ";Events",
    nBins, xMin, xMax, -1, true, true, "kernel file");
  if (!hKernel) {
    std::cerr << "Error: cannot load kernel histogram." << std::endl;
    return -1;
  }

  printf("Reading data/observed from %s [%s/%s]...\n", dataFileName.c_str(), "OutputTree", varDataName.c_str());
  // The signal+noise histogram is loaded with the specified quarter to simulate a realistic measurement scenario.
  TH1D *hData = RootIOUtils::LoadTreeBranchHistogramOrNull(
    dataFileName, "OutputTree", varDataName, "hData", "hData;" + varDataName + ";Events",
    nBins, xMin, xMax, dataQuarter, true, true, "data file");
  if (!hData) {
    std::cerr << "Error: cannot load data histogram." << std::endl;
    return -1;
  }

  printf("  Kernel entries: %.0f\n", hKernel->GetEntries());
  printf("  Data entries:   %.0f\n", hData->GetEntries());

  // Normalize both to unit area
  if (hKernel->Integral() > 0)
    hKernel->Scale(1.0 / hKernel->Integral());
  if (hData->Integral() > 0)
    hData->Scale(1.0 / hData->Integral());

  // Create output ROOT file
  TFile *outputFile = RootIOUtils::OpenFileOrNull(outputFileName, "RECREATE", "output file");
  if (outputFile == nullptr) {
    return -1;
  }
  TTimeStamp *currentTime = new TTimeStamp();
  GeneralInfoManager man(outputFile, "InfoDir", false);
  man.AddSourceFile(kernelFileName, currentTime);
  man.AddSourceFile(dataFileName, currentTime);
  man.AddCutParameter("Iterations", iterations, currentTime);
  man.AddCutParameter("DataQuarter", dataQuarter, currentTime);
  man.AddCutParameter("XMin", xMin, currentTime);
  man.AddCutParameter("XMax", xMax, currentTime);
  man.AddCutParameter("BinsPerGeV", binsPerGeV, currentTime);
  man.AddCutParameter("DampingMin", dampingMin, currentTime);
  man.AddCutParameter("DampingMax", dampingMax, currentTime);
  man.AddCutParameter("LogDamping", logDamping ? 1.0 : 0.0, currentTime);

  // Store input histograms
  hData->SetName("hMeasured");
  hKernel->SetName("hKernel");
  outputFile->cd();
  hData->Write();
  hKernel->Write();

  // Create regularization scan tree
  TTree *regTree = new TTree("RegularizationTree", "Damping factor scan results");
  int treeIter = 0;
  double treeDamping = 0;
  regTree->Branch("Iteration", &treeIter, "Iteration/I");
  regTree->Branch("DampingFactor", &treeDamping, "DampingFactor/D");

  // Wiener deconvolution sweep over damping factors
  printf("\nRunning Wiener FFT deconvolution sweep...\n");
  std::vector<TH1D *> deconvolvedHistograms;
  std::vector<double> dampingFactorsUsed;

  for (int iter = 1; iter <= iterations; ++iter) {
    // Compute damping factor (linear or logarithmic)
    double frac = (iterations > 1) ? static_cast<double>(iter - 1) / (iterations - 1) : 0.0;
    double dampingFactor;
    if (logDamping) {
      dampingFactor = std::exp(std::log(dampingMin) + frac * (std::log(dampingMax) - std::log(dampingMin)));
    } else {
      dampingFactor = dampingMin + frac * (dampingMax - dampingMin);
    }

    printf("  Iteration %d/%d: α = %.3e\n", iter, iterations, dampingFactor);

    // Run Wiener deconvolution
    const std::string deconvName = "hDeconvolved_iter" + std::to_string(iter);
    TH1D *hDeconv = wienerDeconvolveFFT(hData, hKernel, dampingFactor, deconvName.c_str());
    if (!hDeconv) {
      std::cerr << "Failed to deconvolve at iteration " << iter << std::endl;
      continue;
    }

    // Build and store reconvolved histogram with naming expected by validation.
    const std::string reconvName = "hReconvolved_iter" + std::to_string(iter);
    TH1D *hReconv = reconvolveFromDeconvolved(hDeconv, hKernel, reconvName.c_str());
    if (!hReconv) {
      std::cerr << "Failed to reconvolve at iteration " << iter << std::endl;
      delete hDeconv;
      continue;
    }

    hDeconv->SetDirectory(outputFile);
    hDeconv->Write();
    hReconv->SetDirectory(outputFile);
    hReconv->Write();
    deconvolvedHistograms.push_back(hDeconv);
    dampingFactorsUsed.push_back(dampingFactor);

    // Fill regularization tree
    treeIter = iter;
    treeDamping = dampingFactor;
    regTree->Fill();
  }

  // Write regularization tree
  outputFile->cd();
  regTree->Write();

  // Create quick-check canvas with first and last iterations
  if (deconvolvedHistograms.size() >= 1) {
    TCanvas *c1 = new TCanvas("c1", "Wiener FFT Deconvolution Result", 800, 800);
    c1->SetLogy();

    hData->SetMarkerStyle(20);
    hData->SetMarkerColor(kBlue + 1);
    hData->SetLineColor(kBlue + 1);
    hData->Draw("E");

    TH1D *hFirst = deconvolvedHistograms.front();
    hFirst->SetMarkerStyle(24);
    hFirst->SetMarkerColor(kRed + 1);
    hFirst->SetLineColor(kRed + 1);
    hFirst->Draw("E SAME");

    if (deconvolvedHistograms.size() > 1) {
      TH1D *hLast = deconvolvedHistograms.back();
      hLast->SetMarkerStyle(26);
      hLast->SetMarkerColor(kGreen + 2);
      hLast->SetLineColor(kGreen + 2);
      hLast->Draw("E SAME");

      TLegend *legend = new TLegend(0.55, 0.68, 0.88, 0.9);
      legend->AddEntry(hData, "Observed (convolved)", "lep");
      legend->AddEntry(hFirst, Form("Deconvolved (α=%.1e)", dampingFactorsUsed.front()), "lep");
      legend->AddEntry(hLast, Form("Deconvolved (α=%.1e)", dampingFactorsUsed.back()), "lep");
      legend->Draw();
    } else {
      TLegend *legend = new TLegend(0.55, 0.75, 0.88, 0.9);
      legend->AddEntry(hData, "Observed (convolved)", "lep");
      legend->AddEntry(hFirst, "Deconvolved", "lep");
      legend->Draw();
    }

    outputFile->cd();
    c1->Write();
  }

  man.SaveToFile();
  RootIOUtils::CloseAndDeleteFile(outputFile);
  delete currentTime;

  printf("\nOutput written to: %s\n", outputFileName.c_str());
  return 0;
}
