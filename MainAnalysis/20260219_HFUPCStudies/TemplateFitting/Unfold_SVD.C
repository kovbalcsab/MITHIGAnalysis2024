#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMatrixD.h"

#include <iostream>
#include <cmath>
#include <string>

#include "CommandLine.h"

static TH1D *treeToHist(const std::string &fileName, const std::string &treeName, const std::string &varName,
                        const std::string &histName, int nBins, double xMin, double xMax) {
  TFile *f = TFile::Open(fileName.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << fileName << std::endl;
    return nullptr;
  }
  TTree *t = dynamic_cast<TTree *>(f->Get(treeName.c_str()));
  if (!t) {
    std::cerr << "Cannot find tree '" << treeName << "' in " << fileName << std::endl;
    f->Close();
    return nullptr;
  }

  // Check the branch exists
  if (!t->GetBranch(varName.c_str())) {
    std::cerr << "Cannot find branch '" << varName << "' in tree '" << treeName << "'\n";
    std::cerr << "Available branches: ";
    TObjArray *branches = t->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
      TObject *obj = branches->At(i);
      if (obj)
        std::cerr << obj->GetName() << (i + 1 < branches->GetEntries() ? ", " : "\n");
    }
    f->Close();
    return nullptr;
  }

  TH1D *h = new TH1D(histName.c_str(), (histName + ";" + varName + ";Events").c_str(), nBins, xMin, xMax);
  h->SetDirectory(nullptr);

  // Use SetBranchAddress + manual loop — works regardless of ROOT directory state
  double val = 0;
  float valf = 0;
  bool isFloat = false;

  TBranch *br = t->GetBranch(varName.c_str());
  TLeaf *lf = br->GetLeaf(varName.c_str());
  if (!lf)
    lf = static_cast<TLeaf *>(br->GetListOfLeaves()->First());
  std::string typeName = lf ? std::string(lf->GetTypeName()) : "";
  isFloat = (typeName == "Float_t" || typeName == "float");

  t->SetBranchStatus("*", 0);             // disable all branches for speed
  t->SetBranchStatus(varName.c_str(), 1); // enable only what we need

  if (isFloat) {
    t->SetBranchAddress(varName.c_str(), &valf);
  } else {
    t->SetBranchAddress(varName.c_str(), &val);
  }

  Long64_t nEntries = t->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    t->GetEntry(i);
    h->Fill(isFloat ? static_cast<double>(valf) : val);
  }

  printf("  Loaded %lld entries from %s/%s/%s\n", nEntries, fileName.c_str(), treeName.c_str(), varName.c_str());

  t->SetBranchStatus("*", 1); // re-enable all branches
  t->ResetBranchAddresses();
  f->Close();
  delete f;
  return h;
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

    const int requestedNBins = static_cast<int>(std::lround((xMax - xMin) * binsPerGeV));
    if(requestedNBins < 1)
    {
      std::cerr << "Invalid binning: computed number of bins is < 1." << std::endl;
      return -1;
    }

    TH1D* noiseHist = treeToHist(noiseFileName, "OutputTree", varNoiseName, "hNoise", requestedNBins, xMin, xMax);
    TH1D* signalPlusNoiseHist = treeToHist(DataFileName, "OutputTree", varDataName, "hData", requestedNBins, xMin, xMax);

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

    TFile *outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
    if(!outputFile || outputFile->IsZombie())
    {
      std::cerr << "Cannot create output file: " << outputFileName << std::endl;
      return -1;
    }

    signalPlusNoiseHist->SetName("hMeasured");
    noiseHist->SetName("hNoisePDF");
    signalPlusNoiseHist->Write();
    noiseHist->Write();

    TTree *regTree = new TTree("RegularizationTree", "SVD regularization parameter scan");
    int regIteration = 0;
    int regKTerm = 0;
    regTree->Branch("Iteration", &regIteration, "Iteration/I");
    regTree->Branch("KTerm", &regKTerm, "KTerm/I");

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

      regIteration = iter;
      regKTerm = kTerm;
      regTree->Fill();

      std::cout << "Iteration " << iter << " uses k-term = " << kTerm << std::endl;

      const std::string unfoldedName = "hUnfolded_iter" + std::to_string(iter);
      TH1D *unfoldedSignal = dynamic_cast<TH1D *>(unfoldedTmp->Clone(unfoldedName.c_str()));
      if(!unfoldedSignal)
      {
        std::cerr << "Failed to clone unfolded histogram for iteration " << iter << std::endl;
        continue;
      }
      unfoldedSignal->SetDirectory(outputFile);
      unfoldedSignal->Write();

      const std::string refoldedName = "hRefolded_iter" + std::to_string(iter);
      TH1 *refoldedSignal = response.ApplyToTruth(unfoldedSignal, refoldedName.c_str());
      if(!refoldedSignal)
      {
        std::cerr << "Failed to build refolded histogram for iteration " << iter << std::endl;
        continue;
      }
      refoldedSignal->SetDirectory(outputFile);
      refoldedSignal->Write();

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

    outputFile->Close();

    return 0;
}