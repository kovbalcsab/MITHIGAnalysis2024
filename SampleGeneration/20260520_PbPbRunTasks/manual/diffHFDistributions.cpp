#include "utility.h"

int main(int argc, char *argv[]) {
  CommandLine CL(argc, argv);
  std::string DefaultInputName = CL.Get("DefaultInput", "default.root");
  std::string newInputName = CL.Get("NewInput", "new.root");
  std::string OutputName = CL.Get("Output", "diff.root");

  TFile defaultInput(DefaultInputName.c_str(), "READ");
  if (defaultInput.IsOpen() == false) {
    std::cerr << "Failed to open default input file: " << DefaultInputName << std::endl;
    return -1;
  }

  TFile newInput(newInputName.c_str(), "READ");
  if (newInput.IsOpen() == false) {
    std::cerr << "Failed to open new input file: " << newInputName << std::endl;
    return -1;
  }

  std::vector<TH2D *> defaultHists = loadMetaHistograms(defaultInput);
  std::vector<TH2D *> newHists = loadMetaHistograms(newInput);
  std::vector<TH2D *> diffHists;

  TH2D *H36Negative_leadingProb_def = (TH2D *)defaultInput.Get("h36Negative_leadingProb");
  TH2D *H36Positive_leadingProb_def = (TH2D *)defaultInput.Get("h36Positive_leadingProb");
  TH2D *H18NegativeMain_leadingProb_def = (TH2D *)defaultInput.Get("h18NegativeMain_leadingProb");
  TH2D *H18NegativeWrapLow_leadingProb_def = (TH2D *)defaultInput.Get("h18NegativeWrapLow_leadingProb");
  TH2D *H18NegativeWrapHigh_leadingProb_def = (TH2D *)defaultInput.Get("h18NegativeWrapHigh_leadingProb");
  TH2D *H18PositiveMain_leadingProb_def = (TH2D *)defaultInput.Get("h18PositiveMain_leadingProb");
  TH2D *H18PositiveWrapLow_leadingProb_def = (TH2D *)defaultInput.Get("h18PositiveWrapLow_leadingProb");
  TH2D *H18PositiveWrapHigh_leadingProb_def = (TH2D *)defaultInput.Get("h18PositiveWrapHigh_leadingProb");
  std::vector<TH2D *> posHists_leadingProb_def = {H36Positive_leadingProb_def, H18PositiveMain_leadingProb_def,
                                                  H18PositiveWrapLow_leadingProb_def,
                                                  H18PositiveWrapHigh_leadingProb_def};
  std::vector<TH2D *> negHists_leadingProb_def = {H36Negative_leadingProb_def, H18NegativeMain_leadingProb_def,
                                                  H18NegativeWrapLow_leadingProb_def,
                                                  H18NegativeWrapHigh_leadingProb_def};
  TH2D *H36Negative_leadingProb_new = (TH2D *)newInput.Get("h36Negative_leadingProb");
  TH2D *H36Positive_leadingProb_new = (TH2D *)newInput.Get("h36Positive_leadingProb");
  TH2D *H18NegativeMain_leadingProb_new = (TH2D *)newInput.Get("h18NegativeMain_leadingProb");
  TH2D *H18NegativeWrapLow_leadingProb_new = (TH2D *)newInput.Get("h18NegativeWrapLow_leadingProb");
  TH2D *H18NegativeWrapHigh_leadingProb_new = (TH2D *)newInput.Get("h18NegativeWrapHigh_leadingProb");
  TH2D *H18PositiveMain_leadingProb_new = (TH2D *)newInput.Get("h18PositiveMain_leadingProb");
  TH2D *H18PositiveWrapLow_leadingProb_new = (TH2D *)newInput.Get("h18PositiveWrapLow_leadingProb");
  TH2D *H18PositiveWrapHigh_leadingProb_new = (TH2D *)newInput.Get("h18PositiveWrapHigh_leadingProb");
  std::vector<TH2D *> posHists_leadingProb_new = {H36Positive_leadingProb_new, H18PositiveMain_leadingProb_new,
                                                  H18PositiveWrapLow_leadingProb_new,
                                                  H18PositiveWrapHigh_leadingProb_new};
  std::vector<TH2D *> negHists_leadingProb_new = {H36Negative_leadingProb_new, H18NegativeMain_leadingProb_new,
                                                  H18NegativeWrapLow_leadingProb_new,
                                                  H18NegativeWrapHigh_leadingProb_new};

  for (int i = 0; i < defaultHists.size(); i++) {
    TH2D *diffHist = diffHistograms(newHists[i], defaultHists[i]);
    diffHists.push_back(diffHist);
  }
  for (int i = 0; i < posHists_leadingProb_def.size(); i++) {
    diffHists.push_back(diffHistograms(posHists_leadingProb_new[i], posHists_leadingProb_def[i]));
    diffHists.push_back(diffHistograms(negHists_leadingProb_new[i], negHists_leadingProb_def[i]));
  }

  TFile output(OutputName.c_str(), "RECREATE");
  output.cd();
  for (TH2D *hist : diffHists) {
    hist->Write();
  }
  output.Close();

  return 0;
}
