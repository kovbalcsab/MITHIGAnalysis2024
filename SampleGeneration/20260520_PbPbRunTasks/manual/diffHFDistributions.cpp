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

  for (int i = 0; i < defaultHists.size(); i++) {
    TH2D *diffHist = diffHistograms(newHists[i], defaultHists[i]);
    diffHists.push_back(diffHist);
  }

  TFile output(OutputName.c_str(), "RECREATE");
  output.cd();
  for (TH2D *hist : diffHists) {
    hist->Write();
  }
  output.Close();

  return 0;
}
