// File: loopTreeWithVector.C
// Usage from ROOT prompt: .x loopTreeWithVector.C
// or compiled:          .x loopTreeWithVector.C++

#include <TFile.h>
#include <TTree.h>
#include <algorithm> // for std::min
#include <iostream>
#include <vector>

void loop() {
  TFile *f = TFile::Open("20250701_pO_PhysicsIonPhysics0_393952_PPS_ZBtriggerInc_new.root");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open file!\n";
    return;
  }

  TTree *t = nullptr;
  f->GetObject("Tree", t); // ← update if your tree is named differently
  if (!t) {
    std::cerr << "Tree not found!\n";
    f->Close();
    return;
  }

  bool passBaselineEventSelection = false;

  std::vector<float> *trkPt = nullptr;
  std::vector<float> *trkEta = nullptr;
  std::vector<float> *trkPhi = nullptr;

  std::vector<bool> *highPurity = nullptr;

  std::vector<int> *PPS_station = nullptr;
  std::vector<int> *PPS_x = nullptr;
  std::vector<int> *PPS_y = nullptr;

  t->SetBranchAddress("passBaselineEventSelection", &passBaselineEventSelection);
  t->SetBranchAddress("trkPt", &trkPt);
  t->SetBranchAddress("trkEta", &trkEta);
  t->SetBranchAddress("trkPhi", &trkPhi);
  t->SetBranchAddress("highPurity", &highPurity);
  t->SetBranchAddress("PPS_station", &PPS_station);
  t->SetBranchAddress("PPS_x", &PPS_x);
  t->SetBranchAddress("PPS_y", &PPS_y);

  const Long64_t nEntries = t->GetEntries();
  const Long64_t nMax = std::min<Long64_t>(nEntries, 1000); // cap to 1000 as in your code

  for (Long64_t i = 0; i < nMax; ++i) {
    t->GetEntry(i);

    if (!passBaselineEventSelection) {
      continue;
    }
    std::cout << "Event " << i << ": passBaselineEventSelection = " << passBaselineEventSelection << "\n";

    for (size_t j = 0; j < trkPt->size(); ++j) {
      if (!highPurity->at(j))
        continue;

      std::cout << "  Track " << j << ": Pt = " << trkPt->at(j) << ", Eta = " << trkEta->at(j)
                << ", Phi = " << trkPhi->at(j) << "\n";
    }

    if (PPS_station && PPS_x && PPS_y) {
      const size_t nPPS = PPS_station->size();
      for (size_t k = 0; k < nPPS; ++k) {
        int st = PPS_station->at(k);
        int x = (k < PPS_x->size()) ? PPS_x->at(k) : 0;
        int y = (k < PPS_y->size()) ? PPS_y->at(k) : 0;
      }
    }
  }
  f->Close();
}

int main() {
  loop();
  return 0;
}
