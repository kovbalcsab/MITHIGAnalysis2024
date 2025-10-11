// File: loopTreeWithVector.C
// Run with  .x loopTreeWithVector.C  from a ROOT prompt
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>

void loop() {
  TFile *f = TFile::Open(
      "20250701_pO_PhysicsIonPhysics0_393952_PPS_ZBtriggerInc_new.root");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open file!\n";
    return;
  }

  TTree *t = nullptr;
  f->GetObject("Tree", t); // ← change tree name
  if (!t) {
    std::cerr << "Tree not found!\n";
    return;
  }

  bool passBaselineEventSelection = 0;
  std::vector<float> *trkPt = nullptr;  // vector branch
  std::vector<float> *trkEta = nullptr; // vector branch
  std::vector<float> *trkPhi = nullptr; // vector branch

  t->SetBranchAddress("passBaselineEventSelection",
                      &passBaselineEventSelection);
  t->SetBranchAddress("trkPt", &trkPt);
  t->SetBranchAddress("trkEta", &trkEta);
  t->SetBranchAddress("trkPhi", &trkPhi);

  const Long64_t nEntries = t->GetEntries();
  for (Long64_t i = 0; i < 1000; ++i) {
    t->GetEntry(i);
    if (!passBaselineEventSelection)
      continue; // filter events based on event selection

    for (size_t j = 0; j < trkPt->size(); ++j) {

      if (trkPt->at(j) < 3 || trkPt->at(j) > 5 || fabs(trkEta->at(j)) > 1)
        continue; // filter tracks based on pt
    }
    std::cout << "Event " << i << ": "
              << "passBaselineEventSelection = " << passBaselineEventSelection
              << std::endl;
  }

  //------------------------------------------------------------
  // 6. Clean up
  //------------------------------------------------------------
  f->Close();
}
