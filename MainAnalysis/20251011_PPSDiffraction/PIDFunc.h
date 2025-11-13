#include <vector>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"

float GetPIDScore(
  float p,
  float dedx,
  TF1* fCenter,
  TF1* fSigmaLo,
  TF1* fSigmaHi,
  bool logdedx = true
) {
  float PIDScore;
  if (logdedx) dedx = TMath::Log(dedx);
  float dedxCenter = fCenter->Eval(p);
  if (dedx < dedxCenter)
    PIDScore = (dedx - dedxCenter) / (dedxCenter - fSigmaLo->Eval(p));
  else
    PIDScore = (dedx - dedxCenter) / (fSigmaHi->Eval(p) - dedxCenter);
  return PIDScore;
}

void SavePIDRoot(
  TF1* fdedxPionCenter,
  TF1* fdedxPionSigmaLo,
  TF1* fdedxPionSigmaHi,
  TF1* fdedxKaonCenter,
  TF1* fdedxKaonSigmaLo,
  TF1* fdedxKaonSigmaHi,
  TF1* fdedxProtCenter,
  TF1* fdedxProtSigmaLo,
  TF1* fdedxProtSigmaHi,
  TString output = "../../CommonCode/root/DzeroUPC_dedxMap.root"
) {
  TFile* fout = new TFile(output, "RECREATE");
  fout->cd();
  fdedxPionCenter->SetName("fdedxPionCenter");
  fdedxPionSigmaLo->SetName("fdedxPionSigmaLo");
  fdedxPionSigmaHi->SetName("fdedxPionSigmaHi");
  fdedxKaonCenter->SetName("fdedxKaonCenter");
  fdedxKaonSigmaLo->SetName("fdedxKaonSigmaLo");
  fdedxKaonSigmaHi->SetName("fdedxKaonSigmaHi");
  fdedxProtCenter->SetName("fdedxProtCenter");
  fdedxProtSigmaLo->SetName("fdedxProtSigmaLo");
  fdedxProtSigmaHi->SetName("fdedxProtSigmaHi");
  fdedxPionCenter->Write();
  fdedxPionSigmaLo->Write();
  fdedxPionSigmaHi->Write();
  fdedxKaonCenter->Write();
  fdedxKaonSigmaLo->Write();
  fdedxKaonSigmaHi->Write();
  fdedxProtCenter->Write();
  fdedxProtSigmaLo->Write();
  fdedxProtSigmaHi->Write();
}

std::vector<TF1*> ImportPIDRoot(
  TString input = "../../CommonCode/root/DzeroUPC_dedxMap.root"
) {
  TFile* dedxInput = new TFile(input, "READ");
  TF1* fdedxPionCenter  = (TF1*) dedxInput->Get("fdedxPionCenter");
  TF1* fdedxPionSigmaLo = (TF1*) dedxInput->Get("fdedxPionSigmaLo");
  TF1* fdedxPionSigmaHi = (TF1*) dedxInput->Get("fdedxPionSigmaHi");
  TF1* fdedxKaonCenter  = (TF1*) dedxInput->Get("fdedxKaonCenter");
  TF1* fdedxKaonSigmaLo = (TF1*) dedxInput->Get("fdedxKaonSigmaLo");
  TF1* fdedxKaonSigmaHi = (TF1*) dedxInput->Get("fdedxKaonSigmaHi");
  TF1* fdedxProtCenter  = (TF1*) dedxInput->Get("fdedxProtCenter");
  TF1* fdedxProtSigmaLo = (TF1*) dedxInput->Get("fdedxProtSigmaLo");
  TF1* fdedxProtSigmaHi = (TF1*) dedxInput->Get("fdedxProtSigmaHi");
  std::vector<TF1*> dedxFunctions = {
    fdedxPionCenter, fdedxPionSigmaLo, fdedxPionSigmaHi,
    fdedxKaonCenter, fdedxKaonSigmaLo, fdedxKaonSigmaHi,
    fdedxProtCenter, fdedxProtSigmaLo, fdedxProtSigmaHi
  };
  dedxInput->Close();
  delete dedxInput;
  return dedxFunctions;
}

const float KMatchDtrkPMax = 0.9;
const float KMatchSigmaMax = 3.0;

bool matchKaon(
  float trkP,
  float kaonScore,
  float kaonScoreMin = -KMatchSigmaMax,
  float kaonScoreMax = KMatchSigmaMax,
  float KaonMaxP = KMatchDtrkPMax
) {
  return (
    (kaonScore > kaonScoreMin &&
     kaonScore < kaonScoreMax &&
     trkP < KaonMaxP) ||
    trkP > KaonMaxP
  );
}

bool ApplyPIDCut(
  float trk1P,
  float trk1PionScore,
  float trk1KaonScore,
  float trk1ProtScore,
  float trk2P,
  float trk2PionScore,
  float trk2KaonScore,
  float trk2ProtScore,
  bool matchBothKaons = false,
  bool doProtonReject = false,
  float pionScoreMin = -3,
  float pionScoreMax = 3,
  float kaonScoreMin = -KMatchSigmaMax,
  float kaonScoreMax = KMatchSigmaMax,
  float protScoreMin = -1.,
  float protScoreMax = 999.,
  float PionMaxP = 0.85,
  float KaonMaxP = 1.0,
  float ProtMinP = 0.3,
  float ProtMaxP = 1.1
) {
  bool matchKaon0 = false;
  bool matchKaon1 = false;
  bool matchProton = false;
    
  //Match track1 to kaon
  matchKaon0 = (
  trk1KaonScore > kaonScoreMin &&
  trk1KaonScore < kaonScoreMax &&
  trk1P < KaonMaxP) ||
  trk1P > KaonMaxP;

  //Match track2 to kaon
  matchKaon1 = (
  trk2KaonScore > kaonScoreMin &&
  trk2KaonScore < kaonScoreMax &&
  trk2P < KaonMaxP) ||
  trk2P > KaonMaxP;

  // proton rejection
  if (doProtonReject &&
      ((trk1P > ProtMinP &&
        trk1P < ProtMaxP &&
        trk1ProtScore > protScoreMin &&
        trk1ProtScore < protScoreMax) ||
       (trk2P > ProtMinP &&
        trk2P < ProtMaxP &&
        trk2ProtScore > protScoreMin &&
        trk2ProtScore < protScoreMax))
    ) matchProton = true;

  if (matchBothKaons) {
    return (matchKaon0 && matchKaon1 && !matchProton);
  } else {
    return ((matchKaon0 || matchKaon1) && !matchProton);
  }
}