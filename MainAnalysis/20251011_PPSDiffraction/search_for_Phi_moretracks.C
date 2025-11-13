#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Messenger.h"
#include "CommandLine.h"
#include "PIDFunc.h"

int main(int argc, char** argv) {
    // Parse command line arguments
    CommandLine cl(argc, argv);
    std::string inputFile = cl.Get("inputFile", "data.root");
    std::string outputFile = cl.Get("outputFile", "output.root");

    // Open input ROOT file
    TFile* inFile = TFile::Open(inputFile.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return 1;
    }

    std::vector<TF1*> dedxFunctions = ImportPIDRoot("DzeroUPC_dedxMap.root");
    TF1* fdedxPionCenter  = dedxFunctions[0];
    TF1* fdedxPionSigmaLo = dedxFunctions[1];
    TF1* fdedxPionSigmaHi = dedxFunctions[2];
    TF1* fdedxKaonCenter  = dedxFunctions[3];
    TF1* fdedxKaonSigmaLo = dedxFunctions[4];
    TF1* fdedxKaonSigmaHi = dedxFunctions[5];
    TF1* fdedxProtCenter  = dedxFunctions[6];
    TF1* fdedxProtSigmaLo = dedxFunctions[7];
    TF1* fdedxProtSigmaHi = dedxFunctions[8];

    // Create output histograms
    TH1D* hMass = new TH1D("hMass", "Two track Mass;Mass (GeV/c^{2});Entries", 480, 0.8, 2);
    TH1D* hMassLowPt = new TH1D("hMassLowPt", "Two track Mass;Mass (GeV/c^{2});Entries", 200, 0.98, 1.15);
    TH1D* hMassPPS = new TH1D("hMassPPS", "Two track Mass (PPS);Mass (GeV/c^{2});Entries", 480, 0.8, 2);
    TH1D* hPhiPt = new TH1D("hPhiPt", "Phi p_{T};p_{T} (GeV/c);Entries", 100, 0, 10);
    TH1D* hPhiPtPPS = new TH1D("hPhiPtPPS", "Phi p_{T} (PPS);p_{T} (GeV/c);Entries", 100, 0, 10);
    TH1D* hPhiEta = new TH1D("hPhiEta", "Phi #eta;#eta;Entries", 100, -2.4, 2.4);
    TH1D* hPhiEtaPPS = new TH1D("hPhiEtaPPS", "Phi #eta (PPS);#eta;Entries", 100, -2.4, 2.4);
    TH1D* hPhiPhi = new TH1D("hPhiPhi", "Phi #phi;#phi;Entries", 64, -3.2, 3.2);
    TH1D* hPhiPhiPPS = new TH1D("hPhiPhiPPS", "Phi #phi (PPS);#phi;Entries", 64, -3.2, 3.2);


    TH1I* hSelectedEvents = new TH1I("hSelectedEvents", "Selected Events;Selection Step;Entries", 7, 0, 7);
    hSelectedEvents->GetXaxis()->SetBinLabel(1, "Total Events");
    hSelectedEvents->GetXaxis()->SetBinLabel(2, "After Baseline Selection");
    hSelectedEvents->GetXaxis()->SetBinLabel(3, "After ZB Trigger");
    hSelectedEvents->GetXaxis()->SetBinLabel(4, "After 2 Track & 1 Vtx");
    hSelectedEvents->GetXaxis()->SetBinLabel(5, "After HF+ & ZDC+ Cuts");
    hSelectedEvents->GetXaxis()->SetBinLabel(6, "Opposite Charged Tracks");
    hSelectedEvents->GetXaxis()->SetBinLabel(7, "PID Selection");

    TH1F* hHFPlus = new TH1F("hHFPlus", "HF+ Energy Distribution;HF+ Energy (GeV);Entries", 100, 0, 100);
    TH1F* hZDCPlus = new TH1F("hZDCPlus", "ZDC+ Energy Distribution;ZDC+ Energy (fC);Entries", 100, 0, 5000);
    TH1F* hHFPlusPPS = new TH1F("hHFPlusPPS", "HF+ Energy Distribution (PPS);HF+ Energy (GeV);Entries", 100, 0, 100);
    TH1F* hZDCPlusPPS = new TH1F("hZDCPlusPPS", "ZDC+ Energy Distribution (PPS);ZDC+ Energy (fC);Entries", 100, 0, 5000);
    TH1F* hTrackPt = new TH1F("hTrackPt", "Track p_{T} Distribution;Track p_{T} (GeV/c);Entries", 100, 0, 10);
    TH1F* hTrackEta = new TH1F("hTrackEta", "Track #eta Distribution;Track #eta;Entries", 100, -2.4, 2.4);
    TH1F* hTrackPtPPS = new TH1F("hTrackPtPPS", "Track p_{T} Distribution (PPS);Track p_{T} (GeV/c);Entries", 100, 0, 10);
    TH1F* hTrackEtaPPS = new TH1F("hTrackEtaPPS", "Track #eta Distribution (PPS);Track #eta;Entries", 100, -2.4, 2.4);

    TH2F* hTrackMomvsdedx = new TH2F("hTrackMomvsdedx", "Track Momentum vs dE/dx;Momentum (GeV/c);log(dE/dx (MeV/cm));Entries", 200, 0, 10, 200, 0, 10);
    TH2F* hTrackMomvsdedxK = new TH2F("hTrackMomvsdedxK", "Track Momentum vs dE/dx;Momentum (GeV/c);log(dE/dx (MeV/cm));Entries", 200, 0, 10, 200, 0, 10);
    TH2F* hTrackMomvsdedxPi = new TH2F("hTrackMomvsdedxPi", "Track Momentum vs dE/dx;Momentum (GeV/c);log(dE/dx (MeV/cm));Entries", 200, 0, 10, 200, 0, 10);

    // Initialize messenger for diffractive analysis
    pODiffractiveTreeMessenger *M = new pODiffractiveTreeMessenger(inFile, "Tree");

    int nEntries = M->GetEntries();
    float HF_threshold = 9.0; // is this good?
    int ZDC_threshold = 1300;

    for (int i = 0; i < nEntries; ++i) {
        M->GetEntry(i);
        
        // Count total events
        //if (M->Run != 393952) {
        //    continue; // Process only specific run
        //}
        hSelectedEvents->Fill(0.5);
        
        // Apply baseline event selection (VZ, ClusterCompatibilityFilter, PVFilter)
        if (!M->passBaselineEventSelection) {
            continue;
        }
        hSelectedEvents->Fill(1.5);

        // Apply ZB trigger selection
        if (!M->HLT_OxyZeroBias_v1) {
            continue;
        }
        hSelectedEvents->Fill(2.5);

        // Select events with only 2 tracks in |eta| < 2.4 && pT > 0.05 GeV/c && highPurity
        // Select events with exactly 1 vertex
        if (M->multiplicityEta2p4 != 2 || M->nVtx != 1) {
            continue;
        }
        hSelectedEvents->Fill(3.5);

        // Apply HF+ cut
        if (M->hiHFPlus_pfle1 > HF_threshold) {
            continue;
        }

        // Apply ZDC+ cut
        if (M->ZDCsumPlus > ZDC_threshold) {
            continue;
        }
        hSelectedEvents->Fill(4.5);

        // Calculate two-track invariant mass and check for opposite charge
        ROOT::Math::PtEtaPhiMVector kaons[2];
        int kaonIdxs[2];
        int chargeProduct = 1;
        int iKaon = 0;
        for (size_t iTrkId = 0; iTrkId < M->trkPt->size(); ++iTrkId) {
            if (M->trkPt->at(iTrkId) < 0.05 || fabs(M->trkEta->at(iTrkId)) > 2.4 || !M->highPurity->at(iTrkId)) {
                continue;
            }
            kaons[iKaon] = ROOT::Math::PtEtaPhiMVector(M->trkPt->at(iTrkId), M->trkEta->at(iTrkId), M->trkPhi->at(iTrkId), 0.493677);
            if (M->trkPt->at(iTrkId)*TMath::CosH(M->trkEta->at(iTrkId)) != kaons[iKaon].P()) {
                std::cerr << "Error: Momentum calculation mismatch." << std::endl;
            }
            kaonIdxs[iKaon] = iTrkId;
            chargeProduct *= M->trkCharge->at(iTrkId);
            iKaon++;
        }
        if (iKaon != 2) {
            std::cerr << "Error: Expected 2 kaons, found " << iKaon << std::endl;
            continue;
        }
        if (chargeProduct != -1) {
            continue; // Require opposite charged tracks
        }
        hSelectedEvents->Fill(5.5);

        // PID Selection based on dE/dx scores
        float PIDScoreKaon0 = GetPIDScore(kaons[0].P(),
                                          M->dedxAllLikelihood->at(kaonIdxs[0]),
                                          fdedxKaonCenter, fdedxKaonSigmaLo, fdedxKaonSigmaHi);
        float PIDScoreKaon1 = GetPIDScore(kaons[1].P(),
                                          M->dedxAllLikelihood->at(kaonIdxs[1]),
                                          fdedxKaonCenter, fdedxKaonSigmaLo, fdedxKaonSigmaHi);
        float PIDScorePion0 = GetPIDScore(kaons[0].P(),
                                          M->dedxAllLikelihood->at(kaonIdxs[0]),
                                          fdedxPionCenter, fdedxPionSigmaLo, fdedxPionSigmaHi);
        float PIDScorePion1 = GetPIDScore(kaons[1].P(),
                                          M->dedxAllLikelihood->at(kaonIdxs[1]),
                                          fdedxPionCenter, fdedxPionSigmaLo, fdedxPionSigmaHi);
        float PIDScoreProt0 = GetPIDScore(kaons[0].P(),
                                          M->dedxAllLikelihood->at(kaonIdxs[0]),
                                          fdedxProtCenter, fdedxProtSigmaLo, fdedxProtSigmaHi);
        float PIDScoreProt1 = GetPIDScore(kaons[1].P(),
                                          M->dedxAllLikelihood->at(kaonIdxs[1]),
                                          fdedxProtCenter, fdedxProtSigmaLo, fdedxProtSigmaHi);
        

        bool passPID = ApplyPIDCut(
            kaons[0].P(), PIDScorePion0, PIDScoreKaon0, PIDScoreProt0,
            kaons[1].P(), PIDScorePion1, PIDScoreKaon1, PIDScoreProt1,
            true
        );
        if (!passPID) {
            continue;
        }
        hSelectedEvents->Fill(6.5);

        // Calculate di-kaon system
        auto diKaon = kaons[0] + kaons[1];
        
        hMass->Fill(diKaon.M());
        hPhiPt->Fill(diKaon.Pt());
        hPhiEta->Fill(diKaon.Eta());
        hPhiPhi->Fill(diKaon.Phi());

        if (diKaon.Pt() < 0.2) {
            hMassLowPt->Fill(diKaon.M());
        }
        hTrackMomvsdedx->Fill(kaons[0].P(), TMath::Log(M->dedxAllLikelihood->at(kaonIdxs[0])));
        hTrackMomvsdedx->Fill(kaons[1].P(), TMath::Log(M->dedxAllLikelihood->at(kaonIdxs[1])));

        if (diKaon.M() < 1.1 && diKaon.M() > 0.98 && diKaon.Pt() < 0.2) {
            hTrackMomvsdedxK->Fill(kaons[0].P(), TMath::Log(M->dedxAllLikelihood->at(kaonIdxs[0])));
            hTrackMomvsdedxK->Fill(kaons[1].P(), TMath::Log(M->dedxAllLikelihood->at(kaonIdxs[1])));
        } else if (diKaon.M() > 1.1 && diKaon.M() < 1.3 && diKaon.Pt() < 0.2) {
            hTrackMomvsdedxPi->Fill(kaons[0].P(), TMath::Log(M->dedxAllLikelihood->at(kaonIdxs[0])));
            hTrackMomvsdedxPi->Fill(kaons[1].P(), TMath::Log(M->dedxAllLikelihood->at(kaonIdxs[1])));
        }

        // Fill histograms
        hHFPlus->Fill(M->hiHFPlus_pfle1);
        hZDCPlus->Fill(M->ZDCsumPlus);
        for (size_t iTrkId = 0; iTrkId < M->trkPt->size(); ++iTrkId) {
            if (M->trkPt->at(iTrkId) < 0.05 || fabs(M->trkEta->at(iTrkId)) > 2.4 || !M->highPurity->at(iTrkId)) {
                continue;
            }
            hTrackPt->Fill(M->trkPt->at(iTrkId));
            hTrackEta->Fill(M->trkEta->at(iTrkId));
        }

        bool passPPS = (M->PPSStation0M_x->size() > 0 && M->PPSStation2M_x->size() > 0 && M->PPSStation0M_y->size() > 0 && M->PPSStation2M_y->size() > 0);

        if (passPPS) {
            hHFPlusPPS->Fill(M->hiHFPlus_pfle1);
            hZDCPlusPPS->Fill(M->ZDCsumPlus);
            for (size_t iTrkId = 0; iTrkId < M->trkPt->size(); ++iTrkId) {
                if (M->trkPt->at(iTrkId) < 0.05 || fabs(M->trkEta->at(iTrkId)) > 2.4 || !M->highPurity->at(iTrkId)) {
                    continue;
                }
                hTrackPtPPS->Fill(M->trkPt->at(iTrkId));
                hTrackEtaPPS->Fill(M->trkEta->at(iTrkId));
            }
            hMassPPS->Fill(diKaon.M());
            hPhiPtPPS->Fill(diKaon.Pt());
            hPhiEtaPPS->Fill(diKaon.Eta());
            hPhiPhiPPS->Fill(diKaon.Phi());
        }
    }

    // Write output histograms to file
    TFile* outFile = TFile::Open(outputFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFile << std::endl;
        return 1;
    }

    // Write histograms to output file
    outFile->cd();
    hSelectedEvents->Write();
    hHFPlus->Write();
    hZDCPlus->Write();
    hTrackPt->Write();
    hTrackEta->Write();
    hMass->Write();
    hPhiPt->Write();
    hPhiEta->Write();
    hPhiPhi->Write();
    hHFPlusPPS->Write();
    hZDCPlusPPS->Write();
    hTrackPtPPS->Write();
    hTrackEtaPPS->Write();
    hMassPPS->Write();
    hPhiPtPPS->Write();
    hPhiEtaPPS->Write();
    hPhiPhiPPS->Write();
    hMassLowPt->Write();
    hTrackMomvsdedx->Write();
    hTrackMomvsdedxK->Write();
    hTrackMomvsdedxPi->Write();
    fdedxKaonCenter->Write();
    fdedxKaonSigmaLo->Write();
    fdedxKaonSigmaHi->Write();
    fdedxPionCenter->Write();
    fdedxPionSigmaLo->Write();
    fdedxPionSigmaHi->Write();
    fdedxProtCenter->Write();
    fdedxProtSigmaLo->Write();
    fdedxProtSigmaHi->Write();
    outFile->Close();

    return 0;
}
