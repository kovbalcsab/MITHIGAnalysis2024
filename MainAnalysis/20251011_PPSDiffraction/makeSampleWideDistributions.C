#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Messenger.h"
#include "CommandLine.h"

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

    // Create output histograms
    TH1I* hSelectedEvents = new TH1I("hSelectedEvents", "Selected Events;Selection Step;Entries", 5, 0, 5);
    hSelectedEvents->GetXaxis()->SetBinLabel(1, "Total Events");
    hSelectedEvents->GetXaxis()->SetBinLabel(2, "After Baseline Selection");
    hSelectedEvents->GetXaxis()->SetBinLabel(3, "After ZB Trigger");
    hSelectedEvents->GetXaxis()->SetBinLabel(4, "After 1 Vtx");
    hSelectedEvents->GetXaxis()->SetBinLabel(5, "After HF+ & ZDC+ Cuts");

    TH1F* hHFPlus = new TH1F("hHFPlus", "HF+ Energy Distribution;HF+ Energy (GeV);Entries", 100, 0, 100);
    TH1F* hZDCPlus = new TH1F("hZDCPlus", "ZDC+ Energy Distribution;ZDC+ Energy (fC);Entries", 100, 0, 5000);
    TH1F* hHFPlusPPS = new TH1F("hHFPlusPPS", "HF+ Energy Distribution (PPS);HF+ Energy (GeV);Entries", 100, 0, 100);
    TH1F* hZDCPlusPPS = new TH1F("hZDCPlusPPS", "ZDC+ Energy Distribution (PPS);ZDC+ Energy (fC);Entries", 100, 0, 5000);
    TH1F* hTrackPt = new TH1F("hTrackPt", "Track p_{T} Distribution;Track p_{T} (GeV/c);Entries", 100, 0, 10);
    TH1F* hTrackEta = new TH1F("hTrackEta", "Track #eta Distribution;Track #eta;Entries", 100, -2.4, 2.4);
    TH1F* hTrackPtPPS = new TH1F("hTrackPtPPS", "Track p_{T} Distribution (PPS);Track p_{T} (GeV/c);Entries", 100, 0, 10);
    TH1F* hTrackEtaPPS = new TH1F("hTrackEtaPPS", "Track #eta Distribution (PPS);Track #eta;Entries", 100, -2.4, 2.4);
    TH1F* hMultiplicity = new TH1F("hMultiplicity", "Event Multiplicity;Number of Tracks;Entries", 50, 0, 50);
    TH1F* hMultiplicityPPS = new TH1F("hMultiplicityPPS", "Event Multiplicity (PPS);Number of Tracks;Entries", 50, 0, 50);
    TH2F* hTrackMomvsdedx = new TH2F("hTrackMomvsdedx", "Track Momentum vs dE/dx;Momentum (GeV/c);log(dE/dx (MeV/cm));Entries", 200, 0, 10, 200, 0, 10);

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

        // Select events with exactly 1 vertex
        if (M->nVtx != 1) {
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

        // Fill histograms
        hHFPlus->Fill(M->hiHFPlus_pfle1);
        hZDCPlus->Fill(M->ZDCsumPlus);
        int multiplicity = 0;
        for (size_t iTrkId = 0; iTrkId < M->trkPt->size(); ++iTrkId) {
            if (M->trkPt->at(iTrkId) < 0.05 || fabs(M->trkEta->at(iTrkId)) > 2.4 || !M->highPurity->at(iTrkId)) {
                continue;
            }
            hTrackPt->Fill(M->trkPt->at(iTrkId));
            hTrackEta->Fill(M->trkEta->at(iTrkId));
            multiplicity++;
            // Fill momentum vs dE/dx histogram
            hTrackMomvsdedx->Fill((M->trkPt->at(iTrkId))*TMath::CosH(M->trkEta->at(iTrkId)), TMath::Log(M->dedxAllLikelihood->at(iTrkId)));
        }
        hMultiplicity->Fill(multiplicity);

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
            hMultiplicityPPS->Fill(multiplicity);
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
    hMultiplicity->Write();
    hTrackMomvsdedx->Write();
    hHFPlusPPS->Write();
    hZDCPlusPPS->Write();
    hTrackPtPPS->Write();
    hTrackEtaPPS->Write();
    hMultiplicityPPS->Write();
    outFile->Close();

    return 0;
}
