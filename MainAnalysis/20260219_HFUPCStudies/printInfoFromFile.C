#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>
#include <string>

using namespace std;
#include "CommandLine.h" // Yi's Commandline bundle
#include "Messenger.h"   // Yi's Messengers for reading data files
#include "ProgressBar.h" // Yi's fish progress bar
#include "utilities.h"   // Yen-Jie's random utility functions
#include "InfoManager.h"

int main(int argc, char *argv[]){
    CommandLine CL(argc, argv);
    string InputFileName = CL.Get("Input");

    TFile* inputFile = TFile::Open(InputFileName.c_str());
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "\033[1;31mError: Could not open file " << InputFileName << "\033[0m" << std::endl;
        return 1;
    }
    GeneralInfoManager man(inputFile, "InfoDir", true);

    std::cout << "Opened file: " << InputFileName << std::endl;
    if (!man.IsInitialized()) {
        std::cerr << "\033[1;31mError: GeneralInfoManager is not initialized properly.\033[0m" << std::endl;
        return 1;
    }

    if (man.NoDuplicates()) {
        std::cout << "\033[1;32mNo duplicates found in the information.\033[0m" << std::endl;
    } else {
        std::cout << "\033[1;31mDuplicates found in the information.\033[0m" << std::endl;
    }

    if (man.CutConsistencyCheck()) {
        std::cout << "\033[1;32mCut parameters are consistent.\033[0m" << std::endl;
    } else {
        std::cout << "\033[1;31mInconsistent cut parameters found.\033[0m" << std::endl;
    }

    man.PrintInfo();

    inputFile->Close();
    return 0;
}