#!/bin/bash

rm ExecuteScanPFTowerLocations
rm ExecuteMakeHFPFBinEdges
rm ExecuteCheckHFPFBinning
rm ExecutePlotHFDistributionsSides
rm ExecuteHFDistributionExtraction
rm ExecuteRecalculateMeanAndStd
rm ExecuteDiffHFDistributions
rm ExtractEmptyBXFullMaps
rm PlotForestEtaPhiMaps

rm -rf ../../CommonCode/binary/
rm -rf ../../CommonCode/library/

#cd /home/$USER/CMSSW_13_2_13/src
#cd /home/$USER/CMSSW_14_1_7/src
cd /home/$USER/CMSSW_15_1_0_patch5/src
cmsenv

cd -
echo "CMSSW environment is set up"
cd ../../
source SetupAnalysis.sh
cd CommonCode/
make
cd ..
cd SampleGeneration/20260520_PbPbRunTasks
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/RooUnfold
make
