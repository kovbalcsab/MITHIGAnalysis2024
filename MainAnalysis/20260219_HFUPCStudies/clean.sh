rm ExtractEmptyBXDistributons
rm ExtractEmptyZBDistributons
rm PrintInfoFromFile
rm ExtractSignalMCDistributons
rm PlotComparisonMultiVar
rm PlotForestEtaPhiMaps
rm ScalingFit
rm ScalingCDF
rm plotScalingResult
rm PlotCDFResult
rm ExtractEmptyBXFullMaps
rm -rf ../../CommonCode/binary/
rm -rf ../../CommonCode/library/
rm .DS_Store

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
cd MainAnalysis/20260219_HFUPCStudies/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/RooUnfold
make
