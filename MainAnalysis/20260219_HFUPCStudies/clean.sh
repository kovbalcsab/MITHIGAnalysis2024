rm ExtractEmptyBXDistributons
rm ExtractEmptyZBDistributons
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
make
