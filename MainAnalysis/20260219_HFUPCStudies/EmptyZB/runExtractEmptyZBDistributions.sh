#!/bin/bash

INPUTDATAPATH=/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_HiForest_260218_HIPhysicsRawPrime0-10_HIRun2025A_highrZB_399766.root

DATE=$(date +"%Y%m%d")

# Create output folder if it doesn't exist
mkdir -p EmptyZB/output_EmptyZB

# This script runs the ExtractEmptyZBDistributons.C macro to extract the empty bunch crossing distributions from the specified input file.
./ExtractEmptyZBDistributons --Input $INPUTDATAPATH --Output EmptyZB/output_EmptyZB/emptyZBOutput_${DATE}_Trigger0_nTrkFilter1_nVtxFilter1_ZDCFilter1.root --TriggerChoice 0 --nTrkFilter 1 --nVtxFilter 1 --ZDCFilter 1 &
sleep 1
./ExtractEmptyZBDistributons --Input $INPUTDATAPATH --Output EmptyZB/output_EmptyZB/emptyZBOutput_${DATE}_Trigger0_nTrkFilter1_nVtxFilter0_ZDCFilter1.root --TriggerChoice 0 --nTrkFilter 1 --nVtxFilter 0 --ZDCFilter 1 &
sleep 1
./ExtractEmptyZBDistributons --Input $INPUTDATAPATH --Output EmptyZB/output_EmptyZB/emptyZBOutput_${DATE}_Trigger-1_nTrkFilter1_nVtxFilter1_ZDCFilter1.root --TriggerChoice -1 --nTrkFilter 1 --nVtxFilter 1 --ZDCFilter 1 &
sleep 1
./ExtractEmptyZBDistributons --Input $INPUTDATAPATH --Output EmptyZB/output_EmptyZB/emptyZBOutput_${DATE}_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root --TriggerChoice -1 --nTrkFilter 1 --nVtxFilter 0 --ZDCFilter 1 &
wait