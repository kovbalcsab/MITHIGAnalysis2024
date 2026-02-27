#!/bin/bash

INPUTDATAPATH=/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1.root

DATE=$(date +"%Y%m%d")

# Create output folder if it doesn't exist
mkdir -p EmptyBX/output_EmptyBX

# This script runs the ExtractEmptyBXDistributons.C macro to extract the empty bunch crossing distributions from the specified input file.
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger0_nTrkFilter0.root --TriggerChoice 0 --nTrkFilter 0 &
sleep 1 # Sleep for a short time to avoid potential issues with simultaneous file access
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger-1_nTrkFilter0.root --TriggerChoice -1 --nTrkFilter 0 &
sleep 1
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger1_nTrkFilter0.root --TriggerChoice 1 --nTrkFilter 0 &
sleep 1

./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger0_nTrkFilter1.root --TriggerChoice 0 --nTrkFilter 1 &
sleep 1
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger-1_nTrkFilter1.root --TriggerChoice -1 --nTrkFilter 1 &
sleep 1
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger1_nTrkFilter1.root --TriggerChoice 1 --nTrkFilter 1 &
sleep 1

./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger0_nTrkFilter-1.root --TriggerChoice 0 --nTrkFilter -1 &
sleep 1
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger-1_nTrkFilter-1.root --TriggerChoice -1 --nTrkFilter -1 &
sleep 1
./ExtractEmptyBXDistributons --Input $INPUTDATAPATH --Output EmptyBX/output_EmptyBX/emptyBXOutput_${DATE}_Trigger1_nTrkFilter-1.root --TriggerChoice 1 --nTrkFilter -1 &
sleep 1
