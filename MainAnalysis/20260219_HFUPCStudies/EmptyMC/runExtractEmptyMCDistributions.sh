#!/bin/bash

INPUTDATAPATH=/eos/cms/store/group/phys_heavyions/wangj/Forest2025pp/Dzero_260311-ydiffmva_RelValNuGun_151X_mcRun3_2025_realistic_v4_STD_RegeneratedGS_2025_noPU-v2.root

DATE=$(date +"%Y%m%d")

# Create output folder if it doesn't exist
mkdir -p EmptyMC/output_EmptyMC

# This script runs the ExtractEmptyMCDistributons.C macro to extract the empty Monte Carlo distributions from the specified input file.
./ExtractEmptyMCDistributons --Input $INPUTDATAPATH --Output EmptyMC/output_EmptyMC/emptyMCOutput_${DATE}.root &