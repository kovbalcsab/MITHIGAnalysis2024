#!/usr/bin/env bash
set -euo pipefail

./plotScalingResult --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
                     --VarFitName "HFEMaxPlus_forest" \
                     --VarTargetName "HFEMaxPlus_forest" \
                     --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
                     --ParameterFile "TemplateFitting/output_ScalingFit/ScalingFitResult_EMC_EBX_plus_fit_parameters.txt" \
                     --PlotFileName "TemplateFitting/output_ScalingFit/ScalingFitResultPlot_EMC_EBX_plus.pdf" \
                     --XMin 0 \
                     --XMax 25 \
                     --NBins 100

./plotScalingResult --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
                     --VarFitName "HFEMaxMinus_forest" \
                     --VarTargetName "HFEMaxMinus_forest" \
                     --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
                     --ParameterFile "TemplateFitting/output_ScalingFit/ScalingFitResult_EMC_EBX_minus_fit_parameters.txt" \
                     --PlotFileName "TemplateFitting/output_ScalingFit/ScalingFitResultPlot_EMC_EBX_minus.pdf" \
                     --XMin 0 \
                     --XMax 25 \
                     --NBins 100

                     
