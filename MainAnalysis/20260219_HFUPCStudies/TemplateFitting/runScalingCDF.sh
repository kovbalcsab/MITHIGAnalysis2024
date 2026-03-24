#!/bin/bash
set -euo pipefail

# Alternative:
# --FileToFit "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root"
./ScalingCDF --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
             --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
             --VarFitName "HFEMaxPlus_forest" \
             --VarTargetName "HFEMaxPlus_forest" \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingCDFResult_EMC_EBX_plus.root" \
             --NBinsTarget 10000 \
             --NBinsFit 10000 \
             --EMinTarget 1 \
             --EMaxTargetROI 200 \
             --EMinFit 1 \
             --EMaxFitROI 200 \
             --MapFitMin 1 \
             --MapFitMax 6.5 \
             --PolyOrder 2 \
             --PlotYMax 25 \
             --PlotXMax 10 \
             --ForceZeroConstant 1

# Alternative:
# --FileToFit "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root"
./PlotCDFResult --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
             --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
             --VarFitName "HFEMaxPlus_forest" \
             --VarTargetName "HFEMaxPlus_forest" \
             --NBinsTarget 200 \
             --NBinsFit 200 \
             --EMaxTarget 50 \
             --EMaxFit 50 \
             --PlotYMax 1 \
             --PlotXMax 20 \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingCDFResultPlot_EMC_EBX_plus.pdf" \
             --FitResultFileName "TemplateFitting/output_ScalingFit/ScalingCDFResult_EMC_EBX_plus.root" 


# Alternative:
# --FileToFit "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root"
./ScalingCDF --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
             --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
             --VarFitName "HFEMaxMinus_forest" \
             --VarTargetName "HFEMaxMinus_forest" \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingCDFResult_EMC_EBX_minus.root" \
             --NBinsTarget 10000 \
             --NBinsFit 10000 \
             --EMinTarget 1 \
             --EMaxTargetROI 200 \
             --EMinFit 1 \
             --EMaxFitROI 200 \
             --MapFitMin 1 \
             --MapFitMax 6.5 \
             --PolyOrder 2 \
             --PlotYMax 25 \
             --PlotXMax 10 \
             --ForceZeroConstant 1
             

# Alternative:
# --FileToFit "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root"
./PlotCDFResult --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
             --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
             --VarFitName "HFEMaxMinus_forest" \
             --VarTargetName "HFEMaxMinus_forest" \
             --NBinsTarget 200 \
             --NBinsFit 200 \
             --EMaxTarget 50 \
             --EMaxFit 50 \
             --PlotYMax 1 \
             --PlotXMax 20 \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingCDFResultPlot_EMC_EBX_minus.pdf" \
             --FitResultFileName "TemplateFitting/output_ScalingFit/ScalingCDFResult_EMC_EBX_minus.root" 
