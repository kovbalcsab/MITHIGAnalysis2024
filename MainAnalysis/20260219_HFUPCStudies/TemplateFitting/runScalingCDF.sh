#!/bin/bash

./ScalingCDF --FileToFit "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root" \
             \#--FileToFit "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root" \
             --TargetFile "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root" \
             --VarFitName "HFEMaxPlus_forest" \
             --VarTargetName "HFEMaxPlus_forest" \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingCDFResult.root" \
             --NbinsTarget 10000 \
             --NBinsFit 10000 \
             --EMaxTarget 250 \
             --EMaxFit 250 \
             --PlotYMax 25 \
             --PlotXMax 10 

./PlotCDFResult --FileToFit "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root" \
             \#--FileToFit "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root" \
             --TargetFile "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root" \
             --VarFitName "HFEMaxPlus_forest" \
             --VarTargetName "HFEMaxPlus_forest" \
             --NBinsTarget 200 \
             --NBinsFit 200 \
             --EMaxTarget 50 \
             --EMaxFit 50 \
             --PlotYMax 3 \
             --PlotXMax 20 \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingCDFResultPlot.pdf" \
             --FitResultFileName "TemplateFitting/output_ScalingFit/ScalingCDFResult.root" 
