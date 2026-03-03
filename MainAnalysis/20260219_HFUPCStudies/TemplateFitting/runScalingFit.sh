#!/bin/bash

./ScalingFit --FileToFit "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root" \
             --TargetFile "EmptyZB/output_EmptyZB/emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root" \
             --VarFitName "HFEMaxPlus_forest" \
             --VarTargetName "HFEMaxPlus_forest" \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingFitResult.root" \
             --aBase 2.2 \
             --BinsPerGeV 2
