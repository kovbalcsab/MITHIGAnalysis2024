#!/usr/bin/env bash
set -euo pipefail

./ScalingFit --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
             --VarFitName "HFEMaxPlus_forest" \
             --VarTargetName "HFEMaxPlus_forest" \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingFitResult_EMC_EBX_plus.root" \
             --aBase 2 \
             --BinsPerGeV 4 \
             --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
             --UseKeysPdf true \
             --FitNumCPU 4 \
             --FitStrategy 2 \
             --FitOffset true \
             --ClipScaledE true \
             --UseAbsJacobian true

./ScalingFit --TargetFile "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root" \
             --VarFitName "HFEMaxMinus_forest" \
             --VarTargetName "HFEMaxMinus_forest" \
             --OutputFileName "TemplateFitting/output_ScalingFit/ScalingFitResult_EMC_EBX_minus.root" \
             --aBase 2 \
             --BinsPerGeV 4 \
             --FileToFit "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root" \
             --UseKeysPdf true \
             --FitNumCPU 4 \
             --FitStrategy 2 \
             --FitOffset true \
             --ClipScaledE true \
             --UseAbsJacobian true
