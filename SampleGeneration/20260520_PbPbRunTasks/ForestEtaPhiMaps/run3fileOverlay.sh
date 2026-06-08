#!/bin/bash
DATE=$(date +%Y%m%d)

DATA1="output_manual/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1_reforest_Lumimasked_20260604/MergedOutput_forestEtaPhi.root"
DATA2="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404471_20260530/MergedOutput_forestEtaPhi.root"
DATA3="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404471_2025hfcalib_20260530/MergedOutput_forestEtaPhi.root"
DATA4="output_manual/HiForest_260218_HIEmptyBX_HIRun2023A_PromptReco_v2_20260526/MergedOutput_forestEtaPhi.root"
#DATA26="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404186_20260526/MergedOutput_forestEtaPhi.root"
#DATA25="output_manual/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1_20260525/MergedOutput_forestEtaPhi.root"
#DATA23="output_manual/HiForest_260218_HIEmptyBX_HIRun2023A_PromptReco_v2_20260526/MergedOutput_forestEtaPhi.root"
DORATIO=1
TAG="Overlay4_CalibComp_doRatio${DORATIO}_${DATE}"

mkdir -p output_manual/$TAG/plots
./PlotForestEtaPhiMaps --Input "${DATA1},${DATA2},${DATA3},${DATA4}" \
  --Output "output_manual/$TAG/plots" --Labels "2025_(2025_calib.),2026_(2026_calib.),2026_(2025_calib.),2023_(2023_calib.)" --DoRatio $DORATIO --DoSelfNormalize 1 --DoBinomialRatio 0 --RatioDistanceAroundUnity 0.75
cp ForestEtaPhiMaps/run3fileOverlay.sh output_manual/$TAG/plots/remake.sh
