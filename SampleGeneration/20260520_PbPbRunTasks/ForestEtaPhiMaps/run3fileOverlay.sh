#!/bin/bash
DATE=$(date +%Y%m%d)

DATA26="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404186_20260526/MergedOutput_forestEtaPhi.root"
DATA25="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404359_20260527/MergedOutput_forestEtaPhi.root"
DATA23="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404395_20260528/MergedOutput_forestEtaPhi.root"
#DATA26="output_manual/HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404186_20260526/MergedOutput_forestEtaPhi.root"
#DATA25="output_manual/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1_20260525/MergedOutput_forestEtaPhi.root"
#DATA23="output_manual/HiForest_260218_HIEmptyBX_HIRun2023A_PromptReco_v2_20260526/MergedOutput_forestEtaPhi.root"
TAG="Overlay3BunchNumbers_$DATE"

mkdir -p output_manual/$TAG/plots
./PlotForestEtaPhiMaps --Input $DATA26 \
  --Input2 $DATA25 --Input3 $DATA23 \
  --Output "output_manual/$TAG/plots" --Labels "Spurious_404186,Stable_404359,Stable_404395"
