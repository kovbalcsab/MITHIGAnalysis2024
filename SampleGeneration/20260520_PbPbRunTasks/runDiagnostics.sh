#!/bin/bash

DATE=$(date +%Y%m%d)
filelist="eos_first10.txt"
TAG="test_$DATE"
DEFFILE="output_manual/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1/MergedOutput_HFDist_Recalc.root"

make

DODETAILS_RUN=1
DOCOARSE_RUN=1
DO_PLOTS=1

rm -rf "output_manual/$TAG/plots"
if [ $DODETAILS_RUN -eq 1 ]; then
  echo "Producing detailed HF distributions..."

  # Produce detailed HF distributions
  bash manual/runHFDistributionExtraction.sh $filelist $TAG
  wait
fi

if [ $DOCOARSE_RUN -eq 1 ]; then
  bash ForestEtaPhiMaps/runForestEtaPhiMaps.sh $filelist $TAG
  wait
fi

if [ $DO_PLOTS -eq 1 ]; then
  mkdir -p "output_manual/$TAG/plots"

  ./ExecutePlotHFDistributionsSides \
    --Input "./output_manual/$TAG/MergedOutput_HFDist_Recalc.root" \
    --OutputDir "./output_manual/$TAG/plots" \
    --OutputPrefix normal \
    --OutputFormat png \
    --IgnoreEmptyWrapHigh true
  #  --HistogramSuffix _diff \

  ./ExecuteDiffHFDistributions --NewInput "./output_manual/$TAG/MergedOutput_HFDist_Recalc.root" \
    --DefaultInput $DEFFILE \
    --Output "./output_manual/$TAG/output_Merged_Recalc_diff_to_23.root"
  wait

  ./ExecutePlotHFDistributionsSides \
    --Input "./output_manual/$TAG/output_Merged_Recalc_diff_to_23.root" \
    --OutputDir "./output_manual/$TAG/plots" \
    --OutputPrefix diff \
    --OutputFormat png \
    --HistogramSuffix _diff \
    --IgnoreEmptyWrapHigh true

  ./PlotForestEtaPhiMaps --Input "output_manual/$TAG/MergedOutput_forestEtaPhi.root" \
    --Input2 "ForestEtaPhiMaps/Output_EmptyBX23Full/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1_20260522/output_Merged.root" \
    --Output "output_manual/$TAG/plots"
fi
