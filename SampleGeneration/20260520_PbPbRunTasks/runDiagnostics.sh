#!/bin/bash

DATE=$(date +%Y%m%d)
filelist="filelist_HIEmptyBX_404471_50k.txt"
TAG="HIEmptyBX_404471_50k_Promptv1_20260608_$DATE"
#TAG="HIExpress_HLT_HIL1NotBptxOR_v15_run404156_20260525"
#DEFDIR="HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404395_20260528"
DEFDIR="HIExpressRawPrime_HLT_HIL1NotBptxOR_v15_run404471_20260530"
#DEFFILE="output_manual/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1_20260525/MergedOutput_HFDist_Recalc.root"
DEFFILE="output_manual/$DEFDIR/MergedOutput_HFDist_Recalc.root"
#DEFFILE_COARSE="output_manual/HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1_20260525/MergedOutput_forestEtaPhi.root"
DEFFILE_COARSE="output_manual/$DEFDIR/MergedOutput_forestEtaPhi.root"
TRIGGERCHOICE_DETAIL=1
TRIGGERCHOICE_COARSE=0
UseZDC=1
BXSEL=0

make

DODETAILS_RUN=0
DOCOARSE_RUN=0
DO_PLOTS=1

rm -rf "output_manual/$TAG/plots"
if [ $DODETAILS_RUN -eq 1 ]; then
  echo "Producing detailed HF distributions..."

  # Produce detailed HF distributions
  bash manual/runHFDistributionExtraction.sh $filelist $TAG $TRIGGERCHOICE_DETAIL
  wait
fi

if [ $DOCOARSE_RUN -eq 1 ]; then
  bash ForestEtaPhiMaps/runForestEtaPhiMaps.sh $filelist $TAG $TRIGGERCHOICE_COARSE $UseZDC $BXSEL
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
  #This is a misslabeling, comparison is to 25 empty BX

  ./ExecutePlotHFDistributionsSides \
    --Input "./output_manual/$TAG/output_Merged_Recalc_diff_to_23.root" \
    --OutputDir "./output_manual/$TAG/plots" \
    --OutputPrefix diff \
    --OutputFormat png \
    --HistogramSuffix _diff \
    --IgnoreEmptyWrapHigh true

  ./PlotForestEtaPhiMaps --Input "output_manual/$TAG/MergedOutput_forestEtaPhi.root,$DEFFILE_COARSE" \
    --Output "output_manual/$TAG/plots"
  ./PlotForestEtaPhiLeading --Input "output_manual/$TAG/MergedOutput_forestEtaPhi.root" \
    --Output "output_manual/$TAG/plots/LeadingCompPlusMinus.pdf"
fi
