#!/bin/bash

./ExecutePlotHFDistributionsSides \
  --Input ./output_manual/EmptyBXFull23_HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1/output_Merged_Recalc.root \
  --OutputDir output_manual/EmptyBXFull23_HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1/pdf_plots \
  --OutputPrefix EmptyBxFull23 \
  --OutputFormat png \
  --IgnoreEmptyWrapHigh true
#  --HistogramSuffix _diff \
