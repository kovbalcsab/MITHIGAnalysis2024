#!/bin/bash

./ExecutePlotHFDistributionsSides \
  --Input ./output/HFDistributions_First_100_mergeSide_recalc.root \
  --OutputDir output/pdf_plots \
  --OutputPrefix HFDistributions_First_100_mergeSide_recalc \
  --OutputFormat png \
  --HistogramSuffix _diff \
  --IgnoreEmptyWrapHigh true
