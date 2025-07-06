#!/bin/bash
source clean.sh

INPUT="/data00/OOsamples/Skims20250704/skim_HiForest_250514_Starlight_DoubleDiffraction_OO_5362GeV_1505PR47944.root,/data00/OOsamples/Skims20250704/skim_HiForest_250514_Starlight_SingleDiffraction_OO_5362GeV_1505PR47944.root,/data00/OOsamples/Skims20250704/skim_HiForest_250520_Hijing_MinimumBias_b015_OO_5362GeV_250518.root"
OUTPUT=/data00/bakovacs/OOsamples/Skims/testShuffle_1M_v2.root

./Shuffle --Input $INPUT \
   --Output $OUTPUT \
   --DoGenLevel false \
   --Year 2025 \
   --IsData false \
   --IsPP false \
   --Fractions 0.06,0.17 \
   --DebugMode true \
   --nEvents 1000000