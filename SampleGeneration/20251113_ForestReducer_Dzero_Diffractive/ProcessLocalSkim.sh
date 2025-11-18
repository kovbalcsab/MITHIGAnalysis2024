#!/bin/bash

FILEPATH=${1}
COUNTER=${2}
OUTPUT=${3}

file="$FILEPATH"

echo "Processing $file"
./Execute --Input "$file" \
    --Output "${OUTPUT}/output_${COUNTER}.root" \
    --Year 2025 \
    --ApplyTriggerRejection 0 \
    --ApplyEventRejection false \
    --ApplyZDCGapRejection false \
    --ApplyDRejection no \
    --ZDCMinus1nThreshold 1000 \
    --ZDCPlus1nThreshold 1000 \
    --IsData true \
    --HideProgressBar true
wait

sleep 0.1
wait
