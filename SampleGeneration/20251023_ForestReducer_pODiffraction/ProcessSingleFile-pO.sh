#!/bin/bash

FILEPATH=${1}
COUNTER=${2}
OUTPUT=${3}
ISDATA=${4}
APPLYTRIGGERREJECTION=${5}
APPLYEVENTREJECTION=${6}
SAMPLETYPE=${7}

./Execute --Input "$FILEPATH" \
   --Output ${OUTPUT}/output_${COUNTER}.root \
   --IsData $ISDATA \
   --CollisionSystem pO \
   --Fraction 1.0 \
   --ApplyTriggerRejection $APPLYTRIGGERREJECTION \
   --ApplyEventRejection $APPLYEVENTREJECTION \
   --sampleType $SAMPLETYPE \
   --DebugMode false \
   --HideProgressBar true
wait

sleep 0.2
wait
