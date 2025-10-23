#!/bin/bash

FILEPATH=${1}
COUNTER=${2}
OUTPUT=${3}
ISDATA=${4}
APPLYTRIGGERREJECTION=${5}
APPLYEVENTREJECTION=${6}
SAMPLETYPE=${7}
SERVER=${8}
MAXCORES=${9}

mkdir -p "${OUTPUT}/temp_inputs/"
FILE="${OUTPUT}/temp_inputs/job_${COUNTER}.root"
rm $FILE &> /dev/null
xrdcp -N --parallel $MAXCORES -t 2 $SERVER$FILEPATH $FILE
wait

./Execute --Input "$FILE" \
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
rm $FILE
wait
