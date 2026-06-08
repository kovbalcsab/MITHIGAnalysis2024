#!/bin/bash

EXECUTABLE=${1}
SERVER=${2}
SERVERPATH=${3}
COUNTER=${4}
OUTPUT=${5}
ADCTHRESHOLD=${6:-19}
TRIGGERCHOICE=${7:-0}

MAXCORES=${8:-20}

mkdir -p "${OUTPUT}/temp_inputs/"
FILEPATH="${OUTPUT}/temp_inputs/job_${COUNTER}.root"
rm $FILEPATH &>/dev/null
xrdcp -C auto -N --parallel $MAXCORES -t 2 $SERVER$SERVERPATH $FILEPATH
wait

echo "Processing $FILEPATH"
./$EXECUTABLE --Input "$FILEPATH" \
  --Output "${OUTPUT}/output_${COUNTER}.root" --ADCThreshold $ADCTHRESHOLD --TriggerChoice $TRIGGERCHOICE
wait

sleep 1
rm $FILEPATH
wait
