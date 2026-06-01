#!/bin/bash

EXECUTABLE=${1}
SERVER=${2}
SERVERPATH=${3}
COUNTER=${4}
OUTPUT=${5}
TriggerChoice=${6:-0}
nTrkFilter=${7:-1}

ZDCM_THRESH=${8:-1000}
ZDCP_THRESH=${9:-1100}
UseZDC=${10:-1}
BXSEL=${11:-0}
MAXCORES=${12:-20}

mkdir -p "${OUTPUT}/temp_inputs/"
FILEPATH="${OUTPUT}/temp_inputs/job_${COUNTER}.root"
rm $FILEPATH &>/dev/null
xrdcp -C auto -N --parallel $MAXCORES -t 2 $SERVER$SERVERPATH $FILEPATH
wait

echo "Processing $FILEPATH"
./$EXECUTABLE --Input "$FILEPATH" \
  --Output "${OUTPUT}/output_${COUNTER}.root" \
  --ZDCMinus1nThreshold $ZDCM_THRESH \
  --ZDCPlus1nThreshold $ZDCP_THRESH \
  --TriggerChoice $TriggerChoice \
  --nTrkFilter $nTrkFilter \
  --UseZDC $UseZDC \
  --BXSel $BXSEL
wait

sleep 1
rm $FILEPATH
wait
