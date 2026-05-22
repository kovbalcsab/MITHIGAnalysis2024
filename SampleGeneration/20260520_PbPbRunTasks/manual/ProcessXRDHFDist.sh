#!/bin/bash

EXECUTABLE=${1}
SERVER=${2}
SERVERPATH=${3}
COUNTER=${4}
OUTPUT=${5}
DO18MERGING=${6}

MAXCORES=${7:-20}

mkdir -p "${OUTPUT}/temp_inputs/"
FILEPATH="${OUTPUT}/temp_inputs/job_${COUNTER}.root"
rm $FILEPATH &>/dev/null
xrdcp -C auto -N --parallel $MAXCORES -t 2 $SERVER$SERVERPATH $FILEPATH
wait

echo "Processing $FILEPATH"
./$EXECUTABLE --Input "$FILEPATH" \
  --OutputRoot "${OUTPUT}/output_${COUNTER}.root" --do18BinMerging $DO18MERGING &
wait

sleep 1
rm $FILEPATH
wait
