#!/bin/bash

DATE=$(date +%Y%m%d)

TriggerChoice=$3
UseZDC=$4
nTrkFilter=0
ZDCM_THRESH=1000
ZDCP_THRESH=1100

MAXCORES=40
XRDSERV="root://eoscms.cern.ch/"
FILELIST=$1
TAG=$2

OUTPUTPATH="output_manual/${TAG}/"
EXECUTABLE=ExtractEmptyBXFullMaps
MERGEOUTPUT="${OUTPUTPATH}/MergedOutput_forestEtaPhi.root"
rm $MERGEOUTPUT &>/dev/null

wait_for_slot() {
  while (($(jobs -r | wc -l) >= MAXCORES)); do
    # Wait a bit before checking again
    sleep 1
  done
}

COUNTER=0
if [[ ! -s "$FILELIST" ]]; then
  echo "File list $FILELIST is empty or does not exist. Please provide a valid file list."
  exit 1
fi

echo "File list: $FILELIST"

COUNTER=1
while IFS= read -r line || [[ -n "$line" ]]; do
  if [[ -z "$line" ]]; then
    continue
  fi

  wait_for_slot

  bash ForestEtaPhiMaps/ProcessXRDPath.sh $EXECUTABLE $XRDSERV $line $COUNTER $OUTPUTPATH $TriggerChoice $nTrkFilter $ZDCM_THRESH $ZDCP_THRESH $UseZDC $MAXCORES &

  ((COUNTER++))
  sleep 1
done <"$FILELIST"

wait
hadd -f $MERGEOUTPUT ${OUTPUTPATH}/output_*.root
rm ${OUTPUTPATH}/output_*.root
wait
echo "Merged output saved to $MERGEOUTPUT"
