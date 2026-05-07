#!/bin/bash

DATE=$(date +%Y%m%d)

TriggerChoice=-2
nTrkFilter=0
ZDCM_THRESH=1000
ZDCP_THRESH=1100

MAXCORES=40
XRDSERV="root://eoscms.cern.ch/"
TAG="HiForest_260218_RelValNuGun_151X_mcRun3_2025_realistic_v4_STD_RegeneratedGS_2025_noPU-v2"

OUTPUTPATH="ForestEtaPhiMaps/Output_emptyMC/${TAG}_${DATE}"
EXECUTABLE=ExtractEmptyBXFullMaps
FILELIST="ForestEtaPhiMaps/filelist_emptyMC.txt"
MERGEOUTPUT="${OUTPUTPATH}/MergedOutput.root"
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
rm -rf $OUTPUTPATH &>/dev/null
mkdir -p $OUTPUTPATH

COUNTER=1
while IFS= read -r line || [[ -n "$line" ]]; do
  if [[ -z "$line" ]]; then
    continue
  fi

  wait_for_slot

  bash ForestEtaPhiMaps/ProcessXRDPath.sh $EXECUTABLE $XRDSERV $line $COUNTER $OUTPUTPATH $TriggerChoice $nTrkFilter $ZDCM_THRESH $ZDCP_THRESH $MAXCORES &

  ((COUNTER++))
  sleep 1
done <"$FILELIST"
