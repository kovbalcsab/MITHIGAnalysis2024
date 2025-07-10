#!/bin/bash

# ============================================================
# OO data
# ============================================================
source clean.sh

# Function to monitor active processes
wait_for_slot() {
    while (( $(jobs -r | wc -l) >= MAXCORES )); do
        # Wait a bit before checking again
        sleep 1
    done
}
MAXCORES=20  # too many parallel cores can cause event loss, increase with caution!

INPUTFILEARRAY=(20250705_OO_394153_PhysicsIonPhysics0_074244
                20250705_OO_394153_PhysicsIonPhysics2_250705_074258
                20250705_OO_394153_PhysicsIonPhysics4_25075_072449
                20250706_OO_394153_PhysicsIonPhysics6_250705_120227
                20250706_OO_394153_PhysicsIonPhysics8_250705_120244
                20250706_OO_394153_PhysicsIonPhysics10_250705_074313)

TRKPTMIN=1
TRIGGER=1
OFFLINEHFANDARRAY=(2 6 8 10 12 14 16)
ONLINEHFARRAY=(-1)

INPUTFOLDER=/data00/bakovacs/OOsamples/Skims
for ONLINEHF in "${ONLINEHFARRAY[@]}"; do
  OUTPUTANALYSISFOLDER=HFSystOutput_TRIG${TRIGGER}_TRKPT${TRKPTMIN}_ONLINEHFAND${ONLINEHF}
  mkdir -p ${OUTPUTANALYSISFOLDER}

  #INPUT=/data00/bakovacs/OOsamples/Skims/20250705_OO_394153_FIRSTLOOK.root
  #OUTPUTANALYSIS=HFSystOutput/20250705_OO_394153_FIRSTLOOK_eventSel_min04_HFAND14_HFAND12.root

  for INPUTFILE in "${INPUTFILEARRAY[@]}"; do
    INPUT="${INPUTFOLDER}/${INPUTFILE}.root"

    for OFFLINEHFAND in "${OFFLINEHFANDARRAY[@]}"; do
      OUTPUTANALYSIS="${OUTPUTANALYSISFOLDER}/${INPUTFILE}_TRIG${TRIGGER}_HFOR${OFFLINEHFAND}OFF_HFAND${ONLINEHF}ON_MINPT${TRKPTMIN}.root"
      ./ExecuteChargedHadronRAA \
        --Input $INPUT \
        --Output $OUTPUTANALYSIS \
        --IsData true \
        --IsPP false \
        --UseTrackWeight true \
        --UseEventWeight false \
        --ApplyEventSelection true \
        --MinTrackPt $TRKPTMIN \
        --OnlineHFAND $ONLINEHF \
        --OfflineHFOR $OFFLINEHFAND \
        --ScaleFactor 1 \
        --TriggerChoice $TRIGGER &

      #./ExecuteChargedHadronRAA \
      #  --Input $INPUT \
      #  --Output $OUTPUTANALYSIS \
      #  --IsData true \
      #  --IsPP false \
      #  --UseTrackWeight true \
      #  --UseEventWeight false \
      #  --ApplyEventSelection true \
      #  --MinTrackPt $TRKPTMIN \
      #  --OfflineHFAND $OFFLINEHFAND \
      #  --ScaleFactor 1 \
      #  --TriggerChoice 12 &

        wait_for_slot
    done
  done
  wait

  MERGEDOUTPUTTAG=mergedOutputs/mergedOutput_TRIG${TRIGGER}_TRKPT${TRKPTMIN}_ONLINEHFAND${ONLINEHF}
  for OFFLINEHFAND in "${OFFLINEHFANDARRAY[@]}"; do
    rm ${MERGEDOUTPUTTAG}_HFOR${OFFLINEHFAND}OFF.root
    hadd ${MERGEDOUTPUTTAG}_HFOR${OFFLINEHFAND}OFF.root ${OUTPUTANALYSISFOLDER}/*HFOR${OFFLINEHFAND}OFF*.root
  done
done