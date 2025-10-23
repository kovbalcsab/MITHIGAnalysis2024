#!/bin/bash
DATE=$(date +%Y%m%d)

source clean.sh

MAXCORES=40  # too many parallel cores can cause event loss, increase with caution!
NFILES=-1 # number of files to cap the processing at, if -1 processess all files
ISDATA=1
APPLYTRIGGERREJECTION=0 #  trigger = 0 for no rejection, 1 for ZeroBias, 2 for MinBias
APPLYEVENTREJECTION=0
SAMPLETYPE=-1 # -1 for data
INPUT_ON_XRD=1 # 0 for local files, 1 for xrd files
XRDSERV="root://eoscms.cern.ch/" # eos xrootd server, path should start /store/group...

# ============================================================
# OO data, low pT PD, Vipul's filelist for crosscheck
# ============================================================
NAME="${DATE}_Skim_test_pO_sample"
FILELIST="filelists/test_filelist.txt"

# set your output directory here
OUTPUT="/data00/$USER/pOsamples/Skims/output_$NAME/"



# Function to monitor active processes
wait_for_slot() {
    while (( $(jobs -r | wc -l) >= MAXCORES )); do
        # Wait a bit before checking again
        sleep 1
    done
}

rm -rf $OUTPUT &> /dev/null
mkdir -p $OUTPUT

# Loop through each file in the file list
COUNTER=0
while IFS= read -r FILEPATH; do

    if [ $NFILES -gt 0 ] && [ $COUNTER -ge $NFILES ]; then
        break
    fi

    if (( $INPUT_ON_XRD == 1 )); then
        echo ./ProcessSingleFile-pO-xrd.sh $FILEPATH $COUNTER $OUTPUT $ISDATA $APPLYTRIGGERREJECTION $APPLYEVENTREJECTION $SAMPLETYPE $XRDSERV $MAXCORES &
        ./ProcessSingleFile-pO-xrd.sh $FILEPATH $COUNTER $OUTPUT $ISDATA $APPLYTRIGGERREJECTION $APPLYEVENTREJECTION $SAMPLETYPE $XRDSERV $MAXCORES &
    else
        echo ./ProcessSingleFile-pO.sh $FILEPATH $COUNTER $OUTPUT $ISDATA $APPLYTRIGGERREJECTION $APPLYEVENTREJECTION $SAMPLETYPE &
        ./ProcessSingleFile-pO.sh $FILEPATH $COUNTER $OUTPUT $ISDATA $APPLYTRIGGERREJECTION $APPLYEVENTREJECTION $SAMPLETYPE &
    fi

    wait_for_slot
    ((COUNTER++))
done < $FILELIST
wait

echo "Processing COMPLETE"
