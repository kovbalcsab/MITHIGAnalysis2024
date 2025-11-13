#!/bin/bash

# Example file list for J/psi search analysis
MAXCORES=20  # too many parallel cores can cause event loss, increase with caution!
FILELIST_PATH="filelist_test.txt"
OUTPUT="/data00/$USER/MakeSampleWideDistributionOutput"
DATE=$(date +%Y%m%d)
NAME="${DATE}_MakeSampleWideDistribution"

rm -rf $OUTPUT &> /dev/null
mkdir -p $OUTPUT

# Function to monitor active processes
wait_for_slot() {
    while (( $(jobs -r | wc -l) >= MAXCORES )); do
        # Wait a bit before checking again
        sleep 1
    done
}

COUNTER=0
while IFS= read -r FILEPATH; do

    OUTPUT_FILE="${OUTPUT}/${NAME}_${COUNTER}.root"
    echo ./makeSampleWideDistributions --inputFile $FILEPATH --outputFile $OUTPUT_FILE &
    ./makeSampleWideDistributions --inputFile $FILEPATH --outputFile $OUTPUT_FILE &

    wait_for_slot
    ((COUNTER++))
done < $FILELIST_PATH
# Wait for all background processes to finish
wait

# Merge output files
rm -f ${OUTPUT}/${NAME}_merged.root &> /dev/null
hadd ${OUTPUT}/${NAME}_merged.root ${OUTPUT}/${NAME}_*.root