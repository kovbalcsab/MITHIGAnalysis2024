#!/bin/bash
DATE=$(date +%Y%m%d)

source clean.sh

MAXCORES=20  # too many parallel cores can cause event loss, increase with caution!
NFILES=-1 # number of files to cap the processing at, if -1 processess all files
DOGENLEVEL=0
ISDATA=1
SAMPLETYPE=-1 # 0 for HIJING 00, 1 for Starlight SD, 2 for Starlight DD, 4 for HIJING alpha-O, 5 for Agantyr, -1 for data
SAVETRIGGERBITS=1 # 0 for not HLT saved, 1 for HLT OO, 2 for HLT pO
DEBUGMODE=1
INCLUDEPPSANDFSC=1

NAME="${DATE}_OO_394153_PhysicsIonPhysics8_250705_120244"
# PATHSAMPLE="/eos/cms/store/group/phys_heavyions/jdlang/Run3_OO_2025Data_QuickForest/OO_394153_PhysicsIonPhysics0/crab_OO_394153_PhysicsIonPhysics0/250705_074244/0000"
# set your output directory here
OUTPUT="/data00/bakovacs/OOsamples/Skims/output_$NAME"
MERGEDOUTPUT="/data00/bakovacs/OOsamples/Skims/$NAME.root"
FILELIST="filelists/file_list_OO_394153_PhysicsIonPhysics8_250705_120244.txt"
rm $MERGEDOUTPUT &> /dev/null

# Function to monitor active processes
wait_for_slot() {
    while (( $(jobs -r | wc -l) >= MAXCORES )); do
        # Wait a bit before checking again
        sleep 1
    done
}

# Check if the filelist is empty
if [[ ! -s "$FILELIST" ]]; then
    echo "No matching files found in Samples directory."
    exit 1
fi

echo "File list: $FILELIST"
rm -rf $OUTPUT &> /dev/null
mkdir -p $OUTPUT

# Loop through each file in the file list
COUNTER=0
while IFS= read -r FILEPATH; do

    if [ $NFILES -gt 0 ] && [ $COUNTER -ge $NFILES ]; then
        break
    fi

    echo ./ProcessSingleFile-OO.sh "$FILEPATH" $COUNTER $OUTPUT $DOGENLEVEL $ISDATA $SAMPLETYPE $SAVETRIGGERBITS $DEBUGMODE $INCLUDEPPSANDFSC  &
    ./ProcessSingleFile-OO.sh "$FILEPATH" $COUNTER $OUTPUT $DOGENLEVEL $ISDATA $SAMPLETYPE $SAVETRIGGERBITS $DEBUGMODE $INCLUDEPPSANDFSC &

    wait_for_slot
    ((COUNTER++))
done < "$FILELIST"
wait

hadd $MERGEDOUTPUT $OUTPUT/output_*.root
echo "All done!"
echo "Merged output file: $MERGEDOUTPUT"
