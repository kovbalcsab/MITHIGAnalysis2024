#!/bin/bash
MAXCORES=40
OUTPUT="output_tmp"
counter=0
filelist=$1
TAG=$2
MERGE18=1
TRIGGERCHOICE=$3
EXECUTABLE=ExecuteHFDistributionExtraction

MERGEDOUTPUT="./output_manual/${TAG}/"
XRDSERV="root://eoscms.cern.ch/"

rm "$MERGEDOUTPUT/MergedOutput_HFDist.root"
rm "$MERGEDOUTPUT/MergedOutput_HFDist_Recalc.root"

# Function to monitor active processes
wait_for_slot() {
  while (($(jobs -r | wc -l) >= MAXCORES)); do
    # Wait a bit before checking again
    sleep 1
  done
}

# Check if the filelist is empty
if [[ ! -s "$filelist" ]]; then
  echo "No matching files found in Samples directory."
  exit 1
fi

echo "File list created successfully: $filelist"
rm -rf $OUTPUT
mkdir -p $OUTPUT
# Loop through each file in the file list
while IFS= read -r file; do
  if [[ -z "$file" ]]; then
    continue
  fi
  echo "Processing $file"
  wait_for_slot
  bash manual/ProcessXRDHFDist.sh $EXECUTABLE $XRDSERV $file $counter $OUTPUT $MERGE18 $TRIGGERCHOICE &
  ((counter++))
done <"$filelist"
wait

mkdir -p "$MERGEDOUTPUT"
hadd "$MERGEDOUTPUT/MergedOutput_HFDist.root" $OUTPUT/output_*.root

./ExecuteRecalculateMeanAndStd --Input "$MERGEDOUTPUT/MergedOutput_HFDist.root" --Output "$MERGEDOUTPUT/MergedOutput_HFDist_Recalc.root"
echo "All done!"
echo "Merged output file: $MERGEDOUTPUT/MergedOutput_HFDist_Recalc.root"
