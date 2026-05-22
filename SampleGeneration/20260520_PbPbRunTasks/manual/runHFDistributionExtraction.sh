#!/bin/bash
MAXCORES=34
OUTPUT="output_tmp"
counter=0
filelist="./eos_first100.txt"
MERGEDOUTPUT="./output/HFDistributions_First_100_mergeSide.root"
MERGEDOUTPUTRECALC="./output/HFDistributions_First_100_mergeSide_recalc.root"
MERGE18=1
rm $MERGEDOUTPUT

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
  echo "Processing $file"
  ./ExecuteHFDistributionExtraction --Input "$file" \
    --OutputRoot "$OUTPUT/output_$counter.root" --do18BinMerging $MERGE18 &
  ((counter++))
  wait_for_slot
done <"$filelist"
wait

hadd $MERGEDOUTPUT $OUTPUT/output_*.root
echo "All done!"
echo "Merged output file: $MERGEDOUTPUT"

./ExecuteRecalculateMeanAndStd --Input $MERGEDOUTPUT --Output $MERGEDOUTPUTRECALC
