#!/bin/bash
MAXCORES=40
OUTPUT="output_tmp_checkHF2"
counter=0
filelist=$1
TAG=$2
ADCTHRESHOLD=${3:-19}
TRIGGERCHOICE=1
EXECUTABLE=ExecuteCheckHF2LeadingHFDist

MERGEDOUTPUT="./output_manual/${TAG}/"
XRDSERV="root://eoscms.cern.ch/"

rm "$MERGEDOUTPUT/MergedOutput_checkHF2LeadingHFDist.root"

wait_for_slot() {
  while (($(jobs -r | wc -l) >= MAXCORES)); do
    sleep 1
  done
}

if [[ ! -s "$filelist" ]]; then
  echo "No matching files found in Samples directory."
  exit 1
fi

echo "File list created successfully: $filelist"
rm -rf $OUTPUT
mkdir -p $OUTPUT
while IFS= read -r file; do
  if [[ -z "$file" ]]; then
    continue
  fi
  echo "Processing $file"
  wait_for_slot
  bash manual/ProcessXRDCheckHF2.sh $EXECUTABLE $XRDSERV $file $counter $OUTPUT $ADCTHRESHOLD $TRIGGERCHOICE $MAXCORES &
  ((counter++))
done <"$filelist"
wait

mkdir -p "$MERGEDOUTPUT"
hadd "$MERGEDOUTPUT/MergedOutput_checkHF2LeadingHFDist.root" $OUTPUT/output_*.root
mkdir -p "$MERGEDOUTPUT/plots"
./PlotCheckHF2LeadingHFDist --Input "$MERGEDOUTPUT/MergedOutput_checkHF2LeadingHFDist.root" --Output "$MERGEDOUTPUT/plots" --RatioDistanceAroundUnity 0.2
echo "All done!"
echo "Merged output file: $MERGEDOUTPUT/MergedOutput_checkHF2LeadingHFDist.root"
echo "Plots saved to: $MERGEDOUTPUT/plots"
