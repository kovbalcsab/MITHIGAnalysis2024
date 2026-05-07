#!/usr/bin/env bash
set -euo pipefail

if [ "${1:-}" = "" ]; then
  echo "Usage: $0 <SampleSettingCard.json>" >&2
  exit 1
fi

SampleSettingCard="$1"
if [ ! -f "$SampleSettingCard" ]; then
  echo "Config card not found: $SampleSettingCard" >&2
  exit 1
fi

OutputFileDir=$(jq -r '.OutputFileDir' "$SampleSettingCard")
mkdir -p "$OutputFileDir"
cp "$SampleSettingCard" "$OutputFileDir/sampleConfig.json"
SampleSettingCard="$OutputFileDir/sampleConfig.json"

jq -c '.OutputFiles[]' "$SampleSettingCard" | while read -r Entry; do
  DataFile=$(echo "$Entry" | jq -r '.DataFile')
  NoiseFile=$(echo "$Entry" | jq -r '.NoiseFile')
  VarDataName=$(echo "$Entry" | jq -r '.VarDataName')
  VarNoiseName=$(echo "$Entry" | jq -r '.VarNoiseName // .VarDataName')
  OutputFileName=$(echo "$Entry" | jq -r '.OutputFileName')

  XMin=$(echo "$Entry" | jq -r '.XMin // 0')
  XMax=$(echo "$Entry" | jq -r '.XMax // 200')
  BinsPerGeV=$(echo "$Entry" | jq -r '.BinsPerGeV // 4')
  Iterations=$(echo "$Entry" | jq -r '.Iterations // 10')
  DataQuarter=$(echo "$Entry" | jq -r '.DataQuarter // 0')
  LegacyBinEdges=$(echo "$Entry" | jq -r 'if (.BinEdges // null) == null then "" elif (.BinEdges | type) == "array" then (.BinEdges | map(tostring) | join(",")) else (.BinEdges | tostring) end')
  MeasuredBinEdges=$(echo "$Entry" | jq -r 'if (.MeasuredBinEdges // null) == null then "" elif (.MeasuredBinEdges | type) == "array" then (.MeasuredBinEdges | map(tostring) | join(",")) else (.MeasuredBinEdges | tostring) end')
  UnfoldedBinEdges=$(echo "$Entry" | jq -r 'if (.UnfoldedBinEdges // null) == null then "" elif (.UnfoldedBinEdges | type) == "array" then (.UnfoldedBinEdges | map(tostring) | join(",")) else (.UnfoldedBinEdges | tostring) end')

  mkdir -p "$(dirname "$OutputFileName")"

  cmd=(
    ./TemplateFitting/Unfold
    --DataFile "$DataFile"
    --NoiseFile "$NoiseFile"
    --VarDataName "$VarDataName"
    --VarNoiseName "$VarNoiseName"
    --XMin "$XMin"
    --XMax "$XMax"
    --BinsPerGeV "$BinsPerGeV"
    --Iterations "$Iterations"
    --DataQuarter "$DataQuarter"
    --OutputFileName "$OutputFileName"
  )

  if [ -n "$LegacyBinEdges" ]; then
    cmd+=(--BinEdges "$LegacyBinEdges")
  fi
  if [ -n "$MeasuredBinEdges" ]; then
    cmd+=(--MeasuredBinEdges "$MeasuredBinEdges")
  fi
  if [ -n "$UnfoldedBinEdges" ]; then
    cmd+=(--UnfoldedBinEdges "$UnfoldedBinEdges")
  fi

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done
