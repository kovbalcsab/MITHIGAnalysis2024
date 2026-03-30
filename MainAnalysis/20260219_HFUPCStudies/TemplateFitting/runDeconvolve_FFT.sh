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
  KernelFile=$(echo "$Entry" | jq -r '.KernelFile')
  VarDataName=$(echo "$Entry" | jq -r '.VarDataName')
  VarKernelName=$(echo "$Entry" | jq -r '.VarKernelName // .VarDataName')
  OutputFileName=$(echo "$Entry" | jq -r '.OutputFileName')

  XMin=$(echo "$Entry" | jq -r '.XMin // 0')
  XMax=$(echo "$Entry" | jq -r '.XMax // 200')
  BinsPerGeV=$(echo "$Entry" | jq -r '.BinsPerGeV // 4')
  Iterations=$(echo "$Entry" | jq -r '.Iterations // 6')
  DampingMin=$(echo "$Entry" | jq -r '.DampingMin // 0.001')
  DampingMax=$(echo "$Entry" | jq -r '.DampingMax // 1')
  LogDamping=$(echo "$Entry" | jq -r '.LogDamping // 1')
  DataQuarter=$(echo "$Entry" | jq -r '.DataQuarter // 0')

  mkdir -p "$(dirname "$OutputFileName")"

  cmd=(
    ./TemplateFitting/Deconvolve_FFT
    --DataFile "$DataFile"
    --KernelFile "$KernelFile"
    --VarDataName "$VarDataName"
    --VarKernelName "$VarKernelName"
    --XMin "$XMin"
    --XMax "$XMax"
    --BinsPerGeV "$BinsPerGeV"
    --Iterations "$Iterations"
    --DampingMin "$DampingMin"
    --DampingMax "$DampingMax"
    --LogDamping "$LogDamping"
    --DataQuarter "$DataQuarter"
    --OutputFileName "$OutputFileName"
  )

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done
