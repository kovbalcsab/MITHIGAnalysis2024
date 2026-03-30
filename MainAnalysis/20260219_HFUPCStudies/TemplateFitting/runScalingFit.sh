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
  TargetFile=$(echo "$Entry" | jq -r '.TargetFile')
  FileToFit=$(echo "$Entry" | jq -r '.FileToFit')
  VarFitName=$(echo "$Entry" | jq -r '.VarFitName')
  VarTargetName=$(echo "$Entry" | jq -r '.VarTargetName')
  OutputFileName=$(echo "$Entry" | jq -r '.OutputFileName')

  aBase=$(echo "$Entry" | jq -r '.aBase // 1.0')
  bBase=$(echo "$Entry" | jq -r '.bBase // 0.0')
  BinsPerGeV=$(echo "$Entry" | jq -r '.BinsPerGeV // 4')
  UseKeysPdf=$(echo "$Entry" | jq -r '.UseKeysPdf // true')
  FitNumCPU=$(echo "$Entry" | jq -r '.FitNumCPU // 4')
  FitStrategy=$(echo "$Entry" | jq -r '.FitStrategy // 2')
  FitOffset=$(echo "$Entry" | jq -r '.FitOffset // true')
  ClipScaledE=$(echo "$Entry" | jq -r '.ClipScaledE // true')
  UseAbsJacobian=$(echo "$Entry" | jq -r '.UseAbsJacobian // true')
  AMin=$(echo "$Entry" | jq -r '.AMin // 0.0')
  AMax=$(echo "$Entry" | jq -r '.AMax // 5.0')
  BMin=$(echo "$Entry" | jq -r '.BMin // -10.0')
  BMax=$(echo "$Entry" | jq -r '.BMax // 10.0')

  mkdir -p "$(dirname "$OutputFileName")"

  cmd=(
    ./ScalingFit
    --TargetFile "$TargetFile"
    --FileToFit "$FileToFit"
    --VarFitName "$VarFitName"
    --VarTargetName "$VarTargetName"
    --OutputFileName "$OutputFileName"
    --aBase "$aBase"
    --bBase "$bBase"
    --BinsPerGeV "$BinsPerGeV"
    --UseKeysPdf "$UseKeysPdf"
    --FitNumCPU "$FitNumCPU"
    --FitStrategy "$FitStrategy"
    --FitOffset "$FitOffset"
    --ClipScaledE "$ClipScaledE"
    --UseAbsJacobian "$UseAbsJacobian"
    --AMin "$AMin"
    --AMax "$AMax"
    --BMin "$BMin"
    --BMax "$BMax"
  )

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done
