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
cp "$SampleSettingCard" "$OutputFileDir/sampleConfig_validation.json"
SampleSettingCard="$OutputFileDir/sampleConfig_validation.json"

jq -c '.OutputFiles[]' "$SampleSettingCard" | while read -r Entry; do
  InputFile=$(echo "$Entry" | jq -r '.InputFile')
  OutputFile=$(echo "$Entry" | jq -r '.OutputFile')
  MeasuredHist=$(echo "$Entry" | jq -r '.MeasuredHist // "hMeasured"')
  RefoldedPrefix=$(echo "$Entry" | jq -r '.RefoldedPrefix // "hRefolded_iter"')
  UnfoldedPrefix=$(echo "$Entry" | jq -r '.UnfoldedPrefix // "hUnfolded_iter"')
  MaxIterations=$(echo "$Entry" | jq -r '.MaxIterations // 10')
  UnfoldedRatioDenominatorIteration=$(echo "$Entry" | jq -r '.UnfoldedRatioDenominatorIteration // 4')
  TrainingQuarter=$(echo "$Entry" | jq -r '.TrainingQuarter // -1')
  ValidationQuarter=$(echo "$Entry" | jq -r '.ValidationQuarter // -1')
  ValidationDataFile=$(echo "$Entry" | jq -r '.ValidationDataFile // empty')
  ValidationTreeName=$(echo "$Entry" | jq -r '.ValidationTreeName // "OutputTree"')
  ValidationVarName=$(echo "$Entry" | jq -r '.ValidationVarName // empty')

  mkdir -p "$(dirname "$OutputFile")"

  cmd=(
    ./TemplateFitting/unfolding_iter_optimization
    --InputFile "$InputFile"
    --OutputFile "$OutputFile"
    --MeasuredHist "$MeasuredHist"
    --RefoldedPrefix "$RefoldedPrefix"
    --UnfoldedPrefix "$UnfoldedPrefix"
    --MaxIterations "$MaxIterations"
    --UnfoldedRatioDenominatorIteration "$UnfoldedRatioDenominatorIteration"
    --TrainingQuarter "$TrainingQuarter"
  )

  if [ "$ValidationQuarter" -ge 0 ]; then
    cmd+=(--ValidationQuarter "$ValidationQuarter")
    if [ -n "$ValidationDataFile" ]; then
      cmd+=(--ValidationDataFile "$ValidationDataFile" --ValidationTreeName "$ValidationTreeName" --ValidationVarName "$ValidationVarName")
    fi
  fi

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done

