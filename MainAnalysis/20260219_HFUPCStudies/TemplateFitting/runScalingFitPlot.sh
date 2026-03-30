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
  OutputFileName=$(echo "$Entry" | jq -r '.OutputFileName // empty')
  ParameterFile=$(echo "$Entry" | jq -r '.ParameterFile // empty')
  PlotFileName=$(echo "$Entry" | jq -r '.PlotFileName // empty')
  XMin=$(echo "$Entry" | jq -r '.XMin // 0')
  XMax=$(echo "$Entry" | jq -r '.XMax // 25')
  NBins=$(echo "$Entry" | jq -r '.NBins // 100')
  TreeName=$(echo "$Entry" | jq -r '.TreeName // "OutputTree"')

  if [ -z "$ParameterFile" ]; then
    if [ -z "$OutputFileName" ]; then
      echo "Each OutputFiles entry must provide ParameterFile or OutputFileName." >&2
      exit 1
    fi
    ParameterFile="${OutputFileName%.root}_fit_parameters.txt"
  fi

  if [ -z "$PlotFileName" ]; then
    if [ -n "$OutputFileName" ]; then
      PlotFileName="${OutputFileName%.root}_plot.pdf"
    else
      echo "Each OutputFiles entry must provide PlotFileName or OutputFileName." >&2
      exit 1
    fi
  fi

  mkdir -p "$(dirname "$PlotFileName")"

  cmd=(
    ./plotScalingResult
    --TargetFile "$TargetFile"
    --FileToFit "$FileToFit"
    --VarFitName "$VarFitName"
    --VarTargetName "$VarTargetName"
    --ParameterFile "$ParameterFile"
    --PlotFileName "$PlotFileName"
    --TreeName "$TreeName"
    --XMin "$XMin"
    --XMax "$XMax"
    --NBins "$NBins"
  )

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done
