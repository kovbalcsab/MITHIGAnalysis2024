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
  ScalingOutputFileName=$(echo "$Entry" | jq -r '.ScalingOutputFileName')
  PlotOutputFileName=$(echo "$Entry" | jq -r '.PlotOutputFileName')

  NBinsTarget=$(echo "$Entry" | jq -r '.NBinsTarget // 10000')
  NBinsFit=$(echo "$Entry" | jq -r '.NBinsFit // 10000')
  EMinTarget=$(echo "$Entry" | jq -r '.EMinTarget // 1')
  EMaxTargetROI=$(echo "$Entry" | jq -r '.EMaxTargetROI // 200')
  EMinFit=$(echo "$Entry" | jq -r '.EMinFit // 1')
  EMaxFitROI=$(echo "$Entry" | jq -r '.EMaxFitROI // 200')
  MapFitMin=$(echo "$Entry" | jq -r '.MapFitMin // 1')
  MapFitMax=$(echo "$Entry" | jq -r '.MapFitMax // 6.5')
  PolyOrder=$(echo "$Entry" | jq -r '.PolyOrder // 2')
  ForceZeroConstant=$(echo "$Entry" | jq -r '.ForceZeroConstant // 1')
  ScalingPlotYMax=$(echo "$Entry" | jq -r '.ScalingPlotYMax // 25')
  ScalingPlotXMax=$(echo "$Entry" | jq -r '.ScalingPlotXMax // 10')

  PlotNBinsTarget=$(echo "$Entry" | jq -r '.PlotNBinsTarget // 200')
  PlotNBinsFit=$(echo "$Entry" | jq -r '.PlotNBinsFit // 200')
  PlotEMaxTarget=$(echo "$Entry" | jq -r '.PlotEMaxTarget // 50')
  PlotEMaxFit=$(echo "$Entry" | jq -r '.PlotEMaxFit // 50')
  PlotYMax=$(echo "$Entry" | jq -r '.PlotYMax // 1')
  PlotXMax=$(echo "$Entry" | jq -r '.PlotXMax // 20')

  mkdir -p "$(dirname "$ScalingOutputFileName")"
  mkdir -p "$(dirname "$PlotOutputFileName")"

  cmdScaling=(
    ./ScalingCDF
    --FileToFit "$FileToFit"
    --TargetFile "$TargetFile"
    --VarFitName "$VarFitName"
    --VarTargetName "$VarTargetName"
    --OutputFileName "$ScalingOutputFileName"
    --NBinsTarget "$NBinsTarget"
    --NBinsFit "$NBinsFit"
    --EMinTarget "$EMinTarget"
    --EMaxTargetROI "$EMaxTargetROI"
    --EMinFit "$EMinFit"
    --EMaxFitROI "$EMaxFitROI"
    --MapFitMin "$MapFitMin"
    --MapFitMax "$MapFitMax"
    --PolyOrder "$PolyOrder"
    --PlotYMax "$ScalingPlotYMax"
    --PlotXMax "$ScalingPlotXMax"
    --ForceZeroConstant "$ForceZeroConstant"
  )

  cmdPlot=(
    ./PlotCDFResult
    --FileToFit "$FileToFit"
    --TargetFile "$TargetFile"
    --VarFitName "$VarFitName"
    --VarTargetName "$VarTargetName"
    --NBinsTarget "$PlotNBinsTarget"
    --NBinsFit "$PlotNBinsFit"
    --EMaxTarget "$PlotEMaxTarget"
    --EMaxFit "$PlotEMaxFit"
    --PlotYMax "$PlotYMax"
    --PlotXMax "$PlotXMax"
    --OutputFileName "$PlotOutputFileName"
    --FitResultFileName "$ScalingOutputFileName"
  )

  echo "Executing >>>>>>"
  echo "${cmdScaling[*]}"
  "${cmdScaling[@]}"

  echo "Executing >>>>>>"
  echo "${cmdPlot[*]}"
  "${cmdPlot[@]}"
done
