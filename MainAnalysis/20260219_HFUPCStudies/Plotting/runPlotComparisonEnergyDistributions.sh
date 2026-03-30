#!/usr/bin/env bash
set -euo pipefail

if [ "${1:-}" = "" ]; then
  echo "Usage: $0 <PlotConfig.json>" >&2
  exit 1
fi

SampleSettingCard="$1"
if [ ! -f "$SampleSettingCard" ]; then
  echo "Config card not found: $SampleSettingCard" >&2
  exit 1
fi

plottingDir=$(jq -r '.plottingDir // empty' "$SampleSettingCard")
if [ -z "$plottingDir" ]; then
  echo "Config must contain plottingDir." >&2
  exit 1
fi

mkdir -p "$plottingDir"
cp "$SampleSettingCard" "$plottingDir/plottingEnergyDistConfig.json"
SampleSettingCard="$plottingDir/plottingEnergyDistConfig.json"

jq -c '.Plots[]' "$SampleSettingCard" | while read -r Plot; do
  PlotBaseName=$(echo "$Plot" | jq -r '.PlotBaseName // "EnergyComparison.pdf"')
  StatMatrixBaseName=$(echo "$Plot" | jq -r '.StatMatrixBaseName // empty')
  InputFileNames=$(echo "$Plot" | jq -r '.InputFileNames')
  VarNames=$(echo "$Plot" | jq -r '.VarNames')
  IsTreeVar=$(echo "$Plot" | jq -r '.IsTreeVar')
  Labels=$(echo "$Plot" | jq -r '.Labels')

  MinDzeroPT=$(echo "$Plot" | jq -r '.MinDzeroPT // empty')
  MaxDzeroPT=$(echo "$Plot" | jq -r '.MaxDzeroPT // empty')
  MinDzeroY=$(echo "$Plot" | jq -r '.MinDzeroY // empty')
  MaxDzeroY=$(echo "$Plot" | jq -r '.MaxDzeroY // empty')
  IsGammaN=$(echo "$Plot" | jq -r '.IsGammaN // empty')

  xTitle=$(echo "$Plot" | jq -r '.xTitle // empty')
  yTitle=$(echo "$Plot" | jq -r '.yTitle // empty')
  plotTitle=$(echo "$Plot" | jq -r '.plotTitle // empty')
  logY=$(echo "$Plot" | jq -r '.logY // empty')
  xMin=$(echo "$Plot" | jq -r '.xMin // empty')
  xMax=$(echo "$Plot" | jq -r '.xMax // empty')
  yMin=$(echo "$Plot" | jq -r '.yMin // empty')
  yMax=$(echo "$Plot" | jq -r '.yMax // empty')
  normalizeToUnity=$(echo "$Plot" | jq -r '.normalizeToUnity // empty')
  rebinFactor=$(echo "$Plot" | jq -r '.rebinFactor // empty')
  nBins=$(echo "$Plot" | jq -r '.NBins // empty')
  doRatio=$(echo "$Plot" | jq -r '.doRatio // empty')

  OutputDir="$plottingDir"
  if [ -n "$MinDzeroPT" ] && [ -n "$MaxDzeroPT" ] && [ -n "$MinDzeroY" ] && [ -n "$MaxDzeroY" ] && [ -n "$IsGammaN" ]; then
    OutputDir="$plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}"
  fi
  mkdir -p "$OutputDir"

  OutputFileName="$OutputDir/$PlotBaseName"
  if [ -z "$StatMatrixBaseName" ]; then
    StatMatrixOutputFileName="${OutputFileName%.pdf}_Chi2PValueMatrix.pdf"
  else
    StatMatrixOutputFileName="$OutputDir/$StatMatrixBaseName"
  fi
  WassersteinBaseName=$(echo "$Plot" | jq -r '.WassersteinMatrixBaseName // empty')
  JSDistanceBaseName=$(echo "$Plot" | jq -r '.JSDistanceMatrixBaseName // empty')
  HellingerBaseName=$(echo "$Plot" | jq -r '.HellingerMatrixBaseName // empty')
  PoissonDevianceBaseName=$(echo "$Plot" | jq -r '.PoissonDevianceMatrixBaseName // empty')

  if [ -z "$WassersteinBaseName" ]; then
    WassersteinOutputFileName="${OutputFileName%.pdf}_WassersteinMatrix.pdf"
  else
    WassersteinOutputFileName="$OutputDir/$WassersteinBaseName"
  fi
  if [ -z "$JSDistanceBaseName" ]; then
    JSDistanceOutputFileName="${OutputFileName%.pdf}_JSDistanceMatrix.pdf"
  else
    JSDistanceOutputFileName="$OutputDir/$JSDistanceBaseName"
  fi
  if [ -z "$HellingerBaseName" ]; then
    HellingerOutputFileName="${OutputFileName%.pdf}_HellingerMatrix.pdf"
  else
    HellingerOutputFileName="$OutputDir/$HellingerBaseName"
  fi
  if [ -z "$PoissonDevianceBaseName" ]; then
    PoissonDevianceOutputFileName="${OutputFileName%.pdf}_PoissonDevianceMatrix.pdf"
  else
    PoissonDevianceOutputFileName="$OutputDir/$PoissonDevianceBaseName"
  fi

  StatRebin=$(echo "$Plot" | jq -r '.StatRebin // 1')
  StatEpsilon=$(echo "$Plot" | jq -r '.StatEpsilon // 1e-12')
  PoissonFloorMu=$(echo "$Plot" | jq -r '.PoissonFloorMu // 1e-9')

  cmd=(
    ./PlotComparisonEnergyDistributions
    --InputFileNames "$InputFileNames"
    --OutputFileName "$OutputFileName"
    --VarNames "$VarNames"
    --IsTreeVar "$IsTreeVar"
    --Labels "$Labels"
    --StatMatrixOutputFileName "$StatMatrixOutputFileName"
    --WassersteinOutputFileName "$WassersteinOutputFileName"
    --JSDistanceOutputFileName "$JSDistanceOutputFileName"
    --HellingerOutputFileName "$HellingerOutputFileName"
    --PoissonDevianceOutputFileName "$PoissonDevianceOutputFileName"
    --StatRebin "$StatRebin"
    --StatEpsilon "$StatEpsilon"
    --PoissonFloorMu "$PoissonFloorMu"
  )

  [ -n "$MinDzeroPT" ] && cmd+=(--MinDzeroPT "$MinDzeroPT")
  [ -n "$MaxDzeroPT" ] && cmd+=(--MaxDzeroPT "$MaxDzeroPT")
  [ -n "$MinDzeroY" ] && cmd+=(--MinDzeroY "$MinDzeroY")
  [ -n "$MaxDzeroY" ] && cmd+=(--MaxDzeroY "$MaxDzeroY")
  [ -n "$IsGammaN" ] && cmd+=(--IsGammaN "$IsGammaN")
  [ -n "$xTitle" ] && cmd+=(--XTitle "$xTitle")
  [ -n "$yTitle" ] && cmd+=(--YTitle "$yTitle")
  [ -n "$plotTitle" ] && cmd+=(--PlotTitle "$plotTitle")
  [ -n "$logY" ] && cmd+=(--LogY "$logY")
  [ -n "$xMin" ] && cmd+=(--XMin "$xMin")
  [ -n "$xMax" ] && cmd+=(--XMax "$xMax")
  [ -n "$yMin" ] && cmd+=(--YMin "$yMin")
  [ -n "$yMax" ] && cmd+=(--YMax "$yMax")
  [ -n "$normalizeToUnity" ] && cmd+=(--Normalize "$normalizeToUnity")
  [ -n "$rebinFactor" ] && cmd+=(--Rebin "$rebinFactor")
  [ -n "$nBins" ] && cmd+=(--NBins "$nBins")
  [ -n "$doRatio" ] && cmd+=(--DoRatio "$doRatio")

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done
