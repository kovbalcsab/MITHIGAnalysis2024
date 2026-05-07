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

OutputFileDir=$(jq -r '.OutputFileDir // "TemplateFitting/output_unfolding"' "$SampleSettingCard")
mkdir -p "$OutputFileDir"
cp "$SampleSettingCard" "$OutputFileDir/sampleConfig_validation_method_comparison.json"
SampleSettingCard="$OutputFileDir/sampleConfig_validation_method_comparison.json"

jq -c '.OutputFiles[]' "$SampleSettingCard" | while read -r Entry; do
  InputFileNames=$(echo "$Entry" | jq -r '.InputFileNames // empty')
  SelectedIterations=$(echo "$Entry" | jq -r '.SelectedIterations // empty')
  Labels=$(echo "$Entry" | jq -r '.Labels // empty')
  Colors=$(echo "$Entry" | jq -r '.Colors // empty')

  if [ -z "$InputFileNames" ] || [ -z "$SelectedIterations" ] || [ -z "$Labels" ]; then
    echo "InputFileNames, SelectedIterations, and Labels are required in each OutputFiles[] entry." >&2
    exit 1
  fi

  OutputPrefix=$(echo "$Entry" | jq -r '.OutputPrefix // "TemplateFitting/output_unfolding/validation_method_comparison_selected"')
  OutputUnfoldedFileName=$(echo "$Entry" | jq -r '.OutputUnfoldedFileName // empty')
  OutputClosureRatioFileName=$(echo "$Entry" | jq -r '.OutputClosureRatioFileName // empty')
  OutputClosureAbsoluteFileName=$(echo "$Entry" | jq -r '.OutputClosureAbsoluteFileName // empty')
  OutputClosureAbsoluteZoomFileName=$(echo "$Entry" | jq -r '.OutputClosureAbsoluteZoomFileName // empty')
  OutputChi2FileName=$(echo "$Entry" | jq -r '.OutputChi2FileName // empty')

  if [ -z "$OutputUnfoldedFileName" ]; then OutputUnfoldedFileName="${OutputPrefix}_unfolded.pdf"; fi
  if [ -z "$OutputClosureRatioFileName" ]; then OutputClosureRatioFileName="${OutputPrefix}_closure_ratio.pdf"; fi
  if [ -z "$OutputClosureAbsoluteFileName" ]; then OutputClosureAbsoluteFileName="${OutputPrefix}_closure_absolute.pdf"; fi
  if [ -z "$OutputClosureAbsoluteZoomFileName" ]; then OutputClosureAbsoluteZoomFileName="${OutputPrefix}_closure_absolute_x0_25.pdf"; fi
  if [ -z "$OutputChi2FileName" ]; then OutputChi2FileName="${OutputPrefix}_chi2ndf.pdf"; fi

  UnfoldedHistPrefix=$(echo "$Entry" | jq -r '.UnfoldedHistPrefix // "hUnfoldedAbs_iter"')
  ClosureHistPrefix=$(echo "$Entry" | jq -r '.ClosureHistPrefix // "hRefoldedToMeasuredRatio_iter"')
  RefoldedAbsHistPrefix=$(echo "$Entry" | jq -r '.RefoldedAbsHistPrefix // "hRefoldedAbs_iter"')
  MeasuredHistName=$(echo "$Entry" | jq -r '.MeasuredHistName // "hMeasuredForOverlay"')
  Chi2HistName=$(echo "$Entry" | jq -r '.Chi2HistName // "hChi2NDFVsIter"')

  UnfoldedXTitle=$(echo "$Entry" | jq -r '.UnfoldedXTitle // "HF~Energy~[GeV]"')
  UnfoldedYTitle=$(echo "$Entry" | jq -r '.UnfoldedYTitle // "Entries"')
  ClosureXTitle=$(echo "$Entry" | jq -r '.ClosureXTitle // "HF~Energy~[GeV]"')
  ClosureYTitle=$(echo "$Entry" | jq -r '.ClosureYTitle // "Refolded/Measured"')
  ClosureAbsXTitle=$(echo "$Entry" | jq -r '.ClosureAbsXTitle // "HF~Energy~[GeV]"')
  ClosureAbsYTitle=$(echo "$Entry" | jq -r '.ClosureAbsYTitle // "1/N~dN/dE"')
  Chi2XTitle=$(echo "$Entry" | jq -r '.Chi2XTitle // "Regularization~iteration"')
  Chi2YTitle=$(echo "$Entry" | jq -r '.Chi2YTitle // "#chi^{2}/NDF"')
  UnfoldedTitle=$(echo "$Entry" | jq -r '.UnfoldedTitle // "Selected~iteration~unfolded~comparison"')
  ClosureTitle=$(echo "$Entry" | jq -r '.ClosureTitle // "Selected~iteration~closure~ratio~comparison"')
  ClosureAbsTitle=$(echo "$Entry" | jq -r '.ClosureAbsTitle // "Selected~iteration~measured~vs~refolded~comparison"')
  Chi2Title=$(echo "$Entry" | jq -r '.Chi2Title // "Method~comparison~of~#chi^{2}/NDF~scan"')

  NormalizeAbsolute=$(echo "$Entry" | jq -r '.NormalizeAbsolute // true')
  ZoomXMin=$(echo "$Entry" | jq -r '.ZoomXMin // 0')
  ZoomXMax=$(echo "$Entry" | jq -r '.ZoomXMax // 25')
  Chi2RemoveErrorBars=$(echo "$Entry" | jq -r '.Chi2RemoveErrorBars // true')
  DoRatioToFirst=$(echo "$Entry" | jq -r '.DoRatioToFirst // true')
  RatioYMin=$(echo "$Entry" | jq -r '.RatioYMin // 0.0')
  RatioYMax=$(echo "$Entry" | jq -r '.RatioYMax // 2.0')
  LogYUnfolded=$(echo "$Entry" | jq -r '.LogYUnfolded // true')
  RequireAllMethods=$(echo "$Entry" | jq -r '.RequireAllMethods // true')
  XMin=$(echo "$Entry" | jq -r '.XMin // 0')
  XMax=$(echo "$Entry" | jq -r '.XMax // 0')
  YMinUnfolded=$(echo "$Entry" | jq -r '.YMinUnfolded // 0')
  YMaxUnfolded=$(echo "$Entry" | jq -r '.YMaxUnfolded // 0')
  YMinClosure=$(echo "$Entry" | jq -r '.YMinClosure // 0')
  YMaxClosure=$(echo "$Entry" | jq -r '.YMaxClosure // 0')
  YMinClosureAbs=$(echo "$Entry" | jq -r '.YMinClosureAbs // 0')
  YMaxClosureAbs=$(echo "$Entry" | jq -r '.YMaxClosureAbs // 0')
  YMinChi2=$(echo "$Entry" | jq -r '.YMinChi2 // 0')
  YMaxChi2=$(echo "$Entry" | jq -r '.YMaxChi2 // 0')

  mkdir -p "$(dirname "$OutputPrefix")"
  mkdir -p "$(dirname "$OutputUnfoldedFileName")"
  mkdir -p "$(dirname "$OutputClosureRatioFileName")"
  mkdir -p "$(dirname "$OutputClosureAbsoluteFileName")"
  mkdir -p "$(dirname "$OutputClosureAbsoluteZoomFileName")"
  mkdir -p "$(dirname "$OutputChi2FileName")"

  cmd=(
    ./TemplateFitting/plot_validation_method_comparison
    --InputFileNames "$InputFileNames"
    --SelectedIterations "$SelectedIterations"
    --Labels "$Labels"
    --OutputPrefix "$OutputPrefix"
    --OutputUnfoldedFileName "$OutputUnfoldedFileName"
    --OutputClosureRatioFileName "$OutputClosureRatioFileName"
    --OutputClosureAbsoluteFileName "$OutputClosureAbsoluteFileName"
    --OutputClosureAbsoluteZoomFileName "$OutputClosureAbsoluteZoomFileName"
    --OutputChi2FileName "$OutputChi2FileName"
    --UnfoldedHistPrefix "$UnfoldedHistPrefix"
    --ClosureHistPrefix "$ClosureHistPrefix"
    --RefoldedAbsHistPrefix "$RefoldedAbsHistPrefix"
    --MeasuredHistName "$MeasuredHistName"
    --Chi2HistName "$Chi2HistName"
    --UnfoldedXTitle "$UnfoldedXTitle"
    --UnfoldedYTitle "$UnfoldedYTitle"
    --ClosureXTitle "$ClosureXTitle"
    --ClosureYTitle "$ClosureYTitle"
    --ClosureAbsXTitle "$ClosureAbsXTitle"
    --ClosureAbsYTitle "$ClosureAbsYTitle"
    --Chi2XTitle "$Chi2XTitle"
    --Chi2YTitle "$Chi2YTitle"
    --UnfoldedTitle "$UnfoldedTitle"
    --ClosureTitle "$ClosureTitle"
    --ClosureAbsTitle "$ClosureAbsTitle"
    --Chi2Title "$Chi2Title"
    --NormalizeAbsolute "$NormalizeAbsolute"
    --ZoomXMin "$ZoomXMin"
    --ZoomXMax "$ZoomXMax"
    --Chi2RemoveErrorBars "$Chi2RemoveErrorBars"
    --DoRatioToFirst "$DoRatioToFirst"
    --RatioYMin "$RatioYMin"
    --RatioYMax "$RatioYMax"
    --LogYUnfolded "$LogYUnfolded"
    --RequireAllMethods "$RequireAllMethods"
    --XMin "$XMin"
    --XMax "$XMax"
    --YMinUnfolded "$YMinUnfolded"
    --YMaxUnfolded "$YMaxUnfolded"
    --YMinClosure "$YMinClosure"
    --YMaxClosure "$YMaxClosure"
    --YMinClosureAbs "$YMinClosureAbs"
    --YMaxClosureAbs "$YMaxClosureAbs"
    --YMinChi2 "$YMinChi2"
    --YMaxChi2 "$YMaxChi2"
  )

  if [ -n "$Colors" ]; then
    cmd+=(--Colors "$Colors")
  fi

  echo "Executing >>>>>>"
  echo "${cmd[*]}"
  "${cmd[@]}"
done
