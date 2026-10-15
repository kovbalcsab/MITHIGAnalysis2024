#!/usr/bin/env bash
set -euo pipefail

# Edit these inputs directly in the script.
RunFitter="true"
RunPlotter="true"
InputFile="TemplateFitting_2/outputs/01_09/mixedSignalMC_01_09.root"
IsPlus="false"
OutputFile="TemplateFitting_2/outputs/01_09/fitHFTemplatesMC_minus_output.root"
PlotOutputFile="TemplateFitting_2/outputs/01_09/fitHFTemplatesMC_output_minus_replot.png"
PlotXMax="30"
TemplatePaths=(
  "MCSignal/output_MCSignal/24_MC_inclusive/pt2-5_y-2-2_IsGammaN1/SignalMC.root"
  "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root"
)
TemplateNames=(
  "signal"
  "noise"
)

ScriptDir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RepoDir=$(cd "$ScriptDir/.." && pwd)
BinaryFile="$RepoDir/ExecuteFitHFTemplatesMC"
SourceFile="$RepoDir/TemplateFitting_2/fitHFTemplatesMC.cpp"
PlotBinaryFile="$RepoDir/ExecutePlotHFTemplateFitResult"
PlotSourceFile="$RepoDir/TemplateFitting_2/plotHFTemplateFitResult.cpp"
InfoManagerSource="$RepoDir/src/InfoManager.cpp"

if [ "$#" -ne 0 ]; then
  echo "This script uses the input block at the top. Run it without command-line arguments." >&2
  exit 1
fi

if [ "$RunFitter" != "true" ] && [ "$RunFitter" != "false" ]; then
  echo "RunFitter must be either \"true\" or \"false\"." >&2
  exit 1
fi

if [ "$RunPlotter" != "true" ] && [ "$RunPlotter" != "false" ]; then
  echo "RunPlotter must be either \"true\" or \"false\"." >&2
  exit 1
fi

if [ "$RunFitter" != "true" ] && [ "$RunPlotter" != "true" ]; then
  echo "Nothing to do: both RunFitter and RunPlotter are false." >&2
  exit 1
fi

if [ "$RunFitter" = "true" ]; then
  if [ ! -f "$InputFile" ]; then
    echo "Input file not found: $InputFile" >&2
    exit 1
  fi

  if [ "${#TemplatePaths[@]}" -eq 0 ]; then
    echo "TemplatePaths must not be empty." >&2
    exit 1
  fi

  if [ "${#TemplatePaths[@]}" -ne "${#TemplateNames[@]}" ]; then
    echo "TemplatePaths and TemplateNames must have the same length." >&2
    exit 1
  fi

  for TemplateFile in "${TemplatePaths[@]}"; do
    if [ ! -f "$TemplateFile" ]; then
      echo "Template file not found: $TemplateFile" >&2
      exit 1
    fi
  done

  if [ ! -x "$BinaryFile" ] || [ "$SourceFile" -nt "$BinaryFile" ] || [ "$InfoManagerSource" -nt "$BinaryFile" ]; then
    g++ "$SourceFile" "$InfoManagerSource" \
      -I"$RepoDir" \
      -O3 \
      -I"$RepoDir/../../CommonCode/include" \
      -I"$RepoDir/include" \
      -o "$BinaryFile" \
      $(root-config --cflags --glibs) \
      -lRooFit -lRooFitCore
  fi
fi

if [ "$RunPlotter" = "true" ]; then
  if [ ! -x "$PlotBinaryFile" ] || [ "$PlotSourceFile" -nt "$PlotBinaryFile" ] || [ "$InfoManagerSource" -nt "$PlotBinaryFile" ]; then
    g++ "$PlotSourceFile" "$InfoManagerSource" \
      -I"$RepoDir" \
      -O3 \
      -I"$RepoDir/../../CommonCode/include" \
      -I"$RepoDir/include" \
      -o "$PlotBinaryFile" \
      $(root-config --cflags --glibs) \
      -lRooFit -lRooFitCore
  fi
fi

mkdir -p "$(dirname "$OutputFile")"
mkdir -p "$(dirname "$PlotOutputFile")"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/RooUnfold

if [ "$RunFitter" = "true" ]; then
  TemplatePathsCSV=$(
    IFS=,
    echo "${TemplatePaths[*]}"
  )
  TemplateNamesCSV=$(
    IFS=,
    echo "${TemplateNames[*]}"
  )

  fit_cmd=(
    "$BinaryFile"
    --inputFile "$InputFile"
    --isPlus "$IsPlus"
    --outputFile "$OutputFile"
    --templatePaths "$TemplatePathsCSV"
    --templateNames "$TemplateNamesCSV"
  )

  echo "Executing fitter >>>>>>"
  echo "${fit_cmd[*]}"
  "${fit_cmd[@]}"
fi

if [ "$RunPlotter" = "true" ]; then
  if [ ! -f "$OutputFile" ]; then
    echo "Fit result file not found for plotting: $OutputFile" >&2
    exit 1
  fi

  plot_cmd=(
    "$PlotBinaryFile"
    --resultFile "$OutputFile"
    --outputPlot "$PlotOutputFile"
    --xMax "$PlotXMax"
  )

  echo "Executing plotter >>>>>>"
  echo "${plot_cmd[*]}"
  "${plot_cmd[@]}"
fi
