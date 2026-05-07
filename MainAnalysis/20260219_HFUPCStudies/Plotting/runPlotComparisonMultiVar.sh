#!/bin/bash

SampleSettingCard=${1}
plottingDir=$(jq -r '.plottingDir' $SampleSettingCard)
mkdir -p $plottingDir
cp $SampleSettingCard $plottingDir/plottingEvtVarConfig.json
SampleSettingCard=$plottingDir/plottingEvtVarConfig.json

jq -c '.Plots[]' $SampleSettingCard | while read Plot; do
  PlotBaseName=$(echo $Plot | jq -r '.PlotBaseName')
  InputFileNames=$(echo $Plot | jq -r '.InputFileNames')
  VarNames=$(echo $Plot | jq -r '.VarNames')
  IsTreeVar=$(echo $Plot | jq -r '.IsTreeVar')
  Labels=$(echo $Plot | jq -r '.Labels')
  MinDzeroPT=$(echo $Plot | jq -r '.MinDzeroPT')
  MaxDzeroPT=$(echo $Plot | jq -r '.MaxDzeroPT')
  MinDzeroY=$(echo $Plot | jq -r '.MinDzeroY')
  MaxDzeroY=$(echo $Plot | jq -r '.MaxDzeroY')
  IsGammaN=$(echo $Plot | jq -r '.IsGammaN')
  xTitle=$(echo $Plot | jq -r '.xTitle')
  yTitle=$(echo $Plot | jq -r '.yTitle')
  plotTitle=$(echo $Plot | jq -r '.plotTitle')
  logY=$(echo $Plot | jq -r '.logY')
  xMin=$(echo $Plot | jq -r '.xMin')
  xMax=$(echo $Plot | jq -r '.xMax')
  yMin=$(echo $Plot | jq -r '.yMin')
  yMax=$(echo $Plot | jq -r '.yMax')
  normalizeToUnity=$(echo $Plot | jq -r '.normalizeToUnity')
  rebinFactor=$(echo $Plot | jq -r '.rebinFactor')
  nBins=$(echo $Plot | jq -r '.NBins')
  doRatio=$(echo $Plot | jq -r '.doRatio')
  mkdir -p $plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/
  OutputFileName=$plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${PlotBaseName}

  cmd="./PlotComparisonMultiVar --InputFileNames $InputFileNames  --OutputFileName $OutputFileName \
			--VarNames $VarNames \
			--IsTreeVar $IsTreeVar \
			--Labels $Labels"
  [ $MinDzeroPT != null ] && cmd="$cmd --MinDzeroPT $MinDzeroPT"
  [ $MaxDzeroPT != null ] && cmd="$cmd --MaxDzeroPT $MaxDzeroPT"
  [ $MinDzeroY != null ] && cmd="$cmd --MinDzeroY $MinDzeroY"
  [ $MaxDzeroY != null ] && cmd="$cmd --MaxDzeroY $MaxDzeroY"
  [ $IsGammaN != null ] && cmd="$cmd --IsGammaN $IsGammaN"
  [ $xTitle != null ] && cmd="$cmd --XTitle $xTitle"
  [ $yTitle != null ] && cmd="$cmd --YTitle $yTitle"
  [ $plotTitle != null ] && cmd="$cmd --PlotTitle $plotTitle"
  [ $logY != null ] && cmd="$cmd --LogY $logY"
  [ $xMin != null ] && cmd="$cmd --XMin $xMin"
  [ $xMax != null ] && cmd="$cmd --XMax $xMax"
  [ $yMin != null ] && cmd="$cmd --YMin $yMin"
  [ $yMax != null ] && cmd="$cmd --YMax $yMax"
  [ $normalizeToUnity != null ] && cmd="$cmd --Normalize $normalizeToUnity"
  [ $rebinFactor != null ] && cmd="$cmd --Rebin $rebinFactor"
  [ $nBins != null ] && cmd="$cmd --NBins $nBins"
  [ $doRatio != null ] && cmd="$cmd --DoRatio $doRatio"

  echo "Executing >>>>>>"
  echo $cmd

  PlotLogName="${PlotBaseName/.root/.log}"
  PlotLogPath="$plottingDir"
  if [[ $MinDzeroPT != null && $MaxDzeroPT != null ]]; then PlotLogPath+="/pt${MinDzeroPT}-${MaxDzeroPT}"; fi
  if [[ $MinDzeroY != null && $MaxDzeroY != null ]]; then PlotLogPath+="_y${MinDzeroY}-${MaxDzeroY}"; fi
  if [[ $IsGammaN != null ]]; then PlotLogPath+="_IsGammaN${IsGammaN}"; fi

  echo $cmd >$PlotLogPath/$PlotLogName
  ($cmd >>$PlotLogPath/$PlotLogName) &
  sleep 0.1
done

sleep 1
while pgrep -x "PlotComparisonMultiVar" >/dev/null; do
  sleep 2
done
wait
exit 0
