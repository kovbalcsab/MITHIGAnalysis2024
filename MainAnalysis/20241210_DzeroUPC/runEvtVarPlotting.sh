#!/bin/bash

SampleSettingCard=${1}
plottingDir=$(jq -r '.plottingDir' $SampleSettingCard)
mkdir -p $plottingDir
cp $SampleSettingCard $plottingDir/plottingEvtVarConfig.json
SampleSettingCard=$plottingDir/plottingEvtVarConfig.json

jq -c '.Plots[]' $SampleSettingCard | while read Plot; do
	PlotBaseName=$(echo $Plot | jq -r '.PlotBaseName')
	InputFileNames=$(echo $Plot | jq -r '.InputFileNames')
	HistNames=$(echo $Plot | jq -r '.HistNames')
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
	mkdir -p $plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/
	OutputFileName=$plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${PlotBaseName}

	cmd="./PlotMultiEvtVars --InputFileNames $InputFileNames  --OutputFileName $OutputFileName \
			--HistNames $HistNames \
			--Labels $Labels \
			--MinDzeroPT $MinDzeroPT \
			--MaxDzeroPT $MaxDzeroPT \
			--MinDzeroY $MinDzeroY \
			--MaxDzeroY $MaxDzeroY \
			--IsGammaN $IsGammaN"
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

	echo "Executing >>>>>>"
	echo $cmd

	PlotLogName="${PlotBaseName/.root/.log}"

	echo $cmd > $plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${PlotLogName}
	( $cmd >> $plottingDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${PlotLogName} ) &
	sleep 0.1
done

sleep 1
while pgrep -x "PlotMultiEvtVars" > /dev/null; do
  sleep 2
done
wait
exit 0
