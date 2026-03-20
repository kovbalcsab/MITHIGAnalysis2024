#!/bin/bash

# Create output folder if it doesn't exist
mkdir -p UPCDataSignal/output_UPCDataSignal

SampleSettingCard=${1}
OutputFileDir=$(jq -r '.OutputFileDir' $SampleSettingCard)
mkdir -p $OutputFileDir
cp $SampleSettingCard $OutputFileDir/sampleConfig.json
SampleSettingCard=$OutputFileDir/sampleConfig.json

jq -c '.OutputFiles[]' $SampleSettingCard | while read OutFile; do
    FileName=$(echo $OutFile | jq -r '.FileName')
	Input=$(echo $OutFile | jq -r '.Input')
    MinDzeroPT=$(echo $OutFile | jq -r '.DzeroPTMin')
    MaxDzeroPT=$(echo $OutFile | jq -r '.DzeroPTMax')
    MinDzeroY=$(echo $OutFile | jq -r '.DzeroYMin')
    MaxDzeroY=$(echo $OutFile | jq -r '.DzeroYMax')
    TriggerChoice=$(echo $OutFile | jq -r '.TriggerChoice')
    IsGammaN=$(echo $OutFile | jq -r '.IsGammaN')
    DoSystD=$(echo $OutFile | jq -r '.DoSystD')
	
	mkdir -p $OutputFileDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/
	Output=$OutputFileDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${FileName}

	cmd="./ExtractUPCDataSignalDistributons --Input $Input  --Output $Output \
			--MinDzeroPT $MinDzeroPT \
			--MaxDzeroPT $MaxDzeroPT \
			--MinDzeroY $MinDzeroY \
			--MaxDzeroY $MaxDzeroY \
			--IsGammaN $IsGammaN"
			[ $TriggerChoice != null ] && cmd="$cmd --TriggerChoice $TriggerChoice"
            [ $DoSystD != null ] && cmd="$cmd --DoSystD $DoSystD"
			
	echo "Executing >>>>>>"
	echo $cmd

	LogName="${FileName/.root/.log}"

	echo $cmd > $OutputFileDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${LogName}
	( $cmd >> $OutputFileDir/pt${MinDzeroPT}-${MaxDzeroPT}_y${MinDzeroY}-${MaxDzeroY}_IsGammaN${IsGammaN}/${LogName} ) &
	sleep 2
done

sleep 1
while pgrep -x "ExtractUPCDataSignalDistributons" > /dev/null; do
  sleep 2
done
wait
exit 0
