#!/bin/bash

# Create output folder if it doesn't exist
mkdir -p MCSignal/output_MCSignal

SampleSettingCard=${1}
OutputFileDir=$(jq -r '.OutputFileDir' $SampleSettingCard)
mkdir -p $OutputFileDir
cp $SampleSettingCard $OutputFileDir/sampleConfig.json
SampleSettingCard=$OutputFileDir/sampleConfig.json

jq -c '.OutputFiles[]' $SampleSettingCard | while read OutFile; do
    FileName=$(echo $OutFile | jq -r '.FileName')
	Input=$(echo $OutFile | jq -r '.Input')
    DzeroPTMin=$(echo $OutFile | jq -r '.DzeroPTMin')
    DzeroPTMax=$(echo $OutFile | jq -r '.DzeroPTMax')
    DzeroYMin=$(echo $OutFile | jq -r '.DzeroYMin')
    DzeroYMax=$(echo $OutFile | jq -r '.DzeroYMax')
    isResolved=$(echo $OutFile | jq -r '.isResolved')
    useClusterCompatibilityFilter=$(echo $OutFile | jq -r '.useClusterCompatibilityFilter')
    useGenDSignalEventSelection=$(echo $OutFile | jq -r '.useGenDSignalEventSelection')
    useRecoDSignalEventSelection=$(echo $OutFile | jq -r '.useRecoDSignalEventSelection')
    IsGammaN=$(echo $OutFile | jq -r '.IsGammaN')
    DoSystD=$(echo $OutFile | jq -r '.DoSystD')
	
	mkdir -p $OutputFileDir/pt${DzeroPTMin}-${DzeroPTMax}_y${DzeroYMin}-${DzeroYMax}_IsGammaN${IsGammaN}/
	Output=$OutputFileDir/pt${DzeroPTMin}-${DzeroPTMax}_y${DzeroYMin}-${DzeroYMax}_IsGammaN${IsGammaN}/${FileName}

	cmd="./ExtractSignalMCDistributons --Input $Input  --Output $Output \
			--DzeroPTMin $DzeroPTMin \
			--DzeroPTMax $DzeroPTMax \
			--DzeroYMin $DzeroYMin \
			--DzeroYMax $DzeroYMax \
			--IsGammaN $IsGammaN"
			[ $isResolved != null ] && cmd="$cmd --isResolved $isResolved"
            [ $useClusterCompatibilityFilter != null ] && cmd="$cmd --useClusterCompatibilityFilter $useClusterCompatibilityFilter"
            [ $useGenDSignalEventSelection != null ] && cmd="$cmd --useGenDSignalEventSelection $useGenDSignalEventSelection"
            [ $useRecoDSignalEventSelection != null ] && cmd="$cmd --useRecoDSignalEventSelection $useRecoDSignalEventSelection"
			[ $DoSystD != null ] && cmd="$cmd --DoSystD $DoSystD"
			
	echo "Executing >>>>>>"
	echo $cmd

	LogName="${FileName/.root/.log}"

	echo $cmd > $OutputFileDir/pt${DzeroPTMin}-${DzeroPTMax}_y${DzeroYMin}-${DzeroYMax}_IsGammaN${IsGammaN}/${LogName}
	( $cmd >> $OutputFileDir/pt${DzeroPTMin}-${DzeroPTMax}_y${DzeroYMin}-${DzeroYMax}_IsGammaN${IsGammaN}/${LogName} ) &
	sleep 2
done

sleep 1
while pgrep -x "ExtractSignalMCDistributons" > /dev/null; do
  sleep 2
done
wait
exit 0
