#source clean.sh
FitSettingCard=${1}
FitDir=$(jq -r '.FitDir' $FitSettingCard)

# Extract the directory path for fullAnalysis
MicroTreeDir=$(jq -r '.MicroTrees[0].dataInput' $FitSettingCard | grep -o '^[^/]*')
mkdir -p $MicroTreeDir
if [[ $FitSettingCard == *"fitSettings/"* ]]; then
	cp $FitSettingCard $MicroTreeDir/${FitDir}.json
fi

jq -c '.MicroTrees[]' $FitSettingCard | while read MicroTree; do
	dataInput=$(echo $MicroTree | jq -r '.dataInput')
	fitmcInputs=$(echo $MicroTree | jq -r '.fitmcInputs')
	sigswpInputs=$(echo $MicroTree | jq -r '.sigswpInputs')
	KKmcInputs=$(echo $MicroTree | jq -r '.KKmcInputs')
	pipimcInputs=$(echo $MicroTree | jq -r '.pipimcInputs')
	neventsInput=$(echo $MicroTree | jq -r '.neventsInput')
	effmcInput=$(echo $MicroTree | jq -r '.effmcInput')
	doSyst_sig=$(echo $MicroTree | jq -r '.doSyst_sig')
	doSyst_comb=$(echo $MicroTree | jq -r '.doSyst_comb')
	doPkkk=$(echo $MicroTree | jq -r '.doPkkk')
	doPkpp=$(echo $MicroTree | jq -r '.doPkpp')
	RstDir=$(echo $MicroTree | jq -r '.RstDir')
	if [[ $RstDir == "null" ]]; then
		RstDir=$(dirname "$dataInput")
		RstDir=${RstDir}/${FitDir}/
	fi
	mkdir -p $RstDir
	if [[ $FitSettingCard == *"fitSettings/"* ]]; then
		# Create a new JSON with FitDir and the current MicroTree
		newJson=$(jq -n --arg FitDir "$FitDir" --argjson MicroTree "$MicroTree" \
				'{FitDir: $FitDir, MicroTrees: [$MicroTree]}')
		echo "$newJson" > "${RstDir}/fitConfig.json"
		echo "Created JSON: ${RstDir}/fitConfig.json"
		# cp $FitSettingCard $RstDir/fitConfig.json
	fi

  # Build the command dynamically
  cmd="./MassFit --dataInput $dataInput --mcInputs $fitmcInputs"
  [ "$sigswpInputs" != "null" ] && cmd="$cmd --sigswpInputs $sigswpInputs"
  [ "$KKmcInputs" != "null" ] && cmd="$cmd --KKmcInputs $KKmcInputs"
  [ "$pipimcInputs" != "null" ] && cmd="$cmd --pipimcInputs $pipimcInputs"
  [ "$neventsInput" != "null" ] && cmd="$cmd --neventsInput $neventsInput"
  [ "$doSyst_sig" != "null" ] && cmd="$cmd --doSyst_sig $doSyst_sig"
  [ "$doSyst_comb" != "null" ] && cmd="$cmd --doSyst_comb $doSyst_comb"
  [ "$doPkkk" != "null" ] && cmd="$cmd --doPkkk $doPkkk"
  [ "$doPkpp" != "null" ] && cmd="$cmd --doPkpp $doPkpp"
  cmd="$cmd --Output fit.root --RstDir $RstDir"

  echo "Executing >>>>>>"
  echo $cmd
  echo "Executing >>>>>>" > $RstDir/fit.log
  echo $cmd >> $RstDir/fit.log

  # Execute the command and log output
  $cmd >> $RstDir/fit.log

  ##### DEACTIVATED CORRECTEDYIELDS PROCEDURE
  # echo '================================================' >> $RstDir/fit.log
  # echo '= Getting corrected yields:' >> $RstDir/fit.log
  # echo '================================================' >> $RstDir/fit.log
  # cmd="./CorrectedYields --rawYieldInput $RstDir/fit.root --effInput $effmcInput --Output $RstDir/correctedYields.md"

  # echo "Executing >>>>>>"
  # echo $cmd
  # echo "Executing >>>>>>" >> $RstDir/fit.log
  # echo $cmd >> $RstDir/fit.log

  # $cmd >> $RstDir/fit.log
done
