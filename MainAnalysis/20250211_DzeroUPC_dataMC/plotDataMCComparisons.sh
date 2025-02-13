PlotSettingCard=${1}
PlotDir=$(jq -r '.PlotDir' $PlotSettingCard)
InputDir=$(jq -r '.InputDir' $PlotSettingCard)
mkdir -p $PlotDir

mkdir $PlotDir/HFEmaxPlus_vs_EvtMult
mkdir $PlotDir/HFEmaxMinus_vs_EvtMult
mkdir $PlotDir/HFEmaxPlus
mkdir $PlotDir/HFEmaxMinus
mkdir $PlotDir/EvtMult
mkdir $PlotDir/nVtx
mkdir $PlotDir/VX
mkdir $PlotDir/VY
mkdir $PlotDir/VZ
mkdir $PlotDir/trkPt_vs_trkEta
mkdir $PlotDir/trkPt
mkdir $PlotDir/trkEta
mkdir $PlotDir/Dalpha
mkdir $PlotDir/Dchi2cl
mkdir $PlotDir/Ddtheta
mkdir $PlotDir/DsvpvDisErr
mkdir $PlotDir/DsvpvDistance
mkdir $PlotDir/DsvpvSig
mkdir $PlotDir/Dmass
mkdir $PlotDir/Dtrk1Pt_vs_Dtrk2Pt
mkdir $PlotDir/Dpt_vs_Dy
mkdir $PlotDir/Dtrk1Pt
mkdir $PlotDir/Dtrk2Pt
mkdir $PlotDir/Dpt
mkdir $PlotDir/Dy

cp $PlotSettingCard $PlotDir/plotConfig.json
PlotSettingCard=$PlotDir/plotConfig.json

echo "" > $PlotDir/plot.log

jq -c '.Plots[]' $PlotSettingCard | while read Plot; do
	MinDzeroPT=$(echo $Plot | jq -r '.MinDzeroPT')
	MaxDzeroPT=$(echo $Plot | jq -r '.MaxDzeroPT')
	MinDzeroY=$(echo $Plot | jq -r '.MinDzeroY')
	MaxDzeroY=$(echo $Plot | jq -r '.MaxDzeroY')
	IsGammaN=$(echo $Plot | jq -r '.IsGammaN')
	fileType=$(echo $Plot | jq -r '.fileType')
	fileType2=$(echo $Plot | jq -r '.fileType2')
	HFEMax=$(echo $Plot | jq -r '.HFEMax')

  cmd="./PlotDataMCComparisons --PlotDir $PlotDir --InputDir $InputDir --MinDzeroPT $MinDzeroPT --MaxDzeroPT $MaxDzeroPT --MinDzeroY $MinDzeroY --MaxDzeroY $MaxDzeroY --IsGammaN $IsGammaN --fileType $fileType"
  [ "$fileType2" != "null" ] && cmd="$cmd --fileType2 $fileType2"
  [ "$HFEMax" != "null" ] && cmd="$cmd --HFEMax $HFEMax"
  

  echo "Executing >>>>>>"
  echo $cmd

	$cmd >> $PlotDir/plot.log

done
