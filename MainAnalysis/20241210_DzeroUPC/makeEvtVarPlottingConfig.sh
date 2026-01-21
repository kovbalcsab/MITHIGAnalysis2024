#!/bin/bash

PT_Y_BINS=(
  2   5   -2   -1
  2   5   -1    0
  2   5    0    1
  2   5    1    2
  2   5    -2   2
  5   12   -2   2
)

ISGAMMA=(0 1)

JSON_OUTPUT_FILE="evtVarPlotConfigs/CCFInv_20260116_20GeV.json"
PLOTTING_DIR="evtVarPlots_CCFInv_20260116_20GeV"

# Directories containing the input MicroTrees
#INPUTFILES_MicroTreeDirs=(
#    "fullAnalysis_withCCF"
#    "fullAnalysis_noCCF"
#    "fullAnalysis_noCCF"
#    "fullAnalysis_noCCF"
#    "fullAnalysis_withCCF"
#)
INPUTFILES_MicroTreeDirs=(
    "systRapGap20GeV_withCCF"
    "systRapGap20GeV_noCCF"
    "systRapGap20GeV_noCCF"
    "systRapGap20GeV_noCCF"
    "systRapGap20GeV_withCCF"
)

# Input root filenames inside the right bin folders
INPUTFILES_RootFiles=(
    "Data.root"
    "Data.root"
    "MC.root"
    "MC_inclusive.root"
    "MC.root"
)

# Labels for the different samples
SAMPLELABELS=(
    "2025~Data~With~CCF"
    "2025~Data~No~CCF"
    "2023~(Forced)~MC"
    "2023~Inclusive~MC"
    "2025~MC"
)

# Histogram types for each plotted function
HISTNAMES_TYPE=(
    "SignalRecoD"
    "SignalRecoD"
    "SignalRecoD"
    "SignalRecoD"
    "SignalRecoD"
)

# Check if array sizes match
NUM_INPUT_DIRS=${#INPUTFILES_MicroTreeDirs[@]}
NUM_INPUT_FILES=${#INPUTFILES_RootFiles[@]}
NUM_LABELS=${#SAMPLELABELS[@]}
NUM_HIST_TYPES=${#HISTNAMES_TYPE[@]}
if [[ $NUM_INPUT_DIRS -ne $NUM_INPUT_FILES || $NUM_INPUT_DIRS -ne $NUM_LABELS || $NUM_INPUT_DIRS -ne $NUM_HIST_TYPES ]]; then
    echo "Error: Array sizes do not match!"
    exit 1
fi

# Histogram base names to be plotted
HISTNAMES_BASE=(
    "hnTrackInAcceptanceHP"
    "hnVtx"
    "hVX"
    "hVY"
    "hVZ"
    "hHFEmaxMinus"
    "hHFEmaxPlus"
    "hDNumPassingDcuts"
)

HIST_XTITLE=(
    "Number~of~HP~tracks"
    "Vertex~count"
    "Primary~vertex~X~position~(cm)"
    "Primary~vertex~Y~position~(cm)"
    "Primary~vertex~Z~position~(cm)"
    "HF~E_{max}~Minus~(GeV)"
    "HF~E_{max}~Plus~(GeV)"
    "Number~of~D~candidates~passing~D~cuts"
)

HIST_YTITLE=(
    "Normalized~counts"
    "Normalized~counts"
    "Normalized~counts"
    "Normalized~counts"
    "Normalized~counts"
    "Normalized~counts"
    "Normalized~counts"
    "Normalized~counts"
)

HIST_LOGY=(
    "0"
    "1"
    "1"
    "1"
    "1"
    "1"
    "1"
    "1"
)

HIST_XLIMS=(
    0 60
    -0.5 4.5
    -0.2 0.2
    -0.2 0.2
    -30 30
    0 100
    0 100
    0 10
)

HIST_YLIMS=(
    0.0001 0.5
    0.0001 1
    0.0001 1
    0.0001 1
    0.0001 1
    0.0001 1
    0.0001 1
    0.0001 1
)

HIST_NORMALIZED=(
    1
    1
    1
    1
    1
    1
    1
    1
)

HIST_REBINFACTOR=(
    1
    1
    1
    1
    1
    1
    1
    1
)

# Check if histogram parameter array sizes match
NUM_HISTOGRAMS=${#HISTNAMES_BASE[@]}
if [[ $NUM_HISTOGRAMS -ne ${#HIST_XTITLE[@]} 
    || $NUM_HISTOGRAMS -ne ${#HIST_YTITLE[@]} 
    || $NUM_HISTOGRAMS -ne ${#HIST_LOGY[@]} 
    || $NUM_HISTOGRAMS -ne $((${#HIST_XLIMS[@]}/2)) 
    || $NUM_HISTOGRAMS -ne $((${#HIST_YLIMS[@]}/2)) 
    || $NUM_HISTOGRAMS -ne ${#HIST_NORMALIZED[@]} 
    || $NUM_HISTOGRAMS -ne ${#HIST_REBINFACTOR[@]} ]]; then
    echo "Error: Histogram parameter array sizes do not match!"
    exit 1
fi

# Create output JSON configuration file
mkdir -p evtVarPlotConfigs
rm -f $JSON_OUTPUT_FILE
# Initialize or clear the output JSON file
echo '{' >> $JSON_OUTPUT_FILE
echo "  \"plottingDir\": \"$PLOTTING_DIR\"," >> $JSON_OUTPUT_FILE
echo '  "Plots": [' >> $JSON_OUTPUT_FILE

# Begin loop over histogram base names
for ((i=0; i<${#HISTNAMES_BASE[@]}; i++)); do
    #Begin loop over pt-y bins
    for ((j=0; j<${#PT_Y_BINS[@]}/4; j++)); do
        for ((k=0; k<${#ISGAMMA[@]}; k++)); do
            PT_MIN=${PT_Y_BINS[j*4]}
            PT_MAX=${PT_Y_BINS[j*4+1]}
            Y_MIN=${PT_Y_BINS[j*4+2]}
            Y_MAX=${PT_Y_BINS[j*4+3]}

            INPUTFILES=""
            HISTNAMES=""
            LABELS=""

            # Loop over input directories/files to build the combined inputs for this plot
            for ((m=0; m<${#INPUTFILES_MicroTreeDirs[@]}; m++)); do
                INPUTFILE="${INPUTFILES_MicroTreeDirs[m]}/pt${PT_MIN}-${PT_MAX}_y${Y_MIN}-${Y_MAX}_IsGammaN${ISGAMMA[k]}/${INPUTFILES_RootFiles[m]}"
                HISTNAME="${HISTNAMES_BASE[i]}${HISTNAMES_TYPE[m]}"
                LABEL="${SAMPLELABELS[m]}"

                if [ -z "$INPUTFILES" ]; then
                    INPUTFILES="$INPUTFILE"
                    HISTNAMES="$HISTNAME"
                    LABELS="$LABEL"
                else
                    INPUTFILES="$INPUTFILES,$INPUTFILE"
                    HISTNAMES="$HISTNAMES,$HISTNAME"
                    LABELS="$LABELS,$LABEL"
                fi
            done

            # Create JSON plot configuration
            PLOTCONFIG_JSON="{
                \"PlotBaseName\": \"${HISTNAMES_BASE[i]}_EvtVarPlot.pdf\",
                \"InputFileNames\": \"$INPUTFILES\",
                \"HistNames\": \"$HISTNAMES\",
                \"Labels\": \"$LABELS\",
                \"MinDzeroPT\": ${PT_MIN},
                \"MaxDzeroPT\": ${PT_MAX},
                \"MinDzeroY\": ${Y_MIN},
                \"MaxDzeroY\": ${Y_MAX},
                \"IsGammaN\": ${ISGAMMA[k]},
                \"xTitle\": \"${HIST_XTITLE[i]}\",
                \"yTitle\": \"${HIST_YTITLE[i]}\",
                \"logY\": ${HIST_LOGY[i]},
                \"xMin\": ${HIST_XLIMS[i*2]},
                \"xMax\": ${HIST_XLIMS[i*2+1]},
                \"yMin\": ${HIST_YLIMS[i*2]},
                \"yMax\": ${HIST_YLIMS[i*2+1]},
                \"normalizeToUnity\": ${HIST_NORMALIZED[i]},
                \"rebinFactor\": ${HIST_REBINFACTOR[i]}
            },"
            echo "$PLOTCONFIG_JSON" >> $JSON_OUTPUT_FILE
            
        done
    done
done
# Remove trailing comma from the last plot configuration
sed -i '$ s/},$/}/' $JSON_OUTPUT_FILE
echo '  ]' >> $JSON_OUTPUT_FILE
echo '}' >> $JSON_OUTPUT_FILE