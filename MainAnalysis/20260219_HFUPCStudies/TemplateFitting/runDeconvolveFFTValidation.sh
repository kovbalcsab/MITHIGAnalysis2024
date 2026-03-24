#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-HFMinusData}" # HFMinusData(0), HFPlusMC(1), HFPlusData(2)
TRAINING_QUARTER="${2:--1}"
VALIDATION_QUARTER="${3:--1}"

if ! [[ "$TRAINING_QUARTER" =~ ^(-1|[0-3])$ ]]; then
    echo "Training quarter must be -1 or in [0,3]." >&2
    echo "Usage: $0 [HFMinusData|HFPlusData|HFPlusMC] [TrainingQuarter:-1..3] [ValidationQuarter:-1..3]" >&2
    exit 1
fi

if ! [[ "$VALIDATION_QUARTER" =~ ^(-1|[0-3])$ ]]; then
    echo "Validation quarter must be -1 or in [0,3]." >&2
    echo "Usage: $0 [HFMinusData|HFPlusData|HFPlusMC] [TrainingQuarter:-1..3] [ValidationQuarter:-1..3]" >&2
    exit 1
fi

INPUT_FILE=""
OUTPUT_FILE=""
VALIDATION_DATA_FILE=""
VALIDATION_TREE_NAME="OutputTree"
VALIDATION_VAR_NAME=""

case "$MODE" in
    HFMinusData|MinusData|0)
        INPUT_FILE="TemplateFitting/output_unfolding/deconvolveFFT_HFEMaxMinus_forest_Wiener.root"
        OUTPUT_FILE="TemplateFitting/output_unfolding/deconvolveFFT_HFEMaxMinus_forest_validation.root"
        VALIDATION_DATA_FILE="UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root"
        VALIDATION_VAR_NAME="HFEMaxMinus_forest"
        ;;
    HFPlusMC|MC|1)
        INPUT_FILE="TemplateFitting/output_unfolding/deconvolveFFTMC_HFEMaxPlus_forest_Wiener.root"
        OUTPUT_FILE="TemplateFitting/output_unfolding/deconvolveFFTMC_HFEMaxPlus_forest_validation.root"
        VALIDATION_DATA_FILE="MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root"
        VALIDATION_VAR_NAME="HFEMaxPlus_forest"
        ;;
    HFPlusData|PlusData|2)
        INPUT_FILE="TemplateFitting/output_unfolding/deconvolveFFT_HFEMaxPlus_forest_Wiener.root"
        OUTPUT_FILE="TemplateFitting/output_unfolding/deconvolveFFT_HFEMaxPlus_forest_validation.root"
        VALIDATION_DATA_FILE="UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN1/Data.root"
        VALIDATION_VAR_NAME="HFEMaxPlus_forest"
        ;;
    *)
        echo "Usage: $0 [HFMinusData|HFPlusData|HFPlusMC] [TrainingQuarter:-1..3] [ValidationQuarter:-1..3]" >&2
        exit 1
        ;;
esac

CMD=(
    ./TemplateFitting/deconvolve_fft_validation
    --InputFile "$INPUT_FILE"
    --OutputFile "$OUTPUT_FILE"
    --MeasuredHist hMeasured
    --DeconvolvedPrefix hDeconvolved_iter
    --ReconvolvedPrefix hReconvolved_iter
    --RegularizationTree RegularizationTree
    --DampingBranch DampingFactor
    --MaxIterations 6
    --UnfoldedRatioDenominatorIteration 4
    --TrainingQuarter "$TRAINING_QUARTER"
)

if [ "$VALIDATION_QUARTER" -ge 0 ]; then
    CMD+=(
        --ValidationDataFile "$VALIDATION_DATA_FILE"
        --ValidationTreeName "$VALIDATION_TREE_NAME"
        --ValidationVarName "$VALIDATION_VAR_NAME"
        --ValidationQuarter "$VALIDATION_QUARTER"
    )
fi

"${CMD[@]}"
