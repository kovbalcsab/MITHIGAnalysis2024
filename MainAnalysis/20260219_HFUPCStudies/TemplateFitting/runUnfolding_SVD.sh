#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-HFMinusData}" # HFMinusData(0), HFPlusMC(1), HFPlusData(2)
DATA_QUARTER="${2:-0}"

if ! [[ "$DATA_QUARTER" =~ ^[0-3]$ ]]; then
    echo "Data quarter must be in [0,3]." >&2
    echo "Usage: $0 [HFMinusData|HFPlusData|HFPlusMC] [DataQuarter:0-3]" >&2
    exit 1
fi

DATA_FILE=""
NOISE_FILE=""
VAR_NAME=""
OUTPUT_FILE=""

case "$MODE" in
    HFMinusData|MinusData|0)
        DATA_FILE="UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root"
        NOISE_FILE="EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root"
        VAR_NAME="HFEMaxMinus_forest"
        OUTPUT_FILE="TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest_SVD.root"
        ;;
    HFPlusMC|MC|1)
        DATA_FILE="MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root"
        NOISE_FILE="EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root"
        VAR_NAME="HFEMaxPlus_forest"
        OUTPUT_FILE="TemplateFitting/output_unfolding/unfoldingMC_HFEMaxPlus_forest_SVD.root"
        ;;
    HFPlusData|PlusData|2)
        DATA_FILE="UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN1/Data.root"
        NOISE_FILE="EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root"
        VAR_NAME="HFEMaxPlus_forest"
        OUTPUT_FILE="TemplateFitting/output_unfolding/unfolding_HFEMaxPlus_forest_SVD.root"
        ;;
    *)
        echo "Usage: $0 [HFMinusData|HFPlusData|HFPlusMC] [DataQuarter:0-3]" >&2
        exit 1
        ;;
esac

./TemplateFitting/Unfold_SVD \
    --DataFile "$DATA_FILE" \
    --NoiseFile "$NOISE_FILE" \
    --VarDataName "$VAR_NAME" \
    --VarNoiseName "$VAR_NAME" \
    --XMin 0 --XMax 200 --BinsPerGeV 4 \
    --Iterations 6 --KTermMin 60 --KTermMax 200 \
    --DataQuarter "$DATA_QUARTER" \
    --OutputFileName "$OUTPUT_FILE"
