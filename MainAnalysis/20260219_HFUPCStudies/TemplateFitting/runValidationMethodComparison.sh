#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-HFMinusData}" # HFMinusData(0), HFPlusMC(1), HFPlusData(2)

INPUT_BAYES=""
INPUT_SVD=""
INPUT_FFT=""
OUTPUT_PREFIX=""

case "$MODE" in
    HFMinusData|MinusData|0)
        INPUT_BAYES="TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest_iterChi2.root"
        INPUT_SVD="TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest_SVD_validation.root"
        INPUT_FFT="TemplateFitting/output_unfolding/deconvolveFFT_HFEMaxMinus_forest_validation.root"
        OUTPUT_PREFIX="TemplateFitting/output_unfolding/validation_method_comparison_selected_HFEMaxMinusData"
        ;;
    HFPlusMC|MC|1)
        INPUT_BAYES="TemplateFitting/output_unfolding/unfoldingMC_HFEMaxPlus_forest_iterChi2.root"
        INPUT_SVD="TemplateFitting/output_unfolding/unfoldingMC_HFEMaxPlus_forest_SVD_validation.root"
        INPUT_FFT="TemplateFitting/output_unfolding/deconvolveFFTMC_HFEMaxPlus_forest_validation.root"
        OUTPUT_PREFIX="TemplateFitting/output_unfolding/validation_method_comparison_selected_HFEMaxPlusMC"
        ;;
    HFPlusData|PlusData|2)
        INPUT_BAYES="TemplateFitting/output_unfolding/unfolding_HFEMaxPlus_forest_iterChi2.root"
        INPUT_SVD="TemplateFitting/output_unfolding/unfolding_HFEMaxPlus_forest_SVD_validation.root"
        INPUT_FFT="TemplateFitting/output_unfolding/deconvolveFFT_HFEMaxPlus_forest_validation.root"
        OUTPUT_PREFIX="TemplateFitting/output_unfolding/validation_method_comparison_selected_HFEMaxPlusData"
        ;;
    *)
        echo "Usage: $0 [HFMinusData|HFPlusData|HFPlusMC]" >&2
        exit 1
        ;;
esac

OUTPUT_UNFOLDED="${OUTPUT_PREFIX}_unfolded.pdf"
OUTPUT_CLOSURE_RATIO="${OUTPUT_PREFIX}_closure_ratio.pdf"
OUTPUT_CLOSURE_ABS="${OUTPUT_PREFIX}_closure_absolute.pdf"
OUTPUT_CLOSURE_ABS_ZOOM="${OUTPUT_PREFIX}_closure_absolute_x0_25.pdf"
OUTPUT_CHI2="${OUTPUT_PREFIX}_chi2ndf.pdf"

./TemplateFitting/plot_validation_method_comparison \
    --InputFileNames "$INPUT_BAYES","$INPUT_SVD","$INPUT_FFT" \
    --SelectedIterations 4,2,4 \
    --Labels Bayes~Unfold,SVD~Unfold,FFT~Deconvolution \
    --OutputPrefix "$OUTPUT_PREFIX" \
    --OutputUnfoldedFileName "$OUTPUT_UNFOLDED" \
    --OutputClosureRatioFileName "$OUTPUT_CLOSURE_RATIO" \
    --OutputClosureAbsoluteFileName "$OUTPUT_CLOSURE_ABS" \
    --OutputClosureAbsoluteZoomFileName "$OUTPUT_CLOSURE_ABS_ZOOM" \
    --OutputChi2FileName "$OUTPUT_CHI2" \
    --UnfoldedHistPrefix hUnfoldedAbs_iter \
    --ClosureHistPrefix hRefoldedToMeasuredRatio_iter \
    --RefoldedAbsHistPrefix hRefoldedAbs_iter \
    --MeasuredHistName hMeasuredForOverlay \
    --Chi2HistName hChi2NDFVsIter \
    --UnfoldedXTitle HF~Energy~[GeV] \
    --UnfoldedYTitle Entries \
    --ClosureXTitle HF~Energy~[GeV] \
    --ClosureYTitle Refolded/Measured \
    --ClosureAbsXTitle HF~Energy~[GeV] \
    --ClosureAbsYTitle 1/N~dN/dE \
    --Chi2XTitle Regularization~iteration \
    --Chi2YTitle \#chi^{2}/NDF \
    --UnfoldedTitle Selected~iteration~unfolded~comparison \
    --ClosureTitle Selected~iteration~closure~ratio~comparison \
    --ClosureAbsTitle Selected~iteration~measured~vs~refolded~comparison \
    --Chi2Title Method~comparison~of~#chi^{2}/NDF~scan \
    --NormalizeAbsolute true \
    --ZoomXMin 0 \
    --ZoomXMax 25 \
    --Chi2RemoveErrorBars true \
    --DoRatioToFirst true \
    --RatioYMin 0.0 \
    --RatioYMax 2.0 \
    --LogYUnfolded true
