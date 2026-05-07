#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF' >&2
Usage:
  runValidationQuarterChi2.sh [HFMinusData|HFPlusData|HFPlusMC]
  runValidationQuarterChi2.sh <BayesFilesCSV> <SVDFilesCSV> <FFTFilesCSV> [OutputPrefix]

Notes:
  - CSV inputs are comma-separated quarter files (typically q0,q1,q2,q3).
  - The script plots hChi2NDFVsIter (chi2/ndf) from each file and overlays quarters per method.
EOF
}

MODE="${1:-HFPlusData}" # HFMinusData(0), HFPlusMC(1), HFPlusData(2)

BASE_BAYES="TemplateFitting/output_unfolding/validation_quater0"
BASE_SVD="TemplateFitting/output_unfolding/SVD_validation_quater0"
BASE_FFT="TemplateFitting/output_unfolding/FFT_validation_quater0"
OUTPUT_PREFIX=""

declare -a BAYES_FILES
declare -a SVD_FILES
declare -a FFT_FILES

join_by_comma() {
  local IFS=','
  echo "$*"
}

if [ "$#" -ge 3 ] && [[ "$1" == *","* || "$1" == *.root ]]; then
  BAYES_CSV="$1"
  SVD_CSV="$2"
  FFT_CSV="$3"
  OUTPUT_PREFIX="${4:-TemplateFitting/output_unfolding/validation_quarter_chi2ndf_custom}"

  mkdir -p "$(dirname "$OUTPUT_PREFIX")"
  ./TemplateFitting/plot_validation_quarter_chi2 \
    --BayesFiles "$BAYES_CSV" \
    --SVDFiles "$SVD_CSV" \
    --FFTFiles "$FFT_CSV" \
    --Chi2HistName hChi2NDFVsIter \
    --MakeMeasuredOnlyPlots true \
    --BayesMeasuredOnlyHistName hChi2NDFMeasuredOnlyVsIter \
    --SVDMeasuredOnlyHistName hChi2NDFMeasuredOnlyVsIter \
    --FFTMeasuredOnlyHistName hChi2NDFVsIter \
    --OutputPrefix "$OUTPUT_PREFIX" \
    --XTitle Regularization~iteration \
    --YTitle \#chi^{2}/NDF \
    --YTitleMeasuredOnly \#chi^{2}/NDF~\(measured~errors~only\) \
    --BayesTitle Bayes~unfolding:~#chi^{2}/NDF~vs~iteration~for~each~validation~quarter \
    --SVDTitle SVD~unfolding:~#chi^{2}/NDF~vs~iteration~for~each~validation~quarter \
    --FFTTitle FFT~deconvolution:~#chi^{2}/NDF~vs~iteration~for~each~validation~quarter \
    --BayesTitleMeasuredOnly Bayes~unfolding:~#chi^{2}/NDF~\(measured~errors~only\)~vs~iteration~for~each~validation~quarter \
    --SVDTitleMeasuredOnly SVD~unfolding:~#chi^{2}/NDF~\(measured~errors~only\)~vs~iteration~for~each~validation~quarter \
    --FFTTitleMeasuredOnly FFT~deconvolution:~#chi^{2}/NDF~\(measured~errors~only\)~vs~iteration~for~each~validation~quarter \
    --RequireAllQuarters true \
    --LogY false
  exit 0
fi

if [ "$#" -gt 1 ] && [[ "$1" == -* ]]; then
  usage
  exit 1
fi

case "$MODE" in
  HFMinusData|MinusData|0)
    for Q in 0 1 2 3; do
      BAYES_FILES+=("${BASE_BAYES}/unfolding_HFEMaxMinus_forest_iterChi2_vs_q${Q}.root")
      SVD_FILES+=("${BASE_SVD}/unfolding_HFEMaxMinus_forest_SVD_validation_vs_q${Q}.root")
      FFT_FILES+=("${BASE_FFT}/deconvolveFFT_HFEMaxMinus_forest_validation_vs_q${Q}.root")
    done
    OUTPUT_PREFIX="TemplateFitting/output_unfolding/validation_quarter_chi2ndf_HFEMaxMinusData/plot"
    ;;
  HFPlusMC|MC|1)
    for Q in 0 1 2 3; do
      BAYES_FILES+=("${BASE_BAYES}/unfoldingMC_HFEMaxPlus_forest_iterChi2_vs_q${Q}.root")
      SVD_FILES+=("${BASE_SVD}/unfoldingMC_HFEMaxPlus_forest_SVD_validation_vs_q${Q}.root")
      FFT_FILES+=("${BASE_FFT}/deconvolveFFTMC_HFEMaxPlus_forest_validation_vs_q${Q}.root")
    done
    OUTPUT_PREFIX="TemplateFitting/output_unfolding/validation_quarter_chi2ndf_HFEMaxPlusMC/plot"
    ;;
  HFPlusData|PlusData|2)
    for Q in 0 1 2 3; do
      BAYES_FILES+=("${BASE_BAYES}/unfolding_HFEMaxPlus_forest_iterChi2_vs_q${Q}.root")
      SVD_FILES+=("${BASE_SVD}/unfolding_HFEMaxPlus_forest_SVD_validation_vs_q${Q}.root")
      FFT_FILES+=("${BASE_FFT}/deconvolveFFT_HFEMaxPlus_forest_validation_vs_q${Q}.root")
    done
    OUTPUT_PREFIX="TemplateFitting/output_unfolding/validation_quarter_chi2ndf_HFEMaxPlusData/plot"
    ;;
  *)
    usage
    exit 1
    ;;
esac

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

./TemplateFitting/plot_validation_quarter_chi2 \
  --BayesFiles "$(join_by_comma "${BAYES_FILES[@]}")" \
  --SVDFiles "$(join_by_comma "${SVD_FILES[@]}")" \
  --FFTFiles "$(join_by_comma "${FFT_FILES[@]}")" \
  --Chi2HistName hChi2NDFVsIter \
  --MakeMeasuredOnlyPlots true \
  --BayesMeasuredOnlyHistName hChi2NDFMeasuredOnlyVsIter \
  --SVDMeasuredOnlyHistName hChi2NDFMeasuredOnlyVsIter \
  --FFTMeasuredOnlyHistName hChi2NDFVsIter \
  --OutputPrefix "$OUTPUT_PREFIX" \
  --XTitle Regularization~iteration \
  --YTitle \#chi^{2}/NDF \
  --YTitleMeasuredOnly \#chi^{2}/NDF~\(measured~errors~only\) \
  --BayesTitle Bayes~unfolding:~#chi^{2}/NDF~vs~iteration~for~each~validation~quarter \
  --SVDTitle SVD~unfolding:~#chi^{2}/NDF~vs~iteration~for~each~validation~quarter \
  --FFTTitle FFT~deconvolution:~#chi^{2}/NDF~vs~iteration~for~each~validation~quarter \
  --BayesTitleMeasuredOnly Bayes~unfolding:~#chi^{2}/NDF~\(measured~errors~only\)~vs~iteration~for~each~validation~quarter \
  --SVDTitleMeasuredOnly SVD~unfolding:~#chi^{2}/NDF~\(measured~errors~only\)~vs~iteration~for~each~validation~quarter \
  --FFTTitleMeasuredOnly FFT~deconvolution:~#chi^{2}/NDF~\(measured~errors~only\)~vs~iteration~for~each~validation~quarter \
  --RequireAllQuarters true \
  --LogY false
