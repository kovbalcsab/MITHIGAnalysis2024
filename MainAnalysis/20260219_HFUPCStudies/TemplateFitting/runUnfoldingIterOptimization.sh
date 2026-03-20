#!/usr/bin/env bash

./TemplateFitting/unfolding_iter_optimization \
    --InputFile TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest.root \
    --OutputFile TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest_iterChi2.root \
    --MeasuredHist hMeasured \
    --RefoldedPrefix hRefolded_iter \
    --UnfoldedPrefix hUnfolded_iter \
    --MaxIterations 10 \
    --UnfoldedRatioDenominatorIteration 4
