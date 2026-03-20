./TemplateFitting/Unfold \
    --DataFile UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root \
    --NoiseFile   EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root \
    --VarDataName  HFEMaxMinus_forest \
    --VarNoiseName    HFEMaxMinus_forest \
    --XMin 0 --XMax 200 --BinsPerGeV 4 \
    --Iterations 10 \
    --OutputFileName   TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest.root