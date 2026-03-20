./TemplateFitting/Unfold_SVD \
    --DataFile UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root \
    --NoiseFile   EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root \
    --VarDataName  HFEMaxMinus_forest \
    --VarNoiseName    HFEMaxMinus_forest \
    --XMin 0 --XMax 200 --BinsPerGeV 4 \
    --Iterations 20 --KTermMin 30 --KTermMax 200 \
    --OutputFileName   TemplateFitting/output_unfolding/unfolding_HFEMaxMinus_forest_SVD.root