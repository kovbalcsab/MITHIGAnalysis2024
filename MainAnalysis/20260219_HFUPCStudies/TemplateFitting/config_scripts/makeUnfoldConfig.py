import json

SampleFileTag = "data_unfolding_ybinned"
SampleFileName = f"TemplateFitting/SampleConfigs/Unfold_{SampleFileTag}.json"

OutputFileDir = f"TemplateFitting/output_unfolding/{SampleFileTag}"

#OutputFilesArray = [
#    {
#        "DataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root",
#        "NoiseFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
#        "VarDataName": "HFEMaxMinus_forest",
#        "VarNoiseName": "HFEMaxMinus_forest",
#        "OutputFileName": f"{OutputFileDir}/unfolding_HFEMaxMinus_forest.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 20,
#        "DataQuarter": 0
#    },
#    {
#        "DataFile": "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root",
#        "NoiseFile": "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root",
#        "VarDataName": "HFEMaxPlus_forest",
#        "VarNoiseName": "HFEMaxPlus_forest",
#        "OutputFileName": f"{OutputFileDir}/unfoldingMC_HFEMaxPlus_forest.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 20,
#        "DataQuarter": 0
#    },
#    {
#        "DataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN1/Data.root",
#        "NoiseFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
#        "VarDataName": "HFEMaxPlus_forest",
#        "VarNoiseName": "HFEMaxPlus_forest",
#        "OutputFileName": f"{OutputFileDir}/unfolding_HFEMaxPlus_forest.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 20,
#        "DataQuarter": 0
#    }
#]

PtYBins = [ [2,5,-2,-1],
            [2,5,-1,0],
            [2,5,0,1],
            [2,5,1,2],
            [2,5,-2,2]  ]

OutputFilesArray = []

for pt_y_bins in PtYBins:
#    OutputFilesArray.append({
#        "DataFile": f"MCSignal/output_MCSignal/test/pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}_IsGammaN1/SignalMC.root",
#        "NoiseFile": "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root",
#        "VarDataName": "HFEMaxPlus_forest",
#        "VarNoiseName": "HFEMaxPlus_forest",
#        "OutputFileName": f"{OutputFileDir}/unfoldingMC_HFEMaxPlus_forest_pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 10,
#        "DataQuarter": -1
#    })
    OutputFilesArray.append({
        "DataFile": f"UPCDataSignal/output_UPCDataSignal/test/pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}_IsGammaN0/Data.root",
        "NoiseFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
        "VarDataName": "HFEMaxMinus_forest",
        "VarNoiseName": "HFEMaxMinus_forest",
        "OutputFileName": f"{OutputFileDir}/unfolding_HFEMaxMinus_forest_pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}.root",
        "XMin": 0,
        "XMax": 200,
        "BinsPerGeV": 4,
        "Iterations": 10,
        "DataQuarter": -1
    })
    OutputFilesArray.append({
        "DataFile": f"UPCDataSignal/output_UPCDataSignal/test/pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}_IsGammaN1/Data.root",
        "NoiseFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
        "VarDataName": "HFEMaxPlus_forest",
        "VarNoiseName": "HFEMaxPlus_forest",
        "OutputFileName": f"{OutputFileDir}/unfolding_HFEMaxPlus_forest_pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}.root",
        "XMin": 0,
        "XMax": 200,
        "BinsPerGeV": 4,
        "Iterations": 10,
        "DataQuarter": -1
    })

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)

