import json

SampleFileTag = "FFT_MC_deconvolution_ybinned"
SampleFileName = f"TemplateFitting/SampleConfigs/DeconvolveFFT_{SampleFileTag}.json"

OutputFileDir = f"TemplateFitting/output_unfolding/FFT_{SampleFileTag}"

#OutputFilesArray = [
#    {
#        "DataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root",
#        "KernelFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
#        "VarDataName": "HFEMaxMinus_forest",
#        "VarKernelName": "HFEMaxMinus_forest",
#        "OutputFileName": f"{OutputFileDir}/deconvolveFFT_HFEMaxMinus_forest_Wiener.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 10,
#        "DampingMin": 0.05,
#        "DampingMax": 10,
#        "LogDamping": 1,
#        "DataQuarter": 0
#    },
#    {
#        "DataFile": "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root",
#        "KernelFile": "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root",
#        "VarDataName": "HFEMaxPlus_forest",
#        "VarKernelName": "HFEMaxPlus_forest",
#        "OutputFileName": f"{OutputFileDir}/deconvolveFFTMC_HFEMaxPlus_forest_Wiener.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 10,
#        "DampingMin": 0.05,
#        "DampingMax": 10,
#        "LogDamping": 1,
#        "DataQuarter": 0
#    },
#    {
#        "DataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN1/Data.root",
#        "KernelFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
#        "VarDataName": "HFEMaxPlus_forest",
#        "VarKernelName": "HFEMaxPlus_forest",
#        "OutputFileName": f"{OutputFileDir}/deconvolveFFT_HFEMaxPlus_forest_Wiener.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 10,
#        "DampingMin": 0.05,
#        "DampingMax": 10,
#        "LogDamping": 1,
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
#        "DataFile": f"UPCDataSignal/output_UPCDataSignal/test/pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}_IsGammaN1/Data.root",
#        "KernelFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
#        "VarDataName": "HFEMaxPlus_forest",
#        "VarKernelName": "HFEMaxPlus_forest",
#        "OutputFileName": f"{OutputFileDir}/deconvolveFFT_HFEMaxPlus_forest_pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 10,
#        "DampingMin": 0.05,
#        "DampingMax": 10,
#        "LogDamping": 1,
#        "DataQuarter": -1
#    })
#    OutputFilesArray.append({
#        "DataFile": f"UPCDataSignal/output_UPCDataSignal/test/pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}_IsGammaN0/Data.root",
#        "KernelFile": "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root",
#        "VarDataName": "HFEMaxMinus_forest",
#        "VarKernelName": "HFEMaxMinus_forest",
#        "OutputFileName": f"{OutputFileDir}/deconvolveFFT_HFEMaxMinus_forest_pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}.root",
#        "XMin": 0,
#        "XMax": 200,
#        "BinsPerGeV": 4,
#        "Iterations": 10,
#        "DampingMin": 0.05,
#        "DampingMax": 10,
#        "LogDamping": 1,
#        "DataQuarter": -1
#    })
    OutputFilesArray.append({
        "DataFile": f"MCSignal/output_MCSignal/test/pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}_IsGammaN1/SignalMC.root",
        "KernelFile": "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root",
        "VarDataName": "HFEMaxPlus_forest",
        "VarKernelName": "HFEMaxPlus_forest",
        "OutputFileName": f"{OutputFileDir}/deconvolveFFTMC_HFEMaxPlus_forest_pt{pt_y_bins[0]}-{pt_y_bins[1]}_y{pt_y_bins[2]}-{pt_y_bins[3]}.root",
        "XMin": 0,
        "XMax": 200,
        "BinsPerGeV": 4,
        "Iterations": 10,
        "DampingMin": 0.05,
        "DampingMax": 10,
        "LogDamping": 1,
        "DataQuarter": -1
    })

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)

