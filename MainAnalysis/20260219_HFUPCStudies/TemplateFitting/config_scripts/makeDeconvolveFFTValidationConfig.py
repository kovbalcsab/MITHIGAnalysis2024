import json

SampleFileTag = "validation_quater0"
SampleFileName = f"TemplateFitting/SampleConfigs/DeconvolveFFTValidation_{SampleFileTag}.json"

OutputFileDir = f"TemplateFitting/output_unfolding/FFT_{SampleFileTag}"

base = {
    "MeasuredHist": "hMeasured",
    "DeconvolvedPrefix": "hDeconvolved_iter",
    "ReconvolvedPrefix": "hReconvolved_iter",
    "RegularizationTree": "RegularizationTree",
    "DampingBranch": "DampingFactor",
    "MaxIterations": 10,
    "UnfoldedRatioDenominatorIteration": 4,
    "TrainingQuarter": 0,
    "ValidationTreeName": "OutputTree"
}

OutputFilesArray = []
validation_quarter_array = [0, 1, 2, 3]

for validation_quarter in validation_quarter_array:
    OutputFilesArray.extend([
        {
            **base,
            "ValidationQuarter": validation_quarter,
            "InputFile": f"{OutputFileDir}/deconvolveFFT_HFEMaxMinus_forest_Wiener.root",
            "OutputFile": f"{OutputFileDir}/deconvolveFFT_HFEMaxMinus_forest_validation_vs_q{validation_quarter}.root",
            "ValidationDataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root",
            "ValidationVarName": "HFEMaxMinus_forest"
        },
        {
            **base,
            "ValidationQuarter": validation_quarter,
            "InputFile": f"{OutputFileDir}/deconvolveFFTMC_HFEMaxPlus_forest_Wiener.root",
            "OutputFile": f"{OutputFileDir}/deconvolveFFTMC_HFEMaxPlus_forest_validation_vs_q{validation_quarter}.root",
            "ValidationDataFile": "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root",
            "ValidationVarName": "HFEMaxPlus_forest"
        },
        {
            **base,
            "ValidationQuarter": validation_quarter,
            "InputFile": f"{OutputFileDir}/deconvolveFFT_HFEMaxPlus_forest_Wiener.root",
            "OutputFile": f"{OutputFileDir}/deconvolveFFT_HFEMaxPlus_forest_validation_vs_q{validation_quarter}.root",
            "ValidationDataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN1/Data.root",
            "ValidationVarName": "HFEMaxPlus_forest"
        }
    ])

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)
