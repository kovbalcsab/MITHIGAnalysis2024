import json

SampleFileTag = "validation_quater0"
SampleFileName = f"TemplateFitting/SampleConfigs/UnfoldSVDValidation_{SampleFileTag}.json"

OutputFileDir = f"TemplateFitting/output_unfolding/SVD_{SampleFileTag}"

base = {
    "MeasuredHist": "hMeasured",
    "RefoldedPrefix": "hRefolded_iter",
    "UnfoldedPrefix": "hUnfolded_iter",
    "RegularizationTree": "RegularizationTree",
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
            "InputFile": f"{OutputFileDir}/unfolding_HFEMaxMinus_forest_SVD.root",
            "OutputFile": f"{OutputFileDir}/unfolding_HFEMaxMinus_forest_SVD_validation_vs_q{validation_quarter}.root",
            "ValidationDataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN0/Data.root",
            "ValidationVarName": "HFEMaxMinus_forest"
        },
        {
            **base,
            "ValidationQuarter": validation_quarter,
            "InputFile": f"{OutputFileDir}/unfoldingMC_HFEMaxPlus_forest_SVD.root",
            "OutputFile": f"{OutputFileDir}/unfoldingMC_HFEMaxPlus_forest_SVD_validation_vs_q{validation_quarter}.root",
            "ValidationDataFile": "MCSignal/output_MCSignal/test/pt2-5_y-2-2_IsGammaN1/SignalMC.root",
            "ValidationVarName": "HFEMaxPlus_forest"
        },
        {
            **base,
            "ValidationQuarter": validation_quarter,
            "InputFile": f"{OutputFileDir}/unfolding_HFEMaxPlus_forest_SVD.root",
            "OutputFile": f"{OutputFileDir}/unfolding_HFEMaxPlus_forest_SVD_validation_vs_q{validation_quarter}.root",
            "ValidationDataFile": "UPCDataSignal/output_UPCDataSignal/test/pt2-5_y-2-2_IsGammaN1/Data.root",
            "ValidationVarName": "HFEMaxPlus_forest"
        }
    ])

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)
