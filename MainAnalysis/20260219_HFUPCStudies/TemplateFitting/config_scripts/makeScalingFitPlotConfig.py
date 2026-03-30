import json

SampleFileTag = "noiseFit_MC_to_data"
SampleFileName = f"TemplateFitting/SampleConfigs/ScalingFitPlot_{SampleFileTag}.json"

OutputFileDir = "TemplateFitting/output_ScalingFit"

TargetFile = "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root"
FileToFit = "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root"

OutputFilesArray = [
    {
        "TargetFile": TargetFile,
        "FileToFit": FileToFit,
        "VarFitName": "HFEMaxPlus_forest",
        "VarTargetName": "HFEMaxPlus_forest",
        "OutputFileName": f"{OutputFileDir}/ScalingFitResult_EMC_EBX_plus.root",
        "TreeName": "OutputTree",
        "XMin": 0,
        "XMax": 50,
        "NBins": 200
    },
    {
        "TargetFile": TargetFile,
        "FileToFit": FileToFit,
        "VarFitName": "HFEMaxMinus_forest",
        "VarTargetName": "HFEMaxMinus_forest",
        "OutputFileName": f"{OutputFileDir}/ScalingFitResult_EMC_EBX_minus.root",
        "TreeName": "OutputTree",
        "XMin": 0,
        "XMax": 50,
        "NBins": 200
    }
]

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)
