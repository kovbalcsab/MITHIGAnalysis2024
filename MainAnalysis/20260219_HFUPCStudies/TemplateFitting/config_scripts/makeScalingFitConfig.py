import json

SampleFileTag = "noiseFit_MC_to_data"
SampleFileName = f"TemplateFitting/SampleConfigs/ScalingFit_{SampleFileTag}.json"

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
        "aBase": 2.0,
        "bBase": 0.0,
        "BinsPerGeV": 4,
        "UseKeysPdf": True,
        "FitNumCPU": 4,
        "FitStrategy": 2,
        "FitOffset": True,
        "ClipScaledE": True,
        "UseAbsJacobian": True,
        "AMin": 0.0,
        "AMax": 5.0,
        "BMin": 0,
        "BMax": 10.0
    },
    {
        "TargetFile": TargetFile,
        "FileToFit": FileToFit,
        "VarFitName": "HFEMaxMinus_forest",
        "VarTargetName": "HFEMaxMinus_forest",
        "OutputFileName": f"{OutputFileDir}/ScalingFitResult_EMC_EBX_minus.root",
        "aBase": 2.0,
        "bBase": 0.0,
        "BinsPerGeV": 4,
        "UseKeysPdf": True,
        "FitNumCPU": 4,
        "FitStrategy": 2,
        "FitOffset": True,
        "ClipScaledE": True,
        "UseAbsJacobian": True,
        "AMin": 0.0,
        "AMax": 5.0,
        "BMin": 0,
        "BMax": 10.0
    }
]

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)

