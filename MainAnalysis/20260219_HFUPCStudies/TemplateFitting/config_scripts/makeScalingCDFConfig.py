import json

SampleFileTag = "noise_MC_to_data"
SampleFileName = f"TemplateFitting/SampleConfigs/ScalingCDF_{SampleFileTag}.json"

OutputFileDir = "TemplateFitting/output_ScalingFit"

TargetFile = "EmptyBX/output_EmptyBX/emptyBXOutput_20260226_Trigger0_nTrkFilter1.root"
FileToFit = "EmptyMC/output_EmptyMC/emptyMCOutput_20260320.root"

base = {
    "TargetFile": TargetFile,
    "FileToFit": FileToFit,
    "NBinsTarget": 1000000,
    "NBinsFit": 1000000,
    "EMinTarget": 1,
    "EMaxTargetROI": 200,
    "EMinFit": 1,
    "EMaxFitROI": 200,
    "MapFitMin": 1,
    "MapFitMax": 6.8,
    "PolyOrder": 2,
    "ForceZeroConstant": 1,
    "ScalingPlotYMax": 50,
    "ScalingPlotXMax": 15,
    "PlotNBinsTarget": 200,
    "PlotNBinsFit": 200,
    "PlotEMaxTarget": 100,
    "PlotEMaxFit": 100,
    "PlotYMax": 1,
    "PlotXMax": 10
}

OutputFilesArray = [
    {
        **base,
        "VarFitName": "HFEMaxPlus_forest",
        "VarTargetName": "HFEMaxPlus_forest",
        "ScalingOutputFileName": f"{OutputFileDir}/ScalingCDFResult_EMC_EBX_plus.root",
        "PlotOutputFileName": f"{OutputFileDir}/ScalingCDFResultPlot_EMC_EBX_plus.pdf"
    },
    {
        **base,
        "VarFitName": "HFEMaxMinus_forest",
        "VarTargetName": "HFEMaxMinus_forest",
        "ScalingOutputFileName": f"{OutputFileDir}/ScalingCDFResult_EMC_EBX_minus.root",
        "PlotOutputFileName": f"{OutputFileDir}/ScalingCDFResultPlot_EMC_EBX_minus.pdf"
    }
]

data = {"OutputFileDir": OutputFileDir, "OutputFiles": OutputFilesArray}

with open(SampleFileName, "w") as f:
    json.dump(data, f, indent=4)

