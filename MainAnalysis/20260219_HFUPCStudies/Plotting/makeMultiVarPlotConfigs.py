import json

PtYBins = [ [2,5,-2,-1],
            [2,5,-1,0],
            [2,5,0,1],
            [2,5,1,2],
            [2,5,-2,2]  ]

IsGammaNArray = [1]

ConfigFileTag="test"
PlottingDir=f"Plotting/MultiVarPlots/{ConfigFileTag}"
PlotConfigFileName = f"Plotting/multiVarPlotConfigs/PlotConfig_{ConfigFileTag}.json"

# Input files
InputArray = [ {"SampleDir": "MCSignal/output_MCSignal/test", "FileNames": "SignalMC.root", "SampleLabel": "test", "AddKinematicRange": 1},
               {"SampleDir": "EmptyBX/output_EmptyBX", "FileNames": "emptyBXOutput_20260226_Trigger-1_nTrkFilter0.root", "SampleLabel": "EmptyBX", "AddKinematicRange": 0},  
               {"SampleDir": "EmptyZB/output_EmptyZB", "FileNames": "emptyZBOutput_20260226_Trigger-1_nTrkFilter1_nVtxFilter0_ZDCFilter1.root", "SampleLabel": "EmptyZB", "AddKinematicRange": 0}  ]

# Histogram list
HistogramsArray = [ {"HistNameBase": "HFEMaxPlus_forest", "XTitle": "HF~E_{max}^{+}~[GeV]", "YTitle": "Normalized~Entries", "LogY": 1, "XLims": [0, 80], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 1, "IsTreeVar": 1, "NBins": 80, "doRatio": 1},
                    {"HistNameBase": "HFEMaxMinus_forest", "XTitle": "HF~E_{max}^{-}~[GeV]", "YTitle": "Normalized~Entries", "LogY": 1, "XLims": [0, 80], "YLims": [1e-4, 1.2], "Normalize": 1, "RebinFactor": 1, "IsTreeVar": 1, "NBins": 80, "doRatio": 1} ]

PlotArray = []
for IsGammaN in IsGammaNArray:
    for (PtMin, PtMax, YMin, YMax) in PtYBins:
        for HistConfig in HistogramsArray:
            CurrInputArray = []
            CurrLabelArray = []
            for Input in InputArray:
                CurrLabelArray.append(Input["SampleLabel"])
                if (Input["AddKinematicRange"]):
                    CurrInputArray.append(f"{Input['SampleDir']}/pt{PtMin}-{PtMax}_y{YMin}-{YMax}_IsGammaN{IsGammaN}/{Input['FileNames']}")
                else:
                    CurrInputArray.append(f"{Input['SampleDir']}/{Input['FileNames']}")
            PlotArray.append( { "PlotBaseName": f"{HistConfig['HistNameBase']}.pdf",
                                "InputFileNames": ','.join(CurrInputArray),
                                "VarNames": ','.join([HistConfig["HistNameBase"]]*len(CurrInputArray)),
                                "Labels": ','.join(CurrLabelArray),
                                "xTitle": HistConfig["XTitle"],
                                "yTitle": HistConfig["YTitle"],
                                "logY": HistConfig["LogY"],
                                "xMin": HistConfig["XLims"][0],
                                "xMax": HistConfig["XLims"][1],
                                "yMin": HistConfig["YLims"][0],
                                "yMax": HistConfig["YLims"][1],
                                "normalizeToUnity": HistConfig["Normalize"],
                                "rebinFactor": HistConfig["RebinFactor"],
                                "IsTreeVar": ','.join([str(HistConfig["IsTreeVar"])]*len(CurrInputArray)),
                                "NBins": ','.join([str(HistConfig["NBins"])]*len(CurrInputArray)),
                                "doRatio": HistConfig["doRatio"],
                                "MinDzeroPT": PtMin,
                                "MaxDzeroPT": PtMax,
                                "MinDzeroY": YMin,
                                "MaxDzeroY": YMax,
                                "IsGammaN": IsGammaN } )

data={ "plottingDir": PlottingDir,
        "Plots": PlotArray}

with open(PlotConfigFileName, 'w') as json_file:
    json.dump(data, json_file, indent=4)