import json

PtYBins = [ [2,5,-2,-1],
            [2,5,-1,0],
            [2,5,0,1],
            [2,5,1,2],
            [2,5,-2,2]  ]

SampleFileTag="test"
SampleFileName = f"UPCDataSignal/SampleConfigs/{SampleFileTag}.json"

Data="/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_PbPbUPC_HIForward0_Drej-pasor.root"

OutputFileDir=f"UPCDataSignal/output_UPCDataSignal/{SampleFileTag}"

OutputFilesArray=[]
DoSystD=0
TriggerChoice=1

for i, (PtMin, PtMax, YMin, YMax) in enumerate(PtYBins):
    OutputFilesArray.append( { "FileName": "Data.root",
                            "Input": Data,
                            "DzeroPTMin": PtMin,
                            "DzeroPTMax": PtMax,
                            "DzeroYMin": YMin,
                            "DzeroYMax": YMax,
                            "DoSystD": DoSystD,
                            "TriggerChoice": TriggerChoice,
                            "IsGammaN": 1 } )
    OutputFilesArray.append( { "FileName": "Data.root",
                            "Input": Data,
                            "DzeroPTMin": PtMin,
                            "DzeroPTMax": PtMax,
                            "DzeroYMin": YMin,
                            "DzeroYMax": YMax,
                            "DoSystD": DoSystD,
                            "TriggerChoice": TriggerChoice,
                            "IsGammaN": 0 } )


data={ "OutputFileDir": OutputFileDir,
       "OutputFiles": OutputFilesArray}

with open(SampleFileName, 'w') as f:
    json.dump(data, f, indent=4)
