import json

PtYBins = [ [2,5,-2,-1],
            [2,5,-1,0],
            [2,5,0,1],
            [2,5,1,2],
            [2,5,-2,2]  ]

SampleFileTag="test"
SampleFileName = f"MCSignal/SampleConfigs/{SampleFileTag}.json"

MCBeamAPath="/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260212-hfle_HiForest_260218_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_t2.root"
MCBeamBPath=""

OutputFileDir=f"MCSignal/output_MCSignal/{SampleFileTag}"

OutputFilesArray=[]
useClusterCompatibilityFilter=0
useGenDSignalEventSelection=1
useRecoDSignalEventSelection=0
DoSystD=0
isResolved=0

MakeBeamA=True
MakeBeamB=False

if MakeBeamA:
    for i, (PtMin, PtMax, YMin, YMax) in enumerate(PtYBins):
        OutputFilesArray.append( { "FileName": "SignalMC.root",
                                "Input": MCBeamAPath,
                                "DzeroPTMin": PtMin,
                                "DzeroPTMax": PtMax,
                                "DzeroYMin": YMin,
                                "DzeroYMax": YMax,
                                "useClusterCompatibilityFilter": useClusterCompatibilityFilter,
                                "useGenDSignalEventSelection": useGenDSignalEventSelection,
                                "useRecoDSignalEventSelection": useRecoDSignalEventSelection,
                                "DoSystD": DoSystD,
                                "isResolved": isResolved,
                                "IsGammaN": 1 } )

if MakeBeamB:
    for i, (PtMin, PtMax, YMin, YMax) in enumerate(PtYBins):
        OutputFilesArray.append( { "FileName": "SignalMC.root",
                                "Input": MCBeamBPath,
                                "DzeroPTMin": PtMin,
                                "DzeroPTMax": PtMax,
                                "DzeroYMin": YMin,
                                "DzeroYMax": YMax,
                                "useClusterCompatibilityFilter": useClusterCompatibilityFilter,
                                "useGenDSignalEventSelection": useGenDSignalEventSelection,
                                "useRecoDSignalEventSelection": useRecoDSignalEventSelection,
                                "DoSystD": DoSystD,
                                "isResolved": isResolved,
                                "IsGammaN": 0 } )

data={ "OutputFileDir": OutputFileDir,
       "OutputFiles": OutputFilesArray}

with open(SampleFileName, 'w') as f:
    json.dump(data, f, indent=4)
