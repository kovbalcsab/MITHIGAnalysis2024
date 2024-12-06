PTMIN=0
PTMAX=100000
YMIN=-2.4
YMAX=+2.4
TRIGGER=0 # no trigger selection
ISGAMMAN=0
ISDATA=0
DOINCLUSIVE=1
DOPROMPT=0
DONONPROMPT=0
INPUT=/data/NewSkims23_24/20241203_SkimPbPbMC_ForestForcedD0Decay100M_BeamAPtHat0.root
OUTPUTANALYSIS=output.root
OUTPUTFIT=outputfit.root
./ExecuteDzeroUPC --Input $INPUT --MinDzeroPT $PTMIN --MaxDzeroPT $PTMAX --MinDzeroY $YMIN --MaxDzeroY $YMAX --IsGammaN $ISGAMMAN --TriggerChoice $TRIGGER --IsData $ISDATA --Output $OUTPUTANALYSIS --DoInclusive $DOINCLUSIVE --DoPrompt $DOPROMPT --DoNonPrompt $DONONPROMPT
wait
