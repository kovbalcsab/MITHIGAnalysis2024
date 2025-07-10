./ExecuteSystematicPlotter --Input "mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR14_12.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR14_6.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR14_16.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR12_12.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR16_12.root" \
                            --Output HFSyst_HFANDOFF_HFORON \
                            --SystLabels "HFAND12 (ZB + onlineHFOR14),HFAND6 (ZB + onlineHFOR14),HFAND16 (ZB + onlineHFOR14),HFAND12 (ZB + onlineHFOR12),HFAND12 (ZB + onlineHFOR16)" \
                            --ZBPrescale 0,0,0,0,0 \
                            --DoRatio true

./ExecuteSystematicPlotter --Input "mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFAND14_HFAND12OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFAND14_HFAND6OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFAND14_HFAND16OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFAND12_HFAND12OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFAND16_HFAND12OFF.root" \
                            --Output HFSyst_HFANDOFF_HFANDON \
                            --SystLabels "HFAND12 (ZB + onlineHFAND14),HFAND6 (ZB + onlineHFAND14),HFAND16 (ZB + onlineHFAND14),HFAND12 (ZB + onlineHFAND12),HFAND12 (ZB + onlineHFAND16)" \
                            --ZBPrescale 0,0,0,0,0 \
                            --DoRatio true

./ExecuteSystematicPlotter --Input "mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR14_HFOR12OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR14_HFOR6OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR14_HFOR16OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR12_HFOR12OFF.root,mergedOutputs/mergedOutput_TRIG1_TRKPT1_ONLINEHFOR16_HFOR12OFF.root" \
                            --Output HFSyst_HFOROFF_HFORON \
                            --SystLabels "HFOR12 (ZB + onlineHFOR14),HFOR6 (ZB + onlineHFOR14),HFOR16 (ZB + onlineHFOR14),HFOR12 (ZB + onlineHFOR12),HFOR12 (ZB + onlineHFOR16)" \
                            --ZBPrescale 0,0,0,0,0 \
                            --DoRatio true