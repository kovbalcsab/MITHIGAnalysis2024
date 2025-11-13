#!/bin/bash

FOLDER_IDX=(0 1 2 3 4 5)
FILE_IDX=(1 2 3 4 5 6 7 8 9)
# Merging ROOT files using hadd

for folder in "${FOLDER_IDX[@]}"; do
    for file in "${FILE_IDX[@]}"; do
        input_files="/data00/bakovacs/pOsamples/Skims/output_202510*_Skim_v1_pO_sample_000${folder}/output_${file}*.root"
        output_file="/data00/bakovacs/pOsamples/Skims/merged_skims/merged_output_${folder}_${file}.root"
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping merge."
            continue
        fi
        
        # Add output_0.root to 1st file index
        if [ "$file" -eq 1 ]; then
            input_files="/data00/bakovacs/pOsamples/Skims/output_202510*_Skim_v1_pO_sample_000${folder}/output_0.root $input_files"
        fi
        hadd $output_file $input_files
    done
done
echo "Merging completed."
