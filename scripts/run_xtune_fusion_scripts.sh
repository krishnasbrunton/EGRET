#!/bin/bash

tissue=Whole_Blood
plink_dir=plink_results/Whole_Blood/cis_MatrixeQTL_GBAT_transPCO_4_mismatches/fold_0
files=($(ls "$plink_dir" | grep bed))
echo ${#files[@]}
# Define the chunk size
chunk_size=500


num_chunks=$(((${#files[@]} + $chunk_size - 1) / $chunk_size))
echo $num_chunks
# Loop through each chunk
for ((i = 0; i < num_chunks; i++)); do
        start_index=$((i * chunk_size))
        end_index=$((start_index + chunk_size - 1))
	sbatch 14_run_fusion.sh $tissue $chunk_size $end_index
done

