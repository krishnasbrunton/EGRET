#!/bin/bash
tissue=$1
folds=$2
models=$3
output_dir=$4

Rscript 12.4_make_xtune_bed_no_cross_mappable_by_fold.R \
    --tissue $tissue \
    --MatrixeQTL_bed_dir results_FDR_0.1 \
    --GBAT_bed_dir results_FDR_0.1 \
    --transPCO_bed_dir bed_files_FDR_0.1 \
    --output_dir cis_transPCO_FDR_0.1 \
    --exclude_crossmap TRUE \
    --models Cis,MatrixeQTL,GBAT,transPCO \
    --output_dir $output_dir 

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