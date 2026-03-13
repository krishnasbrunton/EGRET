#!/bin/bash
#SBATCH --job-name="transPCO_results"
#SBATCH --output="batch_submissions/transPCO_results.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 02:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
pval_threshold=$2
output_dir=$3
folds=$4

Rscript 11_transPCO_results_analysis.R \
    --tissue $tissue \
    --pval_threshold $pval_threshold \
    --output_dir $output_dir \
    --folds $folds \
    --PCO_association_dir PCO_association_results \
    --module_dir modules
