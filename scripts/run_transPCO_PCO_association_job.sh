#!/bin/bash
#SBATCH --job-name="transPCO_PCO_assoc"
#SBATCH --output="batch_submissions/transPCO_PCO_assoc.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 04:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
module=$2
module_dir=$3
output_dir=$4
folds=$5
gene_info=$6

Rscript 9.3_transPCO_multivariate_association_by_module.R \
    --tissue $tissue \
    --module $module \
    --module_dir $module_dir \
    --output_dir $output_dir \
    --folds $folds \
    --gene_info $gene_info \
    --results_dir PCO_association_results \
    --association_dir association_results
