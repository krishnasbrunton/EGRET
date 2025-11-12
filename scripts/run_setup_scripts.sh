#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="TEST_SETUP.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 08:00:00

source ~./bashrc
conda activate transTWAS_env

tissue="Whole_Blood"

genotypes_file_path="../../../tamariutabartell/gtex/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder"
expression_file_path="../../../tamariutabartell/gtex/${tissue}.v8.normalized_expression.bed.gz"
individuals_file_path="../../TWAS_across_tissues/individuals_per_tissue/${tissue}_individuals.txt"
covariates_file_path="../../../tamariutabartell/o2_files/gtex_covars/Covar_all_${tissue}.txt"
LD_prune_r2="0.9"
plink_path="../plink2"
genotype_output_prefix="GTEX_v8_genotypes_pruned"
folds="5"
gene_info_file_path="../../data/GTEx_V8.txt.gz"
output_dir="test_output"

source setup_genotype_and_expression.sh \
    $genotypes_file_path \
    $expression_file_path \
    $covariates_file_path \
    $tissue \
    $individuals_file_path \
    $LD_prune_r2 \
    $plink_path \
    $genotype_output_prefix \
    $folds \
    $gene_info_file_path \
    $output_dir
