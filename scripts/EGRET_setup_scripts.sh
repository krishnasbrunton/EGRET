#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="TEST_EGRET_SETUP_CROSSMAP.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 08:00:00

home_dir="/expanse/lustre/projects/ddp412/kbrunton/EGRET"
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
FDR="0.1"
num_PCs="10"
models="lasso,enet,blup,xtune"
output_dir="test_output/"

#additional file paths for setup scripts
crossmap_file="$(dirname $home_dir)/mappability_encode/cross_mappability_strength.txt.gz"

Rscript 12.0_make_crossmappable_gene_beds.R \
    --gene_info $gene_info_file_path \
    --crossmap_file $crossmap_file \
    --output_dir $output_dir