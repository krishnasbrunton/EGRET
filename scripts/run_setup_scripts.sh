#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="batch_submissions/Rscript.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 08:00:00


Rscript 0_setup_genotypes.R \
	--plink_path ../../kakamatsu/plink2 \
        --bfile ../../tamariutabartell/gtex/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder \
        --LD_r2 0.9 \
        --LD_chunk 1 \
        --LD_window 100 \
        --out genotype_files/GTEX_v8_genotypes_pruned


tissues=$(ls ../TWAS_across_tissues/weights)

for tissue in $tissues
do 

	Rscript 1_setup_expression.R \
		--expression ../../tamariutabartell/gtex/${tissue}.v8.normalized_expression.bed.gz \
		--out expression_files/${tissue}_expression.bed \
		--individuals ../TWAS_across_tissues/individuals_per_tissue/${tissue}_individuals.txt \
		--covariates ../../tamariutabartell/o2_files/gtex_covars/Covar_all_${tissue}.txt \
		--tissue ${tissue}

	Rscript 2_setup_folds.R \
		--expression expression_files/${tissue}_expression_regressed.txt.gz \
		--folds 5 \
		--tissue ${tissue} \
		--covariates covariate_files/${tissue}_covariates.txt.gz \
		--individuals ../TWAS_across_tissues/individuals_per_tissue/${tissue}_individuals.txt \
		--plink_path ../../kakamatsu/plink2 \
		--bfile genotype_files/GTEX_v8_genotypes_pruned

done
