#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="batch_submissions/Rscript.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 04:00:00



Rscript 9.0_run_WGCNA_clustering.R \
	--expression expression_files/Whole_Blood_expression.txt.gz \
	--out modules_5_PCs/ \
	--individuals ../TWAS_across_tissues/individuals_per_tissue/Whole_Blood_individuals.txt \
	--covariates covariate_files/Whole_Blood_covariates.txt.gz \
	--tissue Whole_Blood \
	--num_PCs 5

#Rscript 9.02_WGCNA_clustering_reimplemented.R \
#        --expression expression_files/Whole_Blood_expression.txt.gz \
#        --out transPCO/Whole_Blood/modules/ \
#        --individuals ../TWAS_across_tissues/individuals_per_tissue/Whole_Blood_individuals.txt \
#        --covariates covariate_files/Whole_Blood_covariates.txt.gz \
#        --tissue Whole_Blood


