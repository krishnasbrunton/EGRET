#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="batch_submissions/Rscript.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00


tissues=$(ls ../TWAS_across_tissues/weights)

for tissue in $tissues
do
	sbatch run_job.sh \
		3_setup_cis_model.R \
		--gene_information ../data/GTEx_V8.txt.gz \
		--tissue ${tissue} \
		--expression_file expression_files/${tissue}_expression_regressed.txt.gz \
		--plink_path ../../kakamatsu/plink2 \
		--bfile genotype_files/GTEX_v8_genotypes_pruned


done

