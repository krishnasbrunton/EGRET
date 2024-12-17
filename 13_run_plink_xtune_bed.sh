#!/bin/bash
#SBATCH --job-name="Plink"
#SBATCH --output="batch_submissions/plinkResults.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 12:00:00

tissue=$1

#for i in {0..5}
#for i in {5} #changed for just fold since many tissues did not finish
#do
	bed_dir=xtune_bed_files/$tissue/cis_MatrixeQTL_1e-06/fold_5
	output_dir=plink_results/$tissue/cis_MatrixeQTL_1e-06/fold_5
	mkdir -p $output_dir

	for gene in $(ls $bed_dir)
	do
		#geneOnly=${gene::-4}  #this was commented out because when I created the bed files I did not put the .bed ending on the files
		../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned_new --extract bed1 ${bed_dir}/$gene --out ${output_dir}/$gene --make-bed
	rm ${output_dir}/$gene.log
	done
#done
