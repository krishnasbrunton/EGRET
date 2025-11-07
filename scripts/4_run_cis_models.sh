#!/bin/bash
#SBATCH --job-name="Running Fusion"
#SBATCH --output="batch_submissions/cis_fusion.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00


tissue=$1
chunk_size=$2
end_index=$3

if [ -n "$tissue" ]; then
    echo $tissue
else
    echo "no tissue defined"
    return
fi
home_dir=/expanse/lustre/projects/ddp412/kbrunton

#Define output directions
wd=working/$tissue/cis
mkdir -p $wd
tmpdir=tmp/$tissue/cis
mkdir -p $tmpdir
weights=FUSION/$tissue/cis
mkdir -p $weights

for gene in $(ls plink_results/${tissue}/cis | grep .bed | head -n $end_index | tail -n $chunk_size)
do
	gene=$(basename $gene .bed)
	gefile=${home_dir}/full_workflow/expression_files/${tissue}_expression.txt.gz
	individuals=../TWAS_across_tissues/individuals_per_tissue/${tissue}_individuals.txt
	covar=../../tamariutabartell/o2_files/gtex_covars/Covar_all_${tissue}.txt
	#covar=covariate_files/${tissue}_covariates.txt.gz

	alldonors=$(zcat $gefile | head -n 1)
	colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $individuals | awk 'BEGIN {ORS=","} {print $1}' | sed 's/,$//')

	../../kakamatsu/plink2 --bfile $home_dir/full_workflow/plink_results/${tissue}/cis/$gene --make-bed --keep $individuals --indiv-sort f $individuals --out $wd/${gene}
	rm $wd/${gene}.log

	rowid=$(zcat $gefile | nl | grep ${gene} | awk '{print $1}')
	ge_EURdonors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind)
	paste --delimiters='\t' <(cut -f1-5 $wd/${gene}.fam) <(echo $ge_EURdonors | sed 's/ /\n/g') > $wd/${gene}.mod.fam
	mv $wd/${gene}.mod.fam $wd/${gene}.fam

	TMP=$tmpdir/${gene}
	OUT=$weights/${gene}

	Rscript FUSION.compute_weights.R --bfile $wd/${gene} --crossval 5 --models lasso,blup,enet,top1 --hsq_p 1 --tmp $TMP --out $OUT --covar $covar --PATH_gcta ../../kakamatsu/fusion_twas-master/gcta_nr_robust --PATH_plink ../../kakamatsu/plink --PATH_gemma ../gemma-0.98.5-linux-static-AMD64

done
