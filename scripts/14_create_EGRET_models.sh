#!/bin/bash
#SBATCH --job-name="Fusion"
#SBATCH --output="batch_submissions/gtex_xtune.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00


tissue=$1
if [ -n "$tissue" ]; then
    echo $tissue
else
    echo "no tissue defined"
    return
fi
home_dir=/expanse/lustre/projects/ddp412/kbrunton

chunk_size=$2
#Define output directions
weights=xtune_fusion_models/${tissue}/cis_MatrixeQTL_GBAT_transPCO_4_mismatches_no_cis
mkdir -p $weights
wd=working/${tissue}/cis_MatrixeQTL_GBAT_transPCO_4_mismatches_no_cis
mkdir -p $wd
output=xtune_fusion_results/${tissue}/cis_MatrixeQTL_GBAT_transPCO_4_mismatches_no_cis
mkdir -p $output

z_matrix_dir=z_matrices/${tissue}/cis_MatrixeQTL_GBAT_transPCO_4_mismatches_no_cis

end_index=$3
plink_dir=plink_results/$tissue/cis_MatrixeQTL_GBAT_transPCO_4_mismatches

for gene in $(ls ${plink_dir}/fold_0/ | grep .bed | head -n $end_index | tail -n $chunk_size)
#for gene in $(find transPCO/Whole_Blood/fold_1/bed_files/ transPCO/Whole_Blood/fold_2/bed_files/ transPCO/Whole_Blood/fold_3/bed_files/ transPCO/Whole_Blood/fold_4/bed_files/ transPCO/Whole_Blood/fold_5/bed_files/ -type f -exec basename {} \; | sort | uniq)
do
gene=${gene::-4}

#Extracting gene expression information
gefile=expression_files/${tissue}_expression.txt.gz
individuals=fold_0_info/${tissue}/train_individuals.txt

alldonors=$(zcat $gefile | head -n 1)
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $individuals | awk 'BEGIN {ORS=","} {print $1}') #donor list always the same
colind2=${colind%,}

for i in {0..5} #changed
do
	mkdir -p $wd/fold_$i
	../../kakamatsu/plink2 --bfile $plink_dir/fold_$i/$gene --make-bed --keep $individuals --indiv-sort f $individuals --out $wd/fold_$i/${gene}
	rm $wd/fold_$i/${gene}.log

	rowid=$(zcat $gefile | awk 'NR > 1 {print $1}' | nl | grep ${gene} | awk '{print $1 + 1}')
	ge_EURdonors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)
	paste --delimiters='\t' <(cut -f1-5 $wd/fold_$i/${gene}.fam) <(echo $ge_EURdonors | sed 's/ /\n/g') > $wd/fold_$i/${gene}.mod.fam
	mv $wd/fold_$i/${gene}.mod.fam $wd/fold_$i/${gene}.fam

done
echo "running xtune fusion"
Rscript xtune_fusion_cis_trans.R --gene $gene --working_dir $wd --models xtune,lasso,enet,blup --output_dir $output --weights_dir $weights --tissue $tissue --PATH_plink ../../kakamatsu/plink2 --PATH_gemma ../full_workflow/gemma-0.98.5-linux-static-AMD64 --covar covariate_files/${tissue}_covariates.txt.gz --z_matrix_dir $z_matrix_dir
done
