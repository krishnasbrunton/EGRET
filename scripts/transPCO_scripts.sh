#!/bin/bash
tissue=$1
FDR=$3
folds=$2
num_PCs=$2
output_dir=$4


Rscript 9.0_run_WGCNA_clustering.R \
	--expression expression_files/Whole_Blood_expression.txt.gz \
	--out modules/ \
	--individuals ../TWAS_across_tissues/individuals_per_tissue/Whole_Blood_individuals.txt \
	--covariates covariate_files/Whole_Blood_covariates.txt.gz \
	--tissue $tissue \
	--num_PCs 5

Rscript 9.1_initialize_chr_genotypes.R \
    --tissue $tissue

directories=("transPCO/$tissue/fold_0/modules" "transPCO/$tissue/fold_1/modules" "transPCO/$tissue/fold_2/modules" "transPCO/$tissue/fold_3/modules" "transPCO/$tissue/fold_4/modules" "transPCO/$tissue/fold_5/modules")

# Initialize a variable to store the max number of files
max_files=0

# Loop through each directory
for dir in "${directories[@]}"; do
    # Count the number of files in the directory (not including subdirectories)
    file_count=$(find "$dir" -type f | wc -l)
        echo $file_count   
    # Check if the current directory has more files than the current max
    if (( file_count > max_files )); then
        max_files=$file_count
        max_module_dir=$dir
    fi
done
echo $max_module_dir
module_path=$max_module_dir/

for file in "$module_path"*; do
    Rscript 9.10_cleaned_up_MatrixeQTL_for_transPCO.R \
        --tissue $tissue \
        --module $(basename $file) \
        --module_dir modules \
        --output_dir association_results
done

modules=$(find transPCO/${tissue}/fold_{0..5}/modules/ -type f -exec basename {} \; | sort | uniq)

for module in $modules
do
    echo $module
    Rscript 9.3_transPCO_multivariate_association_by_module.R \
    --tissue $tissue \
    --module_dir modules \
    --module ${module} \
    --results_dir PCO_association_results \
    --association_dir association_results
done

Rscript 11_transPCO_results_analysis_by_FDR.R \
    --tissue $tissue \
    --FDR $FDR \
    --PCO_association_dir PCO_association_results \
    --output_dir bed_files_FDR_$fdr \
    --module_dir modules