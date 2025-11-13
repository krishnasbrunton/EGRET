#!/bin/bash

tissue=$1                   # tissue for which the analysis will be done on
FDR=$2                      # FDR threshold at which to select variants
folds=$3                    # number of folds
gene_info=$4                # gene info file containing gene id, chr, start location
output_dir=$6               # output directory where results will be stored


Rscript 5_run_GBAT_updated.R \
    --tissue $tissue \
    --FDR $FDR \
    --gene_info $gene_info \
    --output_dir $output_dir

Rscript 5.1_run_GBAT_association.R \
    --tissue $tissue \
    --output_dir $output_dir

Rscript 5.3_make_GBAT_bed_by_FDR.R \
    --tissue $tissue \
    --FDR $FDR \
    --output_dir $output_dir