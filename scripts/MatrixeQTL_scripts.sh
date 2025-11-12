tissue=$1        # tissue for which the analysis will be done on
FDR=$2           # FDR threshold at which to select variants
output_dir=$3    # output directory where results will be stored

Rscript 7_run_MatrixeQTL.R \
    --tissue ${tissue}

Rscript 8_make_MatrixeQTL_bed_by_FDR.R \
    --tissue ${tissue} --FDR ${FDR}


