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

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
pval_threshold=0.000001

Rscript 7_run_MatrixeQTL.R \
    --tissue ${tissue}

Rscript 8_make_MatrixeQTL_bed.R \
    --tissue ${tissue} --pval_threshold ${pval_threshold}


