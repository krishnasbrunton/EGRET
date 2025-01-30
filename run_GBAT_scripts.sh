#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="batch_submissions/Rscript.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 12:00:00

tissue=$1
Rscript 5_run_GBAT.R --tissue $tissue

Rscript 5.1_run_GBAT_association.R --tissue $tissue

Rscript 5.2_make_GBAT_bed.R --tissue $tissue
