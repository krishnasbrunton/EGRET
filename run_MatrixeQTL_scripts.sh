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


#Rscript 7_run_MatrixeQTL.R --tissue Whole_Blood

#Rscript 12.1_make_xtune_bed_cis_MatrixeQTL.R --tissue Whole_Blood
