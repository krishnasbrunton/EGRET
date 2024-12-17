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

Rscript $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}
