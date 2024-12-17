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




#tissues=$(ls ../TWAS_across_tissues/weights)
tissues=("Cells_Cultured_fibroblasts" "Skin_Sun_Exposed_Lower_leg" "Muscle_Skeletal")

for tissue in ${tissues[@]}   # I changed this so that it can use a bash array. If using the ls command, you can remove the @ symbol
do
	echo $tissue
	plink_dir=plink_results/${tissue}/cis/
	files=($(ls "$plink_dir" | grep bed))
	echo ${#files[@]}
	# Define the chunk size
	chunk_size=2000

	# Calculate the number of chunks
	num_chunks=$(((${#files[@]} + $chunk_size - 1) / $chunk_size))
	echo $num_chunks
# Loop through each chunk
	for ((i = 0; i < num_chunks; i++)); do
        	start_index=$((i * chunk_size))
        	end_index=$((start_index + chunk_size - 1))
		sbatch 4_run_cis_models.sh $tissue $chunk_size $end_index

	done
done
