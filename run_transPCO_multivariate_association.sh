#!/bin/bash

#for i in {0..5}; do
for i in {0..5}; do
	#sbatch run_job.sh 9.2_transPCO_multivariate_association.R --tissue Whole_Blood --fold $i --module_dir modules_old --results_dir PCO_association_results_old --association_dir association_results_msigdb
	#Rscript 9.2_transPCO_multivariate_association.R --tissue Whole_Blood --fold $i
	sbatch run_job.sh 9.2_transPCO_multivariate_association.R --tissue Whole_Blood --fold $i --module_dir modules_60_PCs_power_3 --results_dir PCO_association_results_60_PCs_power_3 --association_dir association_results_60_PCs_power_3

done
