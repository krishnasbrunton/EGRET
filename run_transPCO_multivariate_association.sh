#!/bin/bash

for i in {0..5}; do
	sbatch run_job.sh 9.2_transPCO_multivariate_association.R --tissue Whole_Blood --fold $i
	#Rscript 9.2_transPCO_multivariate_association.R --tissue Whole_Blood --fold $i
done
