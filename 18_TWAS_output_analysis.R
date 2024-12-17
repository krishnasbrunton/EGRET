library('data.table')


#TWAS_results_dir <- "TWAS_association_results/Whole_Blood/cis_MatrixeQTL/"
TWAS_results_dir = "TWAS_association_results/Whole_Blood/cis_MatrixeQTL_GBAT_transPCO/"
all_files <- list.files(TWAS_results_dir)
files <- all_files[!grepl("MHC", all_files)]

print(files)
total_significant = 0
for (file in files) {
	results = fread(paste0(TWAS_results_dir,file),header = T)
	results$TWAS.P[which(is.na(results$TWAS.P))] = 1
	results$TWAS.Z[which(is.na(results$TWAS.Z))] = 0
	adjusted_pvals = p.adjust(results$TWAS.P,"bonferroni")
	print(sum(adjusted_pvals < 0.05))
	print(mean(abs(results$TWAS.Z)))
	total_significant = total_significant + sum(adjusted_pvals < 0.05)
	#print(sum(results$TWAS.P < 0.0000001))

}

print(paste("The total number significant is ",total_significant))
