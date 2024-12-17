library('data.table')
library('scales')

#results_yaxis = fread("results_sumstats/Whole_Blood/cis_GBAT_results.txt",header = T)
#results_yaxis = fread("../basil/basil_results_summary_files/Whole_Blood_gtex_allsnps_results.txt",header = T)


#results_yaxis = fread("results_sumstats/Whole_Blood/cis_GBAT_1e-04_0_mismatches.txt",header = T)
#results_yaxis = fread("results_sumstats/Muscle_Skeletal/cis_MatrixeQTL_1e-06_results.txt",header = T)


#results_yaxis = fread("results_sumstats/Whole_Blood/cis_GBAT_BASIL_transPCO_eqtlgen_results.txt",header = T)
#results_yaxis = fread("results_sumstats/Whole_Blood/cis_GBAT_BASIL_lasso_results.txt",header = T)
#results_yaxis = fread("results_sumstats/Whole_Blood/cis_GBAT_BASIL_results.txt",header = T)
#results_xaxis = fread("results_sumstats/Whole_Blood/cis_BASIL_results.txt",header = T)
#results_yaxis = fread("../trans_adapt/results_sumstats/Whole_Blood/cis_GBAT_results.txt",header = T)
results_yaxis = fread("results_sumstats/Whole_Blood/cis_transPCO.txt")

#results_yaxis = fread("results_sumstats/Heart_Atrial_Appendage/cis_MatrixeQTL_1e-06_results.txt",header = T)

#fusion_tf_results = fread("../TWAS_across_tissues/fusion_results_summary/Whole_Blood/gtex_tf_chip_individual_tfs_results.txt",header = T)
#xtune_results_no_z = fread("results_sumstats/Whole_Blood/gtex_tf_chip_no_z_xtune_results.txt",header = T)

#default cis results
#results_xaxis = fread("../xtune/results_sumstats/Whole_Blood/gtex_cis_results.txt",header = T)
#results_xaxis = fread("../TWAS_across_tissues/fusion_results_summary/Heart_Atrial_Appendage/gtex_cis_results.txt",header = T)

#xtune_cis_results = fread("results_sumstats/Whole_Blood/cis_GBAT_BASIL_results.txt",header = T)
results_xaxis = fread("results_sumstats/Whole_Blood/cis_results.txt",header = T)
#results_xaxis = fread("../TWAS_across_tissues/fusion_results_summary/Whole_Blood/gtex_cis_results.txt",header = T)
#results_xaxis = fread("../trans_adapt_10_fold/results_sumstats/Whole_Blood/cis_results.txt",header = T)
#results_xaxis = fread("../trans_adapt_10_fold/results_sumstats/Whole_Blood/cis_GBAT_BASIL_transPCO_eqtlgen_results.txt",header = T)

#results_xaxis$gene = sapply(strsplit(results_xaxis$gene, "\\."), function(x) x[1])

merged_results = merge(results_xaxis,results_yaxis, by = 'gene')
merged_results = replace(merged_results, is.na(merged_results), 0)
merged_results[merged_results < 0] <- 0
print(sort(merged_results$r2.y[which(merged_results$r2.x + 0.05< merged_results$r2.y)]))
print(merged_results[order(merged_results$r2.x - merged_results$r2.y),])
merged_results$color <- ifelse(merged_results$r2.x < merged_results$r2.y, "dark green", "dark blue")
merged_results[, z_score := qnorm(merged_results$'r2 pval.x',lower.tail = FALSE)]
merged_results[, se := merged_results$r2.x/merged_results$z_score]
#print(merged_results$color)
print(mean(merged_results$r2.x,na.rm = T))
print(mean(merged_results$r2.y,na.rm = T))
print(merged_results)
print(sum(merged_results$r2.x > 0.4,na.rm = T))
print(sum(merged_results$r2.y > 0.4,na.rm = T))
print(sum(merged_results$r2.x > merged_results$r2.y))
print(sum(merged_results$r2.y > merged_results$r2.x))
pdf("plot of Whole_Blood cis versus cis_transPCO.pdf")
#plot(merged_results$r2.x,merged_results$r2.y, col = merged_results$color)
plot(merged_results$r2.x, merged_results$r2.y, col = alpha(merged_results$color,0.2), pch = 19,)
abline(a = 0,b = 1)
dev.off()
