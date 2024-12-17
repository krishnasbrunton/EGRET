library(data.table)

snp_info = fread("../data/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim",header = F)
gene_info = fread("../data/GTEx_V8.txt.gz",header = T)

results = fread("MatrixeQTL/Whole_Blood/association_results_fold_0_0.0001_threshold.txt", header = T)
results = results[results$'p-value' < 0.000000001,]
results$snp_chr = snp_info$V1[match(results$SNP,snp_info$V2)]
results$snp_row = match(results$SNP,snp_info$V2)

results$gene_chr = gene_info$'#chrom'[match(results$gene,gene_info$geneId)]

results$gene_chr = sapply(results$gene_chr, function(x) strsplit(x, "chr")[[1]][2])
results$gene_row = match(results$gene,gene_info$geneId)


results = results[results$gene_chr != results$snp_chr,]
min_pval <- min(results$'p-value', na.rm = TRUE)
max_pval <- max(results$'p-value', na.rm = TRUE)

# Inverse scaling (smaller p-values get larger sizes)
point_size <- 1 + (log10(max_pval) - log10(results$'p-value')) / (log10(max_pval) - log10(min_pval)) * 5

alpha_value <- 0.1  # You can adjust this value for more or less transparency

# Define a semi-transparent blue color
point_color <- rgb(0, 0, 1, alpha = alpha_value)

# Plot
pdf("MatrixeQTL trans-band_scatter.pdf")
plot(results$snp_row, results$gene_row, 
     main = "Scatter Plot of SNP Position vs Gene", 
     xlab = "SNP Position", 
     ylab = "Gene", 
     pch = 16, 
     cex = point_size,  # Adjust point size based on p-value
     col = point_color)      # Optional: Add color
dev.off()

