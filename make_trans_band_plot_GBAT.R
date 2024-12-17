library(data.table)

snp_info = fread("../data/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim",header = F)
gene_info = fread("../data/GTEx_V8.txt.gz",header = T)

gene_info[,gene_id_split := sapply(strsplit(geneId, "\\."), `[`, 1)]

results = data.frame(gene = as.character(), snps = as.character())
GBAT_result_files = list.files("GBAT/Whole_Blood/fold_0/results_1e-06/")
for (file in GBAT_result_files) {
	gene_results = fread(paste0("GBAT/Whole_Blood/fold_0/results_1e-06/",file),header = F)
	gene_name = strsplit(file, ".txt")[[1]][1]
	print(gene_name)
	results = rbind(results,data.frame(gene = gene_name,snps = unlist(gene_results$V1)))

}
results$snp_chr = snp_info$V1[match(results$SNP,snp_info$V2)]
results$snp_row = match(results$SNP,snp_info$V2)

results$gene_chr = gene_info$'#chrom'[match(results$gene,gene_info$geneId)]

results$gene_chr = sapply(results$gene_chr, function(x) strsplit(x, "chr")[[1]][2])
results$gene_row = match(results$gene,gene_info$geneId)



results$snp_chr = snp_info$V1[match(results$snps,snp_info$V2)]
results$snp_row = match(results$snps,snp_info$V2)

results$gene_row = match(results$gene,gene_info$gene_id_split)




alpha_value <- 0.1  # You can adjust this value for more or less transparency

# Define a semi-transparent blue color
point_color <- rgb(0, 0, 1, alpha = alpha_value)

# Plot
pdf("GBAT trans-band_scatter.pdf")
plot(results$snp_row, results$gene_row, 
     main = "Scatter Plot of SNP Position vs Gene row", 
     xlab = "SNP Position", 
     ylab = "Gene Position", 
     pch = 16,  
     col = point_color)      # Optional: Add color
dev.off()

