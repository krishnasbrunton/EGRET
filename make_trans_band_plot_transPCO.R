library(data.table)

snp_info = fread("../data/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim",header = F)
gene_info = fread("../data/GTEx_V8.txt.gz",header = T)

num_modules = length(list.files(paste0("transPCO/Whole_Blood/modules/fold_0/")))
sig_results = data.table(module = numeric(), rsIDs = character(), pval = numeric())
pval_threshold = 0.00001
fold = 0
for (module in 1:num_modules) {
                print(module)
       for (chr in 1:22) {
                        if (!file.exists(paste0("transPCO/Whole_Blood/fold_",fold,"/PCO_association_results/module_",module,"_chr_",chr,".txt.gz"))) {
                                next
                        }
                        results = fread(paste0("transPCO/Whole_Blood/fold_",fold,"/PCO_association_results/module_",module,"_chr_",chr,".txt.gz"), header = T)
                        sig_rows = results[results$pval < pval_threshold,]
                        if (nrow(sig_rows) > 0) {
                                sig_results = rbind(sig_results, data.table(module = module, sig_rows))
                        }
                }

}

sig_results = sig_results[sig_results$pval != 0,]

results = sig_results
results$snp_chr = snp_info$V1[match(results$rsIDs,snp_info$V2)]
results$snp_row = match(results$rsIDs,snp_info$V2)


min_pval <- min(results$'pval', na.rm = TRUE)
max_pval <- max(results$'pval', na.rm = TRUE)

# Inverse scaling (smaller p-values get larger sizes)
point_size <- 1 + (log10(max_pval) - log10(results$'pval')) / (log10(max_pval) - log10(min_pval)) * 5

alpha_value <- 0.1  # You can adjust this value for more or less transparency

# Define a semi-transparent blue color
point_color <- rgb(0, 0, 1, alpha = alpha_value)

# Plot
pdf("transPCO trans-band_scatter.pdf")
plot(results$snp_row, results$module, 
     main = "Scatter Plot of SNP Position vs Gene Module", 
     xlab = "SNP Position", 
     ylab = "Gene Module", 
     pch = 16, 
     cex = point_size,  # Adjust point size based on p-value
     col = point_color)      # Optional: Add color
dev.off()

