library(data.table)

genes = fread("../data/GTEx_V8.txt.gz",header = T)
print(unique(genes$geneType))
print(nrow(genes[genes$geneType == c("pseudogene")]))
print(nrow(genes[genes$geneType == c("protein_coding")]))


expression = fread("expression_files/Whole_Blood_expression.txt.gz")

expression$geneType = genes$geneType[match(expression$gene_id,genes$geneId)]

print(unique(expression$geneType))
gene_types = unique(expression$geneType)
pseudo_gene_types = gene_types[grep('pseudogene',unlist( gene_types))]
print(pseudo_gene_types)
print(nrow(expression[!expression$geneType %in% pseudo_gene_types]))
#expression = expression[!expression$geneType %in% pseudo_gene_types]
print(nrow(expression[expression$geneType == c("pseudogene")]))
print(nrow(expression[expression$geneType == c("unprocessed_pseudogene")]))

#results_yaxis = fread("../trans_adapt/results_sumstats/Whole_Blood/cis_MatrixeQTL_1e-06_redo.txt")
results_yaxis = fread("results_sumstats/Whole_Blood/cis_MatrixeQTL_1e-06_0_mismatches.txt")

results_xaxis = fread("results_sumstats/Whole_Blood/cis_results.txt",header = T)
#results_xaxis$gene = sapply(strsplit(results_xaxis$gene, "\\."), `[`, 1)

merged_results = merge(results_xaxis,results_yaxis, by = 'gene')
merged_results = replace(merged_results, is.na(merged_results), 0)
merged_results[merged_results < 0] <- 0
top_genes = merged_results[which(merged_results$r2.x + 0.05< merged_results$r2.y)]
top_genes = top_genes[order(top_genes$r2.y - top_genes$r2.x),]
top_genes$geneType = expression$geneType[match(top_genes$gene, expression$gene_id)]

#top_genes$geneType = expression$geneType[match(top_genes$gene, sapply(strsplit(expression$gene_id, "\\."), `[`, 1))]
print(top_genes)
print(sum(top_genes$geneType %in% pseudo_gene_types))
q()
print(expression$geneType[match(top_genes,expression$gene_id)])
