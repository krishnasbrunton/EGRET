library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--FDR", action="store", default=0.1, type='numeric',
              help="FDR rate to include snps into the model")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

all_gene_info = fread("../data/GTEx_V8.txt.gz", header = T)
all_snp_info = fread("genotype_files/GTEX_v8_genotypes_pruned.bim", header = F)

gene_expression = fread(paste0("expression_files/",tissue,"_expression_regressed.txt.gz"),header = T)

for (fold in 0:5) {
	results = fread(paste0("MatrixeQTL/",tissue,"/association_results_fold_",fold,"_0.001_threshold.txt"),header = T)
	genes = unique(results$gene)
	output_dir = paste0("MatrixeQTL/",tissue,"/fold_",fold,"/results_FDR_",opt$FDR,"/")
	largest_sig_pval = c()
	for (current_gene in genes) {
		subset = results[gene == current_gene,]
		subset$new_FDR = p.adjust(subset$'p-value',n = nrow(all_snp_info), method = 'fdr')
		subset = subset[subset$new_FDR < opt$FDR,]
		largest_sig_pval = c(largest_sig_pval, subset$'p-value'[nrow(subset)])
		gene_info = all_gene_info[match(current_gene,all_gene_info$geneId),]
		chr = strsplit(gene_info$'#chrom','chr')[[1]][2]
		u_bound = max(0,gene_info$chromStart - 500000)
		l_bound = gene_info$chromStart + 500000
		snp_info = all_snp_info[match(subset$SNP,all_snp_info$V2),]
		cis_snps = which(snp_info$V1 == chr & snp_info$V4 > u_bound & snp_info$V4 < l_bound)
		subset = subset[!cis_snps,]
		
		dir.create(output_dir,recursive = T)
		snp_matrix = as.matrix(subset[, c("SNP", "p-value")])
		if (nrow(snp_matrix) == 0) { next }
		fwrite(snp_matrix,paste0(output_dir,current_gene,".txt") ,quote = F, col.names = F, row.names = F,sep = '\t')
		
	}
	print(mean(largest_sig_pval,na.rm = T))

}
