library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--pval_threshold", action="store", default=NA, type='numeric',
              help="the threshold of which snps you want to include in the model")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

gene_expression = fread(paste0("expression_files/",tissue,"_expression_regressed.txt.gz"),header = T)
p_value_threshold = opt$pval_threshold

for (fold in 0:5) {
	results = fread(paste0("MatrixeQTL/",tissue,"/association_results_fold_",fold,"_0.0001_threshold.txt"),header = T)
	results = results[results$'p-value' < p_value_threshold,]
	genes = unique(results$gene)
	for (current_gene in genes) {
		subset = results[gene == current_gene,]
		
		dir.create(paste0("MatrixeQTL/",tissue,"/fold_",fold,"/results_",p_value_threshold,"_threshold/"),recursive = T)
		snp_matrix = as.matrix(subset[, c("SNP", "p-value")])
		fwrite(snp_matrix,paste0("MatrixeQTL/",tissue,"/fold_",fold,"/results_",p_value_threshold,"_threshold/",current_gene,".txt") ,quote = F, col.names = F, row.names = F,sep = '\t')
		
	}

}
