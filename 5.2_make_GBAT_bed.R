library('data.table')

tissue="Whole_Blood"

for (fold in 0:5) {
	pval_threshold = 0.0001
	output_dir = paste0("GBAT/",tissue,"/fold_",fold,"/results_",pval_threshold,"/")
	dir.create(output_dir)
 	association_results = fread(paste0("GBAT/",tissue,"/fold_",fold,"/association_results.txt"),header = T)
        association_results = association_results[which(association_results$'p-value' < pval_threshold),]
        if (nrow(association_results) == 0 ) { next }

        unique_genes = unique(association_results$gene)
        for (gene_name in unique_genes) {
		print(gene_name)
                subset = association_results[gene == gene_name,]
		
                sig_snps = c()
		if (nrow(subset) == 0) {
			next
		}
                for (row in 1:nrow(subset)) {
			if (subset$SNP[row] == gene_name) { next }
                   	#print(subset$SNP[row])     
			load(paste0("FUSION/",tissue,"/cis/",subset$SNP[row],".wgt.RDat"))
			sig_snps = c(sig_snps,snps$V2[wgt.matrix[,colnames(wgt.matrix) == "lasso"] != 0])
                }
                if (length(sig_snps) == 0) {next}
		sig_snps = as.matrix(sig_snps)
		print(sig_snps)
		fwrite(sig_snps,paste0(output_dir,gene_name,".txt"),row.names = F, col.names = F, quote = F, sep = '\t')
        }
}

