library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--pval_threshold", action="store", default=NA, type='numeric',
              help="the threshold of which snps you want to include in the model")
  )

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$tissue)) {
        print("no tissue specified")
        q()
} else {
        tissue = opt$tissue
}

all_gene_info = fread("../data/GTEx_V8.txt.gz", header = T)
for (fold in 0:5) {
	pval_threshold = opt$pval_threshold
	output_dir = paste0("GBAT/",tissue,"/fold_",fold,"/results_",pval_threshold,"/")
	dir.create(output_dir)
 	association_results = fread(paste0("GBAT/",tissue,"/fold_",fold,"/association_results.txt"),header = T)
        association_results = association_results[which(association_results$'p-value' < pval_threshold),]
        if (nrow(association_results) == 0 ) { next }

        unique_genes = unique(association_results$gene)
        for (gene_name in unique_genes) {
		print(gene_name)
                subset = association_results[gene == gene_name,]
		
		if (nrow(subset) == 0) {
			next
		}

		gene_info = all_gene_info[ all_gene_info$geneId %in% gene_name]
		gene_chr = gene_info$'#chrom'
		gene_start = gene_info$chromStart
		sig_snps = data.table(snps = character(), pvals = numeric())
                print('here')
		for (row in 1:nrow(subset)) {
			if (subset$SNP[row] == gene_name) { next }
                   	#print(subset$SNP[row])
		   	trans_gene_info = all_gene_info[all_gene_info$geneId %in% subset$SNP[row]]
			trans_gene_chr = trans_gene_info$'#chrom'
			trans_gene_start = trans_gene_info$chromStart
			if ( (gene_chr == trans_gene_chr) & abs(gene_start - trans_gene_start) < 500000) {
				next
			}
			load(paste0("FUSION/",tissue,"/cis/",subset$SNP[row],".wgt.RDat"))
			snps = snps$V2[wgt.matrix[,colnames(wgt.matrix) == "lasso"] != 0]
			pvals = rep(subset$'p-value'[row],length(snps))
			sig_snps = rbind(sig_snps,data.table(snps,pvals))
		}
                if (nrow(sig_snps) == 0) {next}
		sig_snps = as.matrix(sig_snps)
		print(sig_snps)
		fwrite(sig_snps,paste0(output_dir,gene_name,".txt"),row.names = F, col.names = F, quote = F, sep = '\t')
        }
}

