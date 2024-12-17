library(data.table)
library(optparse)

source("ModifiedPCOMerged.R")
source("ModifiedSigmaOEstimate.R")
source("davies.R")
source("liumod.R")
source("liu.R")
dyn.load( "qfc.so")

option_list = list(
  make_option("--module", action="store", default=NA, type='character',
              help="path of file containing modules")
  )

opt = parse_args(OptionParser(option_list=option_list))
module = opt$module



for (fold in 1:5) {
	for (chr in 1:22) {
		genes_in_module = fread(paste0("transPCO_modules/Whole_Blood/module_",module,"_gene_list.txt"),header = F)
		gene_info = fread("../data/GTEx_V8.txt.gz",header = T)

		genes_in_module = genes_in_module[gene_info$'#chrom'[gene_info$geneId %in% unlist(genes_in_module)] != paste0('chr',chr),]

		gene_expression = fread("gene_expression_cov_regressed/Whole_Blood_expression_regressed.txt.gz",header = T)	
		
		gene_expression_subset = gene_expression[gene_expression$gene_id %in% unlist(genes_in_module),]
		gene_expression_subset = gene_expression_subset[,-1]
		Sigma = cor(t(gene_expression_subset))

		
		Z.mat = fread(paste0("MatrixeQTL/Whole_Blood/fold_",fold,"/association_results_matrix_module_",module,"_chr_",chr,".txt.gz"),header = T)
		snps = Z.mat$SNP
		Z.mat <-setcolorder(Z.mat, unlist(genes_in_module))
		rownames(Z.mat) <- Z.mat$SNP  # Replace SNP_ID with the appropriate column name
		genes_in_module_vector <- unlist(genes_in_module)

		# Subset Z.mat to include only columns that are in 'genes_in_module'
		Z.mat <- Z.mat[, ..genes_in_module_vector, with = FALSE]


		SigmaO = ModifiedSigmaOEstimate(Sigma)
		pvals = ModifiedPCOMerged(Z.mat,Sigma,SigmaO)

		pvals = data.table(rsIDs = snps,pval = pvals)

		results_dir = paste0("transPCO_results/Whole_Blood/fold_",fold,"/")
		dir.create(results_dir,recursive = T)
		fwrite(pvals,paste0(results_dir,"module_",module,"_chr_",chr,".txt.gz"),sep = '\t', col.names = T, row.names = F, quote = F)
		rm(list = c('pvals','genes_in_module','Sigma','Z.mat'))
	}
}
