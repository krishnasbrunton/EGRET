library('data.table')
library('optparse')
library('plink2R')

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--cis_model_dir", action="store",default=NA, type='character',
              help="directory containing previously trained FUSION models"),
  make_option("--gene_info", action="store",default=NA, type='character',
        	  help="file containg gene info such as gene id, gene name, start, chromosome, and gene type"),
  make_option("--output_dir", action="store",default=NA, type='character',
              help="directory where results are stored")
  )
  

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$tissue)) {
	print("no tissue specified")
	q()
} else {
	tissue = opt$tissue
}

gene_info = fread(opt$gene_info, header = T, sep = '\t')

trans_genes = list.files(paste0(opt$cis_model_dir),pattern = ".RDat")
cis_predicted_expression = data.table()

for (trans_gene in trans_genes) {

	gene_name = paste0(strsplit(trans_gene,"\\.")[[1]][1:2],collapse = ".")
	print(gene_name)

    load(paste0(opt$cis_model_dir,trans_gene))
	model_pval = cv.performance[2,colnames(cv.performance) == "lasso"]

	if (is.na(model_pval) | model_pval > 0.01) {
		next
	}

	print(paste("the r2 is",cv.performance[1,1]))


	# based on gene position, create bed file
	# extract genotype snps for specific gene using plink
	# read in genotypes and predict expression for that particular gene


	geno.file = paste0("plink_results/",tissue,"/cis/",gene_name)
	if (!file.exists(paste0(geno.file,".bim"))) {
		next
	}

	genos = read_plink(geno.file,impute="avg")
    genos_snps = scale(genos$bed)

	cv.all = fread(paste0("working/",tissue,"/cis/",gene_name,".fam"),header = F)
    cv.all = cv.all[,c(2,6)]
	

	if (any(is.na(wgt.matrix))) {
		next
	}

	matching_indices = match(rownames(wgt.matrix),colnames(genos_snps))


	matching_indices = matching_indices[!is.na(matching_indices)]
	print(dim(genos_snps))
	genos_snps = genos_snps[,matching_indices]
	print(dim(genos_snps))	

	matching_wgt_indices = match(colnames(genos_snps),rownames(wgt.matrix))
	wgt.matrix = as.matrix(wgt.matrix[matching_wgt_indices,])
	print(dim(wgt.matrix))

	if (is.null(genos_snps) || is.null(wgt.matrix) || nrow(genos_snps) == 0 || nrow(wgt.matrix) == 0) {
		next
	}	
	
	predictions = genos_snps %*% wgt.matrix[,colnames(wgt.matrix) == "lasso"]
	
	if (any(is.na(predictions))) {
                next
        }

	matching_ppl = match(paste0("0:", unlist(cv.all[,1])),rownames(predictions))
	matching_ppl = matching_ppl[!is.na(matching_ppl)]
	predictions = as.matrix(predictions[matching_ppl, ])

	matching_ppl = match(rownames(predictions),paste0("0:", unlist(cv.all[,1])))
	matching_ppl = matching_ppl[!is.na(matching_ppl)]
	cv.all = cv.all[matching_ppl,]
	
	fit = summary(lm(as.matrix(cv.all[,2]) ~ predictions))

	print(paste0("The r2 of the fitted model is : ", fit$adj.r.sq))
	predictions = rbind(gene_name, predictions)
	rownames(predictions)[1] = "gene"
	cis_predicted_expression = rbind(cis_predicted_expression,t(predictions))

}

dir.create(paste0("GBAT/",tissue,"/"), recursive = T)
fwrite(cis_predicted_expression,paste0("GBAT/",tissue,"/cis_predicted_expression.txt"), col.names = T, row.names = F, sep = '\t', quote = F)

