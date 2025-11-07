library('data.table')
library('optparse')
library('plink2R')

option_list = list(
	make_option("--expression", action="store",default=NA, type='character',
		help="Path to gene expression matrix"),
	make_option("--individuals", action="store",default=NA, type='character',
                help="Path to file containing list of individuals"),
        make_option("--tissue", action = "store", default = NA, type = 'character',
		help="tissue that will set up for"),
	make_option("--covariates", action="store",default=NA, type='character',
                help="Path to file containing list of individuals"),
	make_option("--folds", action="store",default=5, type='numeric',
                help="Number of folds"),
	make_option("--plink_path", action="store",default=NA, type='character',
                help="Path to plink"),
	 make_option("--bfile", action="store",default=NA, type='character',
                help="Path to genotype files")
		   
  )

opt = parse_args(OptionParser(option_list=option_list))

tissue = opt$tissue
cross_val_folds = opt$folds

gene_expression = fread(opt$expression,header = T)
covar_file = fread(opt$covariates,header = T)

### SETUP CROSS-VALIDATION INFO ###
all_individuals = fread(opt$individuals,header = F)

N = nrow(all_individuals)

folds = cut(seq(1,N),breaks = cross_val_folds,labels=FALSE)

cross_val_results = matrix(NA,nrow=N,ncol=1)

for (i in 0:cross_val_folds) {
	train_individuals = all_individuals[which(folds != i,arr.ind=TRUE),]
	test_individuals = all_individuals[which(folds == i,arr.ind=TRUE),]

	dir.create(paste0("fold_",i,"_info/",tissue,"/"),recursive = T)

	fwrite(train_individuals, paste0("fold_",i,"_info/",tissue,"/train_individuals.txt"),row.names = F, col.names = F, quote = F, sep = '\t')
	fwrite(test_individuals, paste0("fold_",i,"_info/",tissue,"/test_individuals.txt"),row.names = F, col.names = F, quote = F, sep = '\t')
	train_gene_expression = gene_expression[, c(1,match(unlist(train_individuals), colnames(gene_expression))), with = FALSE]

	fwrite(train_gene_expression,paste0("fold_",i,"_info/",tissue,"/train_expression.txt"),col.names = T,row.names = F,quote = F, sep = '\t')
		
	train_covar = covar_file[match(unlist(train_individuals), covar_file$ID),]
	#write.table(t(train_covar),paste0("fold_",i,"_info/",tissue,"/train_covar.txt"),col.names = F,row.names = T,quote = F, sep = '\t')

	arg = paste0(opt$plink_path ," --bfile ", opt$bfile," --make-bed --keep fold_",i,"_info/",tissue,"/train_individuals.txt --out fold_",i,"_info/",tissue,"/train_genotypes")
        system(arg)

        genotypes = read_plink(paste0("fold_",i,"_info/",tissue,"/train_genotypes"),impute='avg')
        genotypes$bed = scale(genotypes$bed)

        rownames = colnames(genotypes$bed)
        genotypes_transposed = t(genotypes$bed)
        rownames(genotypes_transposed) = rownames
        genotypes_transposed_dt = as.data.table(genotypes_transposed, keep.rownames = "SNP")

        fwrite(genotypes_transposed_dt,paste0("fold_",i,"_info/",tissue,"/train_genotypes.txt"),row.names = F,col.names = T,quote = F, sep = '\t')
}




