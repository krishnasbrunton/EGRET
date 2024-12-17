library("data.table")
library("optparse")

option_list = list(
        make_option("--expression", action="store",default=NA, type='character',
                help="expression file"),
        make_option("--out", action="store",default=NA, type='character',
                help="Path to store output"),
        make_option("--individuals", action="store",default=NA, type='character',
                help="Path to file containing list of individuals to keep"),
	make_option("--covariates", action="store",default=NA, type='character',
                help="Path to file of expresion covariates"),
	make_option("--tissue", action="store",default=NA, type='character',
                help="tissue for analysis")
  )

opt = parse_args(OptionParser(option_list=option_list))

tissue = opt$tissue

expression_file = fread(opt$expression, header = TRUE)
individuals = fread(opt$individuals, header = FALSE)
individuals = unlist(individuals$V1)

expression_file = expression_file[, c("gene_id", individuals), with = FALSE]

all_gene_info = fread("../data/GTEx_V8.txt.gz",header = T)
expression_file$geneType = all_gene_info$geneType[match(expression_file$gene_id,all_gene_info$geneId)]
expression_file = expression_file[expression_file$geneType == 'protein_coding']
print(dim(expression_file))
q()


fwrite(expression_file, paste0("expression_files/",tissue,"_expression.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

covariates = fread(opt$covariates, header = TRUE)

covariates = covariates[match(individuals, covariates$ID)]

fwrite(covariates,paste0("covariate_files/",tissue,"_covariates.txt.gz"),col.names = T, row.names = F, quote = F, sep = '\t')

residual_gene_expression = matrix(nrow = nrow(expression_file), ncol = ncol(expression_file))
colnames(residual_gene_expression) = colnames(expression_file)

for (row in 1:nrow(expression_file)) {
        reg = summary(lm( t(as.matrix(expression_file[row,-c('gene_id')])) ~ as.matrix(covariates[,c(3:ncol(covariates)),with = FALSE]) ))
        residual_gene_expression[row,2:ncol(residual_gene_expression)] = scale(reg$residuals) # scaled residuals
}

residual_gene_expression_df <- as.data.frame(residual_gene_expression)
residual_gene_expression_df$gene_id <- expression_file$gene_id

# Reorder columns to place gene_id first

# Write to file with fwrite
fwrite(residual_gene_expression_df, paste0("expression_files/",tissue,"_expression_regressed.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')





