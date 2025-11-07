library('MatrixEQTL')
library('optparse')

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis")
  )

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$tissue)) {
        print("no tissue specified")
        q()
} else {
        tissue = opt$tissue
}



for (fold in 0:5) {

	expression_file_name = paste0("fold_",fold,"_info/",tissue,"/train_expression.txt")
        gene = SlicedData$new()
        gene$fileDelimiter = "\t"      # the TAB character
        gene$fileOmitCharacters = "NA" # denote missing values
        gene$fileSkipRows = 1          # one row of column labels
        gene$fileSkipColumns = 1       # one column of row labels
        gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        gene$LoadFile( expression_file_name )

	SNP_file_name = paste0("fold_",fold,"_info/",tissue,"/train_genotypes.txt")
	snps = SlicedData$new()
	snps$fileDelimiter = "\t"      # the TAB character
	snps$fileOmitCharacters = "NA" # denote missing values
	snps$fileSkipRows = 1          # one row of column labels
	snps$fileSkipColumns = 1       # one column of row labels
	snps$fileSliceSize = 20000      # read file in pieces of 2,000 rows
	snps$LoadFile( SNP_file_name )

	errorCovariance = numeric()

	dir.create(paste0("MatrixeQTL/",tissue,"/"),recursive = T)
	me = Matrix_eQTL_main(
                        gene = gene,
                        snps = snps,
                        output_file_name = paste0("MatrixeQTL/",tissue,"/association_results_fold_",fold,"_0.0001_threshold.txt"),
                        pvOutputThreshold = 0.00001, 
                        useModel = modelLINEAR,
                        errorCovariance = errorCovariance,
                        verbose = TRUE,
                        pvalue.hist = TRUE,
                        min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE)
	
}


