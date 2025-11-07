library('MatrixEQTL')
library('data.table')
library('optparse')
library('plink2R')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--module", action="store", default=NA, type='character',
              help="path of file containing modules")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue
module = opt$module

	for (fold in 0:5) {
		if (file.exists(paste0("transPCO/",tissue,"/fold_",fold,"/modules/",module))) {
			genes_in_module = fread(paste0("transPCO/",tissue,"/fold_",fold,"/modules/",module),header = F)
	      	} else { next }
		for (chr in 1:22) {
			print(module)
			print(chr)
			module_number = strsplit(module,"_")[[1]][2]
			if (file.exists(paste0("transPCO/",tissue,"/fold_",fold,"/association_results/module_",module_number,"_chr_",chr,".txt.gz"))) {
				next
			}

			gene_expression = fread(paste0("fold_",fold,"_info/",tissue,"/train_expression.txt"),header = T)
			gene_expression_module = gene_expression[match(unlist(genes_in_module$V1), gene_expression$gene_id),]
			dir.create(paste0("transPCO/",tissue,"/fold_",fold,"/expression_files/"),recursive = T)
			fwrite(gene_expression_module,paste0("transPCO/",tissue,"/fold_",fold,"/expression_files/module_",module_number,".txt"),sep = '\t',quote = F, col.names = T, row.names = F)
        		
			expression_file_name = paste0("transPCO/",tissue,"/fold_",fold,"/expression_files/module_",module_number,".txt")
        		gene = SlicedData$new()
        		gene$fileDelimiter = "\t"      # the TAB character
        		gene$fileOmitCharacters = "NA" # denote missing values
        		gene$fileSkipRows = 1          # one row of column labels
        		gene$fileSkipColumns = 1       # one column of row labels
        		gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        		gene$LoadFile( expression_file_name )

			dir.create(paste0("transPCO/",tissue,"/fold_",fold,"/genotype_files/"),recursive = T)
			if (!file.exists( paste0("transPCO/",tissue,"/fold_",fold,"/genotype_files/train_genotypes_chr_",chr,".txt"))) {
				train_individuals_path = paste0("fold_",fold,"_info/",tissue,"/train_individuals.txt")
        			arg = paste0("../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned --make-bed --keep ",train_individuals_path," --chr ",chr," --out transPCO/",tissue,"/fold_",fold,"/genotype_files/train_genotypes_chr_",chr)
        			system(arg)
        	
				genotypes = read_plink(paste0("transPCO/",tissue,"/fold_",fold,"/genotype_files/train_genotypes_chr_",chr),impute='avg')
        			genotypes$bed = scale(genotypes$bed)

        			rownames = colnames(genotypes$bed)
        			genotypes_transposed = t(genotypes$bed)
        			rownames(genotypes_transposed) = rownames
        			genotypes_transposed_dt = as.data.table(genotypes_transposed, keep.rownames = "SNP")

				fwrite(genotypes_transposed_dt,paste0("transPCO/",tissue,"/fold_",fold,"/genotype_files/train_genotypes_chr_",chr,".txt"),row.names = F,col.names = T,quote = F, sep = '\t')
	
			}
			
			SNP_file_name = paste0("transPCO/",tissue,"/fold_",fold,"/genotype_files/train_genotypes_chr_",chr,".txt")
        		snps = SlicedData$new()
        		snps$fileDelimiter = "\t"      # the TAB character
        		snps$fileOmitCharacters = "NA" # denote missing values
        		snps$fileSkipRows = 1          # one row of column labels
        		snps$fileSkipColumns = 1       # one column of row labels
        		snps$fileSliceSize = 20000      # read file in pieces of 2,000 rows
        		snps$LoadFile( SNP_file_name )
	
			errorCovariance = numeric()
			
			dir.create(paste0("transPCO/",tissue,"/fold_",fold,"/association_results/"),recursive = T)
			me = Matrix_eQTL_main(
                        	gene = gene,
                        	snps = snps,
                        	output_file_name = paste0("transPCO/",tissue,"/fold_",fold,"/association_results/module_",module_number,"_chr_",chr,".txt"),
                        	pvOutputThreshold = 1, #now set to 1 for transPCO analysis
                        	useModel = modelLINEAR,
                        	errorCovariance = errorCovariance,
                        	verbose = TRUE,
                        	pvalue.hist = TRUE,
                        	min.pv.by.genesnp = FALSE,
                        	noFDRsaveMemory = FALSE)

        	results = fread(paste0("transPCO/",tissue,"/fold_",fold,"/association_results/module_",module_number,"_chr_",chr,".txt"),header = T)
        	matrix_form = dcast(results, SNP~gene, value.var = "t-stat",fun.aggregate = mean, drop = FALSE)
        	fwrite(matrix_form, paste0("transPCO/",tissue,"/fold_",fold,"/association_results/matrix_module_",module_number,"_chr_",chr,".txt.gz"),sep = "\t", row.names = FALSE, col.names = TRUE)
		file.remove(paste0("transPCO/",tissue,"/fold_",fold,"/association_results/module_",module_number,"_chr_",chr,".txt"))
		rm(list = c("matrix_form", "results","snps","gene","genotypes","genotypes_transposed","genotypes_transposed_dt"))
		}

	}

