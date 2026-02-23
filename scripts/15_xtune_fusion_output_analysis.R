library('data.table')
library('optparse')

option_list = list(
  make_option("--results_dir", action="store", default=NA, type='character',
              help="Path to fusion results directory"),
  make_option("--plink_dir", action="store",default=NA, type='character',
	      help='Path to plink files'),
  make_option("--weight_dir", action="store",default=NA, type='character',
              help='Path to weight files'),
  make_option("--output_file", action="store",default=NA, type='character',
              help='Path to output'),
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue")

  )
opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue


if (!is.na(opt$results_dir)) {
	file_dir = opt$results_dir
} else {
	#manually set fusion file directory
	file_dir = paste0("xtune_fusion_results/",opt$tissue,"/cis_transPCO_old/")
}

if (!is.na(opt$plink_dir)) {
	plink_dir = opt$plink_dir
} else {
	#manually set plink file directory
	plink_dir = paste0("plink_results/",opt$tissue,"/cis_transPCO_old/")
}

if (!is.na(opt$weight_dir)) {
        weight_dir = opt$weight_dir
} else {
	#manually set weight file directory
	weight_dir = paste0("weights/",opt$tissue,"/cis_transPCO_old/")
}
if (!is.na(opt$output_file)) {
	output_file = opt$output_file
} else {
	#manually set output file
	output_file = paste0("results_sumstats/",opt$tissue,"/cis_transPCO_old.txt")
}


fusion_result_files = list.files(file_dir)
pattern_to_exclude = "alphas"
fusion_result_files = fusion_result_files[!grepl(pattern_to_exclude, fusion_result_files)]
print(file_dir)
#fusion_result_files = fusion_result_files[1:10]
model_results = c()
for (gene in fusion_result_files) {
	gene = substr(gene,1,nchar(gene)-4)
  	gene_info = readLines(paste0(file_dir,gene,".txt"))
	print(gene)

  	if (file.size(file = paste0(file_dir,gene,".txt")) == 0 | length(gene_info) == 0) {
    		next
  	}
 	
  	heritability_match = grep("Heritability",gene_info)
  	r2_match = grep("rsq",gene_info)
  	pvalue_match = grep("pval",gene_info)

  	if (length(heritability_match) != 0) {
    		heritability_info = strsplit(gene_info[heritability_match]," ")[[1]]
  	} else {
    		heritability_info = c("NA","NA","NA","NA","NA","NA","NA")
  	}

  	if (length(r2_match) != 0) {
    		r2_info = strsplit(gene_info[r2_match],"\t")[[1]]
  	} else {
    		r2_info = c("NA","NA")
  	}

  	if (length(pvalue_match) != 0) {
    		pvalue_info = strsplit(gene_info[pvalue_match],"\t")[[1]]
    		#print(pvalue_info)
  	} else {
    		pvalue_info = c("NA","NA")
  	}
  	if (file.exists(paste0(plink_dir,"/",gene,".bim"))) {
		num_snps = nrow(fread(file = paste0(plink_dir,"/",gene,".bim")))
  	} else {
		num_snps = "NA"
	}
  
  	if (file.exists(paste0(weight_dir,"/",tissue,".",gene,".wgt.RDat"))) {
    		num_nonzero_snps = get(load(paste0(weight_dir,"/",tissue,".",gene,".wgt.RDat"))[1])
    		num_nonzero_snps = length(which(num_nonzero_snps != 0))
  	} else {
    		num_nonzero_snps = 0
  	}
	model_results = rbind(model_results,c(gene,heritability_info[3],heritability_info[4],heritability_info[7],num_snps,r2_info[2],pvalue_info[2],num_nonzero_snps))
}
colnames(model_results) <- c("gene","h2","se","h2 pval","number of snps","r2","r2 pval","number non-zero snps")

dir.create(paste0("results_sumstats/",opt$tissue,"/"),recursive = T)
#print(model_results)
write.table(model_results,file = output_file,col.names = TRUE,row.names = FALSE, sep = "\t",quote= FALSE)






