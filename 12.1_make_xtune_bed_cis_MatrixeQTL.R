library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

genes = list.files(paste0("bed_files/cis/",tissue,"/"))
allsnps = fread("../data/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim",header = F)
threshold = 0.000001

for (gene in genes) {
	for (fold in 0:5) {
		bed_file = matrix(ncol = 6,nrow = 0)
		gene_name = strsplit(gene,".bed")[[1]][1]
		#cis
		cis_bed = fread(paste0("bed_files/cis/",tissue,"/",gene_name,".bed"),header = F)
		
		if (file.exists(paste0("MatrixeQTL/",tissue,"/fold_",fold,"/results_",threshold,"_threshold/",gene_name,".txt"))) {
                        MatrixeQTL_snps = fread(paste0("MatrixeQTL/",tissue,"/fold_",fold,"/results_",threshold,"_threshold/",gene_name,".txt"),header = F)
                } else {
                        MatrixeQTL_snps = matrix(ncol = 6,nrow = 0)
                }

		MatrixeQTL_matches = allsnps[match(unlist(MatrixeQTL_snps),allsnps$V2),]
		if (nrow(MatrixeQTL_matches) == 0) {
			MatrixeQTL_bed = c()
		} else {
			MatrixeQTL_bed = cbind(MatrixeQTL_matches$V1,MatrixeQTL_matches$V4,MatrixeQTL_matches$V4 + 1,MatrixeQTL_matches$V2,"0","+")
		}
		

		bed_file = rbind(bed_file,cis_bed)
		#bed_file = rbind(bed_file,GBAT_bed)
		#bed_file = rbind(bed_file,BASIL_bed)
		#bed_file = rbind(bed_file,transPCO_bed)
		#bed_file = rbind(bed_file,eqtlgen_bed)
		bed_file = rbind(bed_file,MatrixeQTL_bed)

		dir.create(paste0("xtune_bed_files/",tissue,"/cis_MatrixeQTL_",threshold,"/fold_",fold,"/"),recursive = T)
		fwrite(bed_file,paste0("xtune_bed_files/",tissue,"/cis_MatrixeQTL_",threshold,"/fold_",fold,"/",gene_name),row.names = F , col.names = F,quote = F, sep = '\t')
	}
}
