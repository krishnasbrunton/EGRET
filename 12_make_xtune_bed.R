library('data.table')

genes = list.files("../data/rData/gtex_cis/")
allsnps = fread("../data/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim",header = F)

for (gene in genes) {
	gene_name = strsplit(gene,"\\.")[[1]][1]
	for (fold in 1:5) {
		bed_file = matrix(ncol = 6,nrow = 0)
		#cis
		cis_bed = fread(paste0("../data/rData/gtex_cis/",gene),header = F)
		
		#GBAT
		if (file.exists(paste0("fold_",fold,"/Whole_Blood/GBAT_results/",gene_name,".txt"))) {
			GBAT_snps = fread(paste0("fold_",fold,"/Whole_Blood/GBAT_results/",gene_name,".txt"),header = F)
		} else {
			GBAT_snps = matrix(ncol = 6,nrow = 0)
		}

		if (file.exists(paste0("BASIL_extra_large_cis/Whole_Blood/",gene_name,"/fold_",fold,"/selected_snps.txt")) & (file.size(paste0("BASIL_extra_large_cis/Whole_Blood/",gene_name,"/fold_",fold,"/selected_snps.txt")) > 1) ) {
			BASIL_snps = fread(paste0("BASIL_extra_large_cis/Whole_Blood/",gene_name,"/fold_",fold,"/selected_snps.txt"),header = F,fill = T)
		} else {
			BASIL_snps = matrix(ncol = 6,nrow = 0)
		}

		if (file.exists(paste0("../data/rData/transPCO_nocis/",gene_name,".bed")) ) {
                        transPCO_bed = fread(paste0("../data/rData/transPCO_nocis/",gene_name,".bed"),header = F,fill = T)
                } else {
                        transPCO_bed = matrix(ncol = 6,nrow = 0)
                }

		if (file.exists(paste0("../data/rData/eqtlgen_nocis/",gene_name,".bed")) ) {
                        eqtlgen_bed = fread(paste0("../data/rData/eqtlgen_nocis/",gene_name,".bed"),header = F,fill = T)
                } else {
                        eqtlgen_bed = matrix(ncol = 6,nrow = 0)
                }

		GBAT_matches = allsnps[match(unlist(GBAT_snps),allsnps$V2),]	
		if (nrow(GBAT_matches) == 0) {
			GBAT_bed = c()
		} else {
			GBAT_bed = cbind(GBAT_matches$V1,GBAT_matches$V4,GBAT_matches$V4, GBAT_matches$V2,"0","+")
		}
		
		BASIL_matches = allsnps[match(unlist(BASIL_snps),allsnps$V2),]
		if (nrow(BASIL_matches) == 0) {
			BASIL_bed = c()
		} else {
			BASIL_bed = cbind(BASIL_matches$V1,BASIL_matches$V4,BASIL_matches$V4, BASIL_matches$V2,"0","+")
		}
		

		bed_file = rbind(bed_file,cis_bed)
		bed_file = rbind(bed_file,GBAT_bed)
		#bed_file = rbind(bed_file,BASIL_bed)
		#bed_file = rbind(bed_file,transPCO_bed)
		#bed_file = rbind(bed_file,eqtlgen_bed)

		print(bed_file)	
		dir.create(paste0("xtune_bed_files/Whole_Blood/cis_GBAT/fold_",fold,"/"),recursive = T)
		fwrite(bed_file,paste0("xtune_bed_files/Whole_Blood/cis_GBAT/fold_",fold,"/",gene),row.names = F , col.names = F,quote = F, sep = '\t')
	}
}
