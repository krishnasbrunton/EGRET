library('data.table')

#genes = list.files("bed_files/Whole_Blood/cis/")
genes = list.files("transPCO/Whole_Blood/fold_0/bed_files")
allsnps = fread("../data/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim",header = F)


for (gene in genes) {
	for (fold in 0:5) {
		bed_file = matrix(ncol = 6,nrow = 0)
	
		gene_name = substr(gene, 1, nchar(gene) - 4)
		gene = paste0(gene_name,".bed")
		print(gene_name)
		
		#cis
		cis_bed = fread(paste0("bed_files/Whole_Blood/cis/",gene),header = F)
		
		#GBAT
		GBAT_bed_path = paste0("GBAT/Whole_Blood/fold_",fold,"/results_1e-04/",gene_name,".txt")
		if (file.exists(GBAT_bed_path)) {
			GBAT_snps = fread(GBAT_bed_path,header = F)
		} else {
			GBAT_snps = matrix(ncol = 6,nrow = 0)
		}
		MatrixeQTL_bed_path = paste0("MatrixeQTL/Whole_Blood/fold_",fold,"/results_1e-06_threshold/",gene_name,".txt")
		if (file.exists(MatrixeQTL_bed_path)) {
			MatrixeQTL_snps = fread(MatrixeQTL_bed_path,header = F)
		} else {
			MatrixeQTL_snps = matrix(ncol = 6,nrow = 0)
		}
		transPCO_bed_path = paste0("transPCO/Whole_Blood/fold_",fold,"/bed_files/",gene_name,".txt")
		if (file.exists(paste0(transPCO_bed_path)) ) {
                        transPCO_snps = fread(transPCO_bed_path,header = F)
                } else {
                        transPCO_snps = matrix(ncol = 6,nrow = 0)
                }
		GBAT_matches = allsnps[match(unlist(GBAT_snps),allsnps$V2),]	
		if (nrow(GBAT_matches) == 0) {
			GBAT_bed = c()
		} else {
			GBAT_bed = cbind(GBAT_matches$V1,GBAT_matches$V4,GBAT_matches$V4+1, GBAT_matches$V2,"0","+")
		}
		MatrixeQTL_matches = allsnps[match(unlist(MatrixeQTL_snps),allsnps$V2),]
		if (nrow(MatrixeQTL_matches) == 0) {
			MatrixeQTL_bed = c()
		} else {
			MatrixeQTL_bed = cbind(MatrixeQTL_matches$V1,MatrixeQTL_matches$V4,MatrixeQTL_matches$V4+1, MatrixeQTL_matches$V2,"0","+")
		}
		transPCO_matches = allsnps[match(unlist(transPCO_snps),allsnps$V2),]
		if(nrow(transPCO_matches) == 0) {
			transPCO_bed = c()
		} else {
			transPCO_bed = cbind(transPCO_matches$V1,transPCO_matches$V4,transPCO_matches$V4+1,transPCO_matches$V2, "0","+")
		}
		
		bed_file = rbind(bed_file,cis_bed)
		#bed_file = rbind(bed_file,MatrixeQTL_bed)
		#bed_file = rbind(bed_file,GBAT_bed)
		bed_file = rbind(bed_file,transPCO_bed)


		bed_file_dir = paste0("bed_files/Whole_Blood/cis_transPCO/fold_",fold,"/")
		dir.create(bed_file_dir,recursive = T)
		fwrite(bed_file,paste0(bed_file_dir,gene),row.names = F , col.names = F,quote = F, sep = '\t')

		plink_output_dir = paste0("plink_results/Whole_Blood/cis_transPCO/fold_",fold,"/")
		dir.create(plink_output_dir,recursive = T)
		
		arg = paste0("../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned --extract bed1 ",bed_file_dir,gene, " --make-bed --out ",plink_output_dir,gene_name)
                system(arg)

		next
		if(file.exists(paste0("bed_files/cross_mapped/4_mismatches/",gene))) {
			print("exluding cross mappable genes")
			arg = paste0("../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned --extract bed1 ",bed_file_dir,gene, " --make-bed --exclude bed1 bed_files/cross_mapped/4_mismatches/",gene, " --out ",plink_output_dir,gene_name)
		} else {
			print("no cross mappable genes to exlude")
			arg = paste0("../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned --extract bed1 ",bed_file_dir,gene, " --make-bed --out ",plink_output_dir,gene_name)
		}
		system(arg)
	}
}
