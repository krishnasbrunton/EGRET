library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--pval_threshold", action="store", default=0.000001, type='numeric',
              help="the threshold of which snps you want to include in the model"),
  make_option("--PCO_bed_dir", action="store", default=NA, type='character',
              help="the dir which contains the PCO bed files"),
  make_option("--exclude_crossmap", action="store", default=TRUE, type='logical',
              help="if True, will exclude regions which are cross mappable to that gene")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue
pval_threshold = opt$pval_threshold

#genes = list.files(paste("transPCO/Whole_Blood/fold_1/bed_files_60_PCs_power_3"))
genes = list.files(paste0("bed_files/",tissue,"/cis/"))
allsnps = fread("genotype_files/GTEX_v8_genotypes_pruned.bim",header = F)

output_dir = "cis_MatrixeQTL_GBAT_transPCO_4_mismatches"

crossmap_dir = "4_mismatches"
for (gene in genes) {
	for (fold in 0:5) {
		bed_file = matrix(ncol = 6,nrow = 0)
	
		gene_name = substr(gene, 1, nchar(gene) - 4)
		gene = paste0(gene_name,".bed")
		print(gene_name)
		
		#cis
		cis_bed_path = paste0("bed_files/",tissue,"/cis/",gene)
		if (file.exists(cis_bed_path)) {
                        cis_bed = fread(cis_bed_path,header = F)
                } else {
                        next
                }

		
		#GBAT
		GBAT_bed_path = paste0("GBAT/",tissue,"/fold_",fold,"/results_1e-04/",gene_name,".txt")
		if (file.exists(GBAT_bed_path)) {
			GBAT_snps = fread(GBAT_bed_path,header = F)
		} else {
			GBAT_snps = matrix(ncol = 6,nrow = 0)
		}
		MatrixeQTL_bed_path = paste0("MatrixeQTL/",tissue,"/fold_",fold,"/results_",pval_threshold,"_threshold/",gene_name,".txt")
		if (file.exists(MatrixeQTL_bed_path)) {
			MatrixeQTL_snps = fread(MatrixeQTL_bed_path,header = F)
		} else {
			MatrixeQTL_snps = matrix(ncol = 6,nrow = 0)
		}
		transPCO_bed_path = paste0("transPCO/",tissue,"/fold_",fold,"/",opt$PCO_bed_dir,"/",gene_name,".txt")
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
			transPCO_matches = transPCO_matches[!is.na(transPCO_matches$V2),]
			transPCO_bed = cbind(transPCO_matches$V1,transPCO_matches$V4,transPCO_matches$V4+1,transPCO_matches$V2, "0","+")
		}
		
		bed_file = rbind(bed_file,cis_bed)
		bed_file = rbind(bed_file,MatrixeQTL_bed)
		bed_file = rbind(bed_file,GBAT_bed)
		bed_file = rbind(bed_file,transPCO_bed)

		bed_file_dir = paste0("bed_files/",tissue,"/",output_dir,"/fold_",fold,"/")
		dir.create(bed_file_dir,recursive = T)
		fwrite(bed_file,paste0(bed_file_dir,gene),row.names = F , col.names = F,quote = F, sep = '\t')

		plink_output_dir = paste0("plink_results/",tissue,"/",output_dir,"/fold_",fold,"/")
		dir.create(plink_output_dir,recursive = T)		

		if(file.exists(paste0("bed_files/cross_mapped/",crossmap_dir,"/",gene)) & opt$exclude_crossmap) {
			print("exluding cross mappable genes")
			arg = paste0("../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned --extract bed1 ",bed_file_dir,gene, " --make-bed --exclude bed1 bed_files/cross_mapped/",crossmap_dir,"/",gene, " --out ",plink_output_dir,gene_name)
		} else {
			print("no cross mappable genes to exlude")
			arg = paste0("../../kakamatsu/plink2 --bfile genotype_files/GTEX_v8_genotypes_pruned --extract bed1 ",bed_file_dir,gene, " --make-bed --out ",plink_output_dir,gene_name)
		}

		system(arg)

		if(file.exists(paste0(plink_output_dir,gene_name,".bim"))) {
			plink_results = fread(paste0(plink_output_dir,gene_name,".bim"), header =  F)
		} else { next }
		print(paste0("plink_results/",tissue,"/cis/fold_",fold,"/",gene_name,".bim"))
		if(file.exists(paste0("plink_results/",tissue,"/cis/",gene_name,".bim"))) {
			cis_plink_results = fread(paste0("plink_results/",tissue,"/cis/",gene_name,".bim"), header = F)
		}

		cis_matches = match(cis_plink_results$V2,plink_results$V2)
		plink_results$cis = 0
		plink_results$cis[cis_matches] =  1

		MatrixeQTL_matches = match(unlist(MatrixeQTL_snps),plink_results$V2)
		plink_results$MatrixeQTL = 0
		plink_results$MatrixeQTL[MatrixeQTL_matches] = 1
		
		GBAT_matches = match(unlist(GBAT_snps),plink_results$V2)
		plink_results$GBAT = 0
		plink_results$GBAT[GBAT_matches] = 1

		transPCO_matches = match(unlist(transPCO_snps),plink_results$V2)
		plink_results$transPCO = 0
		plink_results$transPCO[transPCO_matches] = 1

		z_matrix <- data.table(plink_results$cis,plink_results$GBAT,plink_results$MatrixeQTL,plink_results$transPCO)
		z_matrix_dir = paste0("z_matrices/",tissue,"/",output_dir,"/fold_",fold,"/")
		dir.create(z_matrix_dir, recursive = T)
		fwrite(z_matrix,paste0(z_matrix_dir,gene_name,"_z_matrix.txt"),sep = '\t',col.names = F,row.names = F,quote=F)
	}
}
