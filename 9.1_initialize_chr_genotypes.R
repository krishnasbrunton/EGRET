library(data.table)
library(plink2R)
tissue = "Whole_Blood"


for (fold in 0:5) {
	for (chr in 1:22) {

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


}
