library(data.table)

all_gene_info = fread("../data/GTEx_V8.txt.gz", header = T)

cis_sumstats = fread("results_sumstats/Whole_Blood/cis_results.txt",header = T)
MatrixeQTL_sumstats = fread("../trans_adapt/results_sumstats/Whole_Blood/cis_MatrixeQTL_1e-06_redo.txt", header = T)
transPCO_sumstats = fread("../transPCO/results_sumstats/Whole_Blood/cis_transPCO_results.txt", header = T)
GBAT_sumstats = fread("results_sumstats/Whole_Blood/cis_GBAT_1e-04.txt", header = T)

cis_wgt_dir = "FUSION/Whole_Blood/cis/"   #ENSG00000000457.13.wgt.RDat
MatrixeQTL_wgt_dir = "../trans_adapt/weights/Whole_Blood/cis_MatrixeQTL_1e-06/"   #ENSG00000000457.wgt.RDat
transPCO_wgt_dir = "../transPCO/weights/Whole_Blood/cis_transPCO/"    #ENSG00000000457.wgt.RDat
GBAT_wgt_dir = "xtune_fusion_models/Whole_Blood/cis_GBAT_1e-04/"     #ENSG00000000457.13.wgt.RDat

new_wgt_dir = "xtune_fusion_models/Whole_Blood/cis_GBAT/"
dir.create(new_wgt_dir,recursive = T)

cis_sumstats$gene_full_cis = cis_sumstats$gene
cis_sumstats$gene = sapply(strsplit(cis_sumstats$gene, "\\."), function(x) x[1])
GBAT_sumstats$gene_full_GBAT = GBAT_sumstats$gene
GBAT_sumstats$gene = sapply(strsplit(GBAT_sumstats$gene, "\\."), function(x) x[1])

merge_one = merge(cis_sumstats, MatrixeQTL_sumstats,by = 'gene', all.x = T, all.y = T,suffixes = c('.cis','.matrix'))
merge_two = merge(merge_one, transPCO_sumstats,by = 'gene', all.x = T, all.y = T)
merge_three = merge(merge_two,GBAT_sumstats,by = 'gene',all.x = T, all.y = T, suffixes = c('.PCO','.GBAT'))

merge_three[is.na(merge_three)] = 0

r2_columns = c("r2.cis", "r2.GBAT")
merge_three$best_model = apply(merge_three[, ..r2_columns], 1, function(row) {
  # Use which.max to find the index of the maximum r2 value
  colnames(merge_three[, ..r2_columns])[which.max(row)]
})

merge_three$best_r2 = apply(merge_three[, ..r2_columns], 1, max)

print(merge_three)
pos_matrix = data.frame(WGT = as.character(),ID = as.character(),CHR = as.numeric(), P0 = as.numeric(),P1 = as.numeric())

for (row in 1:nrow(merge_three)) {
    print(row)
    print(merge_three$best_r2[row])
    if (merge_three$best_r2[row] < 0.01) {
        next
    }
    if (merge_three$best_model[row] == "r2.cis") {
        print('cis is the best')
        arg = paste0("cp ",cis_wgt_dir,merge_three$gene_full_cis[row], ".wgt.RDat ",new_wgt_dir,merge_three$gene[row],".wgt.RDat ")
    } else if (merge_three$best_model[row] == 'r2.matrix') {
        print('MatrixeQTL is the best')
        arg = paste0("cp ",MatrixeQTL_wgt_dir,merge_three$gene[row],".wgt.RDat ",new_wgt_dir,merge_three$gene[row],".wgt.RDat ")
    } else if (merge_three$best_model[row] == 'r2.PCO') {
        print('transPCO is the best')
        arg = paste0("cp ",transPCO_wgt_dir,merge_three$gene[row],".wgt.RDat ",new_wgt_dir,merge_three$gene[row],".wgt.RDat ")
    } else if (merge_three$best_model[row] == 'r2.GBAT') {
        print('GBAT is the best')
        arg = paste0("cp ",GBAT_wgt_dir,merge_three$gene_full_GBAT[row],".wgt.RDat ",new_wgt_dir,merge_three$gene[row],".wgt.RDat ")
    }
    system(arg)

    gene_info = all_gene_info[grep(merge_three$gene[row],all_gene_info$geneId),]
    gene_chr = strsplit(gene_info$'#chrom',"chr")[[1]][2]
    gene_start = gene_info$chromStart
    gene_end = gene_info$chromEnd

    pos_matrix = rbind(pos_matrix,data.frame(WGT = paste0(merge_three$gene[row],".wgt.RDat"),ID = merge_three$gene[row],CHR = gene_chr,P0 = gene_start, P1 = gene_end))

}

pos_dir = "pos_files/Whole_Blood/"
dir.create(pos_dir,recursive = T)

fwrite(pos_matrix,paste0(pos_dir,"cis_GBAT.pos"),sep = '\t', quote = F, col.names = T, row.names = F)
