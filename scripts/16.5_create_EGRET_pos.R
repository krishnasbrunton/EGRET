library(data.table)
library(optparse)

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

all_gene_info = fread("../data/GTEx_V8.txt.gz", header = T)
sumstats = fread(paste0("results_sumstats/",tissue,"/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1.txt"),header = T)
sumstats$best_r2_pval = pmin(sumstats$p_gw, sumstats$p_cis_part, sumstats$p_trans_part, na.rm = TRUE)
print(sumstats)
wgt_dir = paste0("xtune_fusion_models/",tissue,"/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1/")  

pos_matrix = data.frame(WGT = as.character(),ID = as.character(),CHR = as.numeric(), P0 = as.numeric(),P1 = as.numeric())

for (row in 1:nrow(sumstats)) {
    print(row)
    if (is.na(sumstats$'best_r2_pval'[row]) || sumstats$'best_r2_pval'[row] > 0.01) {
        next
    }

    gene_info = all_gene_info[grep(sumstats$gene[row],all_gene_info$geneId),]
    gene_chr = strsplit(gene_info$'#chrom',"chr")[[1]][2]
    gene_start = gene_info$chromStart
    gene_end = gene_info$chromEnd

    pos_matrix = rbind(pos_matrix,data.frame(WGT = paste0(sumstats$gene[row],".wgt.RDat"),ID = sumstats$gene[row],CHR = gene_chr,P0 = gene_start, P1 = gene_end))

}

pos_dir = paste0("pos_files/",tissue,"/")
dir.create(pos_dir,recursive = T)

fwrite(pos_matrix,paste0(pos_dir,"cis_MatrixeQTL_GBAT_transPCO_FDR_0.1.pos"),sep = '\t', quote = F, col.names = T, row.names = F)
