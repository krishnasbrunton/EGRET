library("data.table")
library("optparse")

option_list = list(
        make_option("--gene_information", action="store",default=NA, type='character',
                help="Path to file containing gene information"),
	make_option("--expression_file", action="store",default=NA, type='character',
                help="Path to file containing gene expression"),
        make_option("--cis_window_size", action="store",default=1000000, type='numeric',
                help="size of cis window to create"),
        make_option("--tissue", action = "store", default = NA, type = 'character',
                help="tissue that will set up for"),
        make_option("--plink_path", action="store",default=NA, type='character',
                help="Path to plink"),
	make_option("--bfile", action="store",default=NA, type='character',
                help="path to location genotype files")
  )

opt = parse_args(OptionParser(option_list=option_list))

all_gene_info = fread(opt$gene_information,header = T)
expression = fread(opt$expression_file,header = T)

genes = unlist(expression$gene_id)

dir.create(paste0("bed_files/",opt$tissue,"/cis/"),recursive = T)
dir.create(paste0("plink_results/",opt$tissue,"/cis/"),recursive = T)

for (gene in genes) {
	gene_info = all_gene_info[grep(gene,all_gene_info$geneId)]
	chr = strsplit(gene_info$'#chrom',split = 'chr')[[1]][2]
	tss = gene_info$chromStart

	print(gene)
	if (chr == 'X' | chr == 'Y') {
  		next
	}

	bounds_u = max(1,tss - floor(opt$cis_window_size/2))
	bounds_l = tss + floor(opt$cis_window_size/2) 
	bed = c(chr,bounds_u,bounds_l,gene,0,"+")
	bed = matrix(bed, nrow = 1, ncol = 6)

	write.table(bed, paste0("bed_files/",opt$tissue,"/cis/",gene,".bed"), row.names = F, col.names = F, sep ="\t", quote = F)

	arg = paste0(opt$plink_path," --bfile ",opt$bfile," --extract bed1 bed_files/",opt$tissue,"/cis/",gene,".bed --make-bed --out plink_results/",opt$tissue,"/cis/",gene)
	system(arg)
	file.remove(paste0("plink_results/",opt$tissue,"/cis/",gene,".log"))
}

