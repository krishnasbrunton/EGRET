library(data.table)

all_gene_info = fread("../data/GTEx_V8.txt.gz",header = T)
cross_mappable_genes = fread("../mappability_encode/cross_mappability_strength.txt.gz",header = F)

num_mismatches_allowed = 0
cross_mappable_genes =cross_mappable_genes[cross_mappable_genes$V3 > num_mismatches_allowed]
focal_genes = unique(cross_mappable_genes$V2)

output_dir = paste0("bed_files/cross_mapped/",num_mismatches_allowed,"_mismatches/")
dir.create(output_dir,recursive = T)

for (gene in focal_genes) {
	#if (file.exists(paste0(output_dir,gene,".bed"))) { next }
	print(gene)
	subset = cross_mappable_genes[cross_mappable_genes$V2 == gene,]
	print(sum(subset$V3)/20000)
	next
	subset = subset[subset$V3 > num_mismatches_allowed,]
	focal_gene_info = all_gene_info[all_gene_info$geneId %in% gene,]
	crossmap_gene_info = all_gene_info[all_gene_info$geneId %in% subset$V1,c('#chrom','chromStart','geneId')]
	
	crossmap_gene_info = crossmap_gene_info[crossmap_gene_info$'#chrom' != focal_gene_info$'#chrom' | crossmap_gene_info$chromStart < focal_gene_info$chromStart - 1000000 | crossmap_gene_info$chromStart > focal_gene_info$chromStart + 1000000, ]



	crossmap_gene_info$chr = sapply(strsplit(as.character(crossmap_gene_info$'#chrom'), "chr"), function(x) x[2])
	u_bounds = pmax(0,crossmap_gene_info$chromStart - 100000)
	l_bounds = crossmap_gene_info$chromStart + 100000
	crossmap_bed = cbind(crossmap_gene_info$chr,u_bounds,l_bounds,crossmap_gene_info$geneId, "0","+")
	fwrite(crossmap_bed,paste0(output_dir,gene,".bed"),quote = F, col.names = F, row.names = F, sep = '\t')
}
