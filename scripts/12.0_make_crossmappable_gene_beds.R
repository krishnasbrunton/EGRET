library(data.table)
library(optparse)

option_list = list(
  make_option("--gene_info", action="store", default=NA, type='character',
              help="Path to gene info file (e.g. GTEx_V8.txt.gz)"),
  make_option("--crossmap_file", action="store", default=NA, type='character',
              help="Path to cross-mappability strength file (cross_mappability_strength.txt.gz)"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Base output directory; BED files written to output_dir/cross_mapped/background_mismatches/"),
  make_option("--cis_window", action="store", default=1000000, type='integer',
              help="Window (bp) around focal gene TSS used to exclude same-region cross-map partners [default: 1000000]"),
  make_option("--crossmap_window", action="store", default=100000, type='integer',
              help="Window (bp) around each cross-mapped gene TSS for BED region bounds [default: 100000]")
)

opt = parse_args(OptionParser(option_list=option_list))

# gene info
all_gene_info = fread(opt$gene_info, header = TRUE)

# cross mappability downloaded from saha et al.
cross_mappable_genes = fread(opt$crossmap_file, header = FALSE)

focal_genes = unique(cross_mappable_genes$V2)

out_dir = file.path(opt$output_dir, "cross_mapped", "background_mismatches")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
num_cross_mapped_genes = length(focal_genes)

for (gene in focal_genes) {
  print(gene)
  gene_rows = cross_mappable_genes[cross_mappable_genes$V2 == gene, ]
  background_mappability = sum(gene_rows$V3) / num_cross_mapped_genes

  print(background_mappability)
  gene_rows = gene_rows[gene_rows$V3 > background_mappability, ]
  if (nrow(gene_rows) == 0) next

  focal_gene_info = all_gene_info[all_gene_info$geneId == gene, ]
  if (nrow(focal_gene_info) == 0) next

  crossmap_gene_info = all_gene_info[all_gene_info$geneId %in% gene_rows$V1, c('#chrom', 'chromStart', 'geneId')]

  crossmap_gene_info = crossmap_gene_info[
    crossmap_gene_info$'#chrom' != focal_gene_info$'#chrom' |
    crossmap_gene_info$chromStart < focal_gene_info$chromStart - opt$cis_window |
    crossmap_gene_info$chromStart > focal_gene_info$chromStart + opt$cis_window, ]

  if (nrow(crossmap_gene_info) == 0) next

  crossmap_gene_info$chr = sapply(strsplit(as.character(crossmap_gene_info$'#chrom'), "chr"), function(x) x[2])

  bed_start = pmax(0, crossmap_gene_info$chromStart - opt$crossmap_window)
  bed_end   = crossmap_gene_info$chromStart + opt$crossmap_window
  crossmap_bed = cbind(crossmap_gene_info$chr, bed_start, bed_end, crossmap_gene_info$geneId, "0", "+")
  fwrite(crossmap_bed, file.path(out_dir, paste0(gene, ".bed")), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
}
