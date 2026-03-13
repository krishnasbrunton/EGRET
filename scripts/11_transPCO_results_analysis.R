library(data.table)
library(optparse)

# Define command-line options
option_list = list(
  make_option("--tissue", type = "character", default = "Whole_Blood",
              help = "Tissue name to process", metavar = "character"),
  make_option("--pval_threshold", type = "numeric", default = 0.00001,
              help = "P-value threshold to filter results", metavar = "numeric"),
  make_option("--output_dir", type = "character", default = NA,
              help = "Base output directory for the pipeline", metavar = "character"),
  make_option("--bed_dir", type = "character", default = "bed_files",
              help = "Subdirectory name for output bed files", metavar = "character"),
  make_option("--PCO_association_dir", type = 'character', default = 'PCO_association_results',
              help = "Subdirectory name where PCO association results are stored"),
  make_option("--module_dir", type = "character", default = "modules",
              help = "Subdirectory name containing modules", metavar = "character"),
  make_option("--folds", type = "integer", default = 5,
              help = "Number of cross-validation folds (iterates fold_0 through fold_N)")
)

# Parse command-line arguments
opt = parse_args(OptionParser(option_list = option_list))

tissue = opt$tissue
pval_threshold = opt$pval_threshold
base_dir = opt$output_dir
PCO_association_dir = opt$PCO_association_dir

# Initialize results data.table
for (fold in 0:opt$folds) {
  sig_results = data.table(module = numeric(), rsIDs = character(), pval = numeric())

  module_dir_path = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), opt$module_dir)
  modules = list.files(module_dir_path)

  PCO_association_path = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), PCO_association_dir)

  # Create bed directory
  bed_dir_path = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), opt$bed_dir)
  dir.create(bed_dir_path, recursive = TRUE)

  for (module in modules) {

    module = strsplit(module, "\\.")[[1]][1]
    print(module)
    for (chr in 1:22) {
      result_file = file.path(PCO_association_path, paste0(module, "_chr_", chr, ".txt.gz"))
      print(result_file)
      if (!file.exists(result_file)) {
        next
      }

      # Read results and filter by p-value
      results = fread(result_file, header = TRUE)
      sig_rows = results[results$pval < pval_threshold, ]
      if (nrow(sig_rows) > 0) {
        sig_results = rbind(sig_results, data.table(module = module, sig_rows))
      }
    }
    print(sig_results)

    # Read genes in the current module
    genes_in_module = fread(file.path(module_dir_path, paste0(module, ".txt")), header = FALSE)
    genes_in_module = unlist(genes_in_module$V1)

    # Identify significant SNPs for each gene in this module
    for (gene in genes_in_module) {
      sig_snps_per_gene = as.matrix(sig_results[sig_results$module == module, c('rsIDs', 'pval')])

      if (length(sig_snps_per_gene) > 0) {
        fwrite(sig_snps_per_gene,
               file.path(bed_dir_path, paste0(gene, ".txt")),
               col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE, append = TRUE)
      }
    }
  }
  fwrite(sig_results,
         file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), paste0("significant_results_", opt$module_dir, ".txt")),
         col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
}
