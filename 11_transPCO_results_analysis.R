library(data.table)
library(optparse)

# Define command-line options
option_list = list(
  make_option("--tissue", type = "character", default = "Whole_Blood", 
              help = "Tissue name to process", metavar = "character"),
  make_option("--pval_threshold", type = "numeric", default = 0.00001, 
              help = "P-value threshold to filter results", metavar = "numeric"),
  make_option("--output_dir", type = "character", default = "bed_files", 
              help = "Output file name for significant results", metavar = "character"),
  make_option("--PCO_association_dir",type = 'character',default = 'PCO_association_results',
	      help = "place where association results are stored"),
  make_option("--module_dir", type = "character", default = "modules",
              help = "dir containing modules", metavar = "character")
)

# Parse command-line arguments
opt = parse_args(OptionParser(option_list = option_list))

# Extract arguments
tissue = opt$tissue
pval_threshold = opt$pval_threshold
output_dir = opt$output_dir
PCO_association_dir = opt$PCO_association_dir

# Initialize results data.table
for (fold in 0:5) {
  sig_results = data.table(module = numeric(), rsIDs = character(), pval = numeric())
  # Determine number of modules
  module_dir = paste0("transPCO/", tissue, "/fold_", fold, "/",opt$module_dir,"/")
  modules = list.files(module_dir)

  PCO_association_path = paste0("transPCO/", tissue, "/fold_", fold, "/",opt$PCO_association_dir,"/")

  # Create bed directory
  bed_dir = paste0("transPCO/", tissue, "/fold_", fold, "/",output_dir,"/")
  dir.create(bed_dir, recursive = TRUE)

  for (module in modules) {

    module = strsplit(module,"\\.")[[1]][1]
    print(module) 
    for (chr in 1:22) {
      result_file = paste0(PCO_association_path, module, "_chr_", chr, ".txt.gz")
      # Skip if result file doesn't exist
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
    genes_in_module = fread(paste0(module_dir, module,".txt"), header = FALSE)
    genes_in_module = unlist(genes_in_module$V1)

    # Identify significant SNPs for each gene in this module
    for (gene in genes_in_module) {
      sig_snps_per_gene = as.matrix(sig_results[sig_results$module == module,c('rsIDs','pval')])

      if (length(sig_snps_per_gene) > 0) {
        fwrite(sig_snps_per_gene, 
               paste0(bed_dir, gene, ".txt"), 
               col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,append = TRUE)
      }
    }
  }
  fwrite(sig_results,paste0("transPCO/",tissue,"/fold_",fold,"/significant_results_",opt$module_dir,".txt"), col.names = T, row.names = F, sep = '\t', quote = F)
}


