library(data.table)
library(optparse)

# Load necessary functions and shared library
source("ModifiedPCOMerged.R")
source("ModifiedSigmaOEstimate.R")
source("davies.R")
source("liumod.R")
source("liu.R")
dyn.load("qfc.so")

# Define and parse command-line options
option_list = list(
  make_option("--tissue", action = "store", default = NA, type = 'character',
              help = "Tissue"),
  make_option("--module", action = "store", default = NA, type = 'character',
              help = "Module file name"),
  make_option("--module_dir", action = "store", default = NA, type = 'character',
              help = "Directory containing modules"),
  make_option("--association_dir", action = "store", default = NA, type = "character",
              help = "Directory with association results"),
  make_option("--results_dir", action = "store", default = NA, type = 'character',
              help = "Directory to put PCO association results")
)

opt = parse_args(OptionParser(option_list = option_list))
tissue = opt$tissue
module = opt$module
module_dir = opt$module_dir
association_dir = opt$association_dir
results_dir = opt$results_dir

# Read gene information and expression data upfront
gene_info = fread("../data/GTEx_V8.txt.gz", header = TRUE)
gene_expression = fread(paste0("expression_files/", tissue, "_expression_regressed.txt.gz"), header = TRUE)

# Process module across folds
for (fold in 0:5) {
  cat(sprintf("Processing Fold: %d, Module: %s\n", fold, module))

  # Define module path
  module_path = file.path("transPCO", tissue, paste0("fold_", fold), module_dir, module)
  if (!file.exists(module_path)) {
    cat(sprintf("Module file not found: %s\n", module_path))
    next
  }

  genes_in_module = fread(module_path, header = FALSE)

  # Define results directory
  fold_results_dir = file.path("transPCO", tissue, paste0("fold_", fold), results_dir)
  dir.create(fold_results_dir, recursive = TRUE, showWarnings = FALSE)

  module_name = sub("\\..*$", "", module)  # Remove file extension

  # Iterate over chromosomes 1–22
  for (chr in 1:22) {
    cat(sprintf("Processing Module: %s, Chromosome: %d, Fold: %d\n", module_name, chr, fold))

    # Filter genes in module by chromosome
    valid_genes = gene_info$geneId[gene_info$`#chrom` != paste0("chr", chr)]
    genes_in_module_filtered = genes_in_module[V1 %in% valid_genes, ]

    if (nrow(genes_in_module_filtered) == 0) {
      cat(sprintf("No valid genes for Module: %s, Chromosome: %d, Fold: %d\n", module_name, chr, fold))
      next
    }

    # Subset gene expression
    gene_expression_subset = gene_expression[gene_id %in% genes_in_module_filtered$V1, ]
    gene_expression_subset = gene_expression_subset[match(genes_in_module_filtered$V1, gene_expression_subset$gene_id), ]
    
    if (nrow(gene_expression_subset) == 0) {
      cat(sprintf("No gene expression data for Module: %s, Chromosome: %d, Fold: %d\n", module_name, chr, fold))
      next
    }

    gene_expression_subset = gene_expression_subset[, -1, with = FALSE]  # Remove gene_id column
    Sigma = cor(t(gene_expression_subset))  # Compute correlation matrix

    # Read Z matrix
    z_matrix_path = file.path("transPCO", tissue, paste0("fold_", fold), association_dir,
                              paste0("matrix_", module_name, "_chr_", chr, ".txt.gz"))

    if (!file.exists(z_matrix_path)) {
      cat(sprintf("Z matrix file not found: %s\n", z_matrix_path))
      next
    }

    Z.mat = fread(z_matrix_path, header = TRUE)
    snps = Z.mat$SNP
    rownames(Z.mat) = Z.mat$SNP  # Set SNP column as row names
    genes_in_module_vector = genes_in_module_filtered$V1

    # Subset Z.mat to include only columns matching genes in module
    Z.mat = Z.mat[, ..genes_in_module_vector, with = FALSE]

    # Estimate SigmaO and compute p-values
    SigmaO = ModifiedSigmaOEstimate(Sigma)
    pvals = ModifiedPCOMerged(Z.mat, Sigma, SigmaO, method = 'liu')

    # Save results
    output_file = file.path(fold_results_dir, paste0(module_name, "_chr_", chr, ".txt.gz"))
    pvals = data.table(rsIDs = snps, pval = pvals)
    fwrite(pvals, file = output_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    if (file.exists(paste0(results_dir, "/", module_name, "_chr_", chr, ".txt.gz"))) {
	file.remove(z_matrix_path)
    }

  }
}
