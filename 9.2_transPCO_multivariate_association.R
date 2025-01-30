library(data.table)
library(optparse)
# Load necessary functions and shared library
#source("ModifiedPCOMerged_acat_version.R")
#source("ModifiedPCOMerged.R")
source("ModifiedPCOMerged.R")
source("ModifiedSigmaOEstimate.R")
source("davies.R")
source("liumod.R")
source("liu.R")
dyn.load("qfc.so")

# Define and parse command-line options
option_list = list(
  make_option("--fold", action = "store", default = NA, type = 'character',
              help = "Fold to run analysis on"),
  make_option("--tissue", action = "store", default = NA, type = 'character',
              help = "Tissue"),
  make_option("--module_dir", action = "store", default = NA, type = 'character',
              help = "Directory containing the modules to run on"),
  make_option("--association_dir", action = "store", default = NA, type = "character",
	      help = "Directory with association results"),
  make_option("--results_dir", action = "store", default = NA, type = 'character',
              help = "Directory to put PCO association results")
  )

opt = parse_args(OptionParser(option_list = option_list))
fold = opt$fold
tissue = opt$tissue
module_dir = opt$module_dir
results_dir = opt$results_dir
association_dir = opt$association_dir

# Define commonly used directories
base_dir = paste0("transPCO/", tissue, "/fold_", fold)
results_dir = paste0(base_dir, "/", results_dir)
modules_dir = paste0(base_dir, "/", module_dir)

# Create results directory if it doesn't exist
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Read gene information and expression data upfront
gene_info = fread("../data/GTEx_V8.txt.gz", header = TRUE)
gene_expression = fread("expression_files/Whole_Blood_expression_regressed.txt.gz", header = TRUE)

# Get module files
modules = list.files(modules_dir)
print(modules)
# Loop through each module and chromosome
for (module in modules) {
  module_name = sub("\\..*$", "", module)  # Remove file extension for module name
  module_name_short = paste(strsplit(module_name,"_")[[1]][1:2],collapse = "_")
  print(module_name_short)
  all_genes_in_module = fread(file.path(modules_dir, module), header = FALSE)

  for (chr in 1:22) {
  #for (chr in c(1,19)) {	
    cat(sprintf("Processing Module: %s, Chromosome: %d\n", module_name, chr))
    
    # Filter genes in module by chromosome
    valid_genes = gene_info$geneId[gene_info$`#chrom` != paste0("chr", chr)]
    genes_in_module = all_genes_in_module[V1 %in% valid_genes, ]
    if (nrow(genes_in_module) == 0) {
      cat(sprintf("No valid genes for Module: %s, Chromosome: %d\n", module_name, chr))
      next
    }
    
    # Subset gene expression
    gene_expression_subset = gene_expression[gene_id %in% genes_in_module$V1, ]
    gene_expression_subset = gene_expression_subset[match(genes_in_module$V1, gene_expression_subset$gene_id), ]
    if (nrow(gene_expression_subset) == 0) {
      cat(sprintf("No gene expression data for Module: %s, Chromosome: %d\n", module_name, chr))
      next
    }
    gene_expression_subset = gene_expression_subset[, -1, with = FALSE]  # Remove gene_id column
    Sigma = cor(t(gene_expression_subset))  # Compute correlation matrix

    # Read Z matrix
    z_matrix_path = paste0(base_dir, "/",association_dir,"/matrix_", module_name, "_chr_", chr, ".txt.gz")
    if (!file.exists(z_matrix_path)) {
      cat(sprintf("Z matrix file not found: %s\n", z_matrix_path))
      next
    }
    Z.mat = fread(z_matrix_path, header = TRUE)
    snps = Z.mat$SNP
    rownames(Z.mat) = Z.mat$SNP  # Set SNP column as row names
    genes_in_module_vector = genes_in_module$V1
    
    # Subset Z.mat to include only columns matching genes in module
    Z.mat = Z.mat[, ..genes_in_module_vector, with = FALSE]
    #Z.mat = Z.mat[1:10,]
    #print(Z.mat)
    #Z.mat = transpose(Z.mat)    
    # Estimate SigmaO and compute p-values
    SigmaO = ModifiedSigmaOEstimate(Sigma)
    pvals = ModifiedPCOMerged(Z.mat, Sigma, SigmaO,method = 'liu')
    #pvals = ModifiedPCOMerged_acat(Z.mat,Sigma)

    # Prepare and save results
    pvals = data.table(rsIDs = snps, pval = pvals)
    fwrite(pvals, file = paste0(results_dir, "/", module_name, "_chr_", chr, ".txt.gz"),
           sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


