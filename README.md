# EGRET (Estimating Genome-wide Regulatory Effects on the Transcriptome)

EGRET is a multivariate linear model designed to identify genome-wide loci that are predictive of gene expression levels. EGRET integrates predictions from existing trans-eQTL mapping approaches (Matrix eQTL, GBAT, and trans-PCO) and determines the optimal weighted combination of regulatory variants that best explains gene expression. Finally, EGRET utilizes a genome-wide summary statistics-based TWAS to identify novel gene-disease associations.

## Setting up the environment
Install dependences
- create yml file

R packages
Required
-   data.table
-   optparse
-   plink2R
-   Matrix eQTL
  
Recomended
-   xtune-lasso

## Setup
To setup expression and genotype files in the format that EGRET requires, we have created a wrapper script to carry out preprocessing of genoytype and expression data. All outputs are stored in the ouput_dir.

| Argument                 | Description                                                                                                   |
| ------------------------ | ------------------------------------------------------------------------------------------------------------- |
| `genotypes_file_path`    | Path to PLINK `.bed/.bim/.fam` files containing genotype data.                                                |
| `expression_file_path`   | Gene expression matrix (rows = genes, columns = individuals). Must include a `gene_id` column in ENSG format. |
| `covariates_file_path`   | Covariate file (rows = individuals, columns = covariates to regress out).                                     |
| `tissue`                 | Name of tissue for analysis. Used for labeling and output organization.                                       |
| `individuals_file_path`  | File containing list of individual IDs (one per line, no header).                                             |
| `LD_prune_r2`            | LD pruning threshold (default: `0.9`).                                                                        |
| `plink_path`             | Path to PLINK2 executable.                                                                                    |
| `genotype_output_prefix` | Prefix for pruned genotype output files.                                                                      |
| `folds`                  | Number of cross-validation folds (e.g., 5).                                                                   |
| `gene_info_file_path`    | File containing gene metadata (gene ID, name, chromosome, start position).                                    |
| `output_dir`             | Directory to store all processed outputs.                                                                     |

```
./setup_genotype_and_expression.sh \
    $genotypes_file_path \
    $expression_file_path \
    $covariates_file_path \
    $tissue \
    $individuals_file_path \
    $LD_prune_r2 \
    $plink_path \
    $genotype_output_prefix \
    $folds \
    $gene_info_file_path \
    $output_dir
```

## Running Matrix eQTL, GBAT, and _trans_-PCO
In EGRET, we implement Matrix eQTL, GBAT, and _trans_-PCO methods to identify _trans_-variants with potential to enhance gene expression models. Each one is run separatly and the outputs of each is used in the final training of the genome-wide model.

### Running Matrix eQTL
We utilize the R package, Matrix eQTL, to conduct pairwise association between all genetic variants and gene expression. As input this script takes the tissue you wish to carry this analysis on. It is assumed that the setup scripts were run so it looks for their respective outputs.

```
./MatrixeQTL_scripts.sh \
    $tissue \
    $FDR \
    $output_dir
```

### Running GBAT
GBAT utilizes previously trained cis-expression models to find genome-wide regulators. This code is designed to take as input models fitted by FUSION. 
* If FUSION models have not already been generated, model weights can be downloaded from http://gusevlab.org/projects/fusion/

```
./GBAT_scripts.sh \
  


```




