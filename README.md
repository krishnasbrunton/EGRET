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
We recommend using these setup scripts in order to ensure the formatting of your input files for EGRET are consistent with what the scripts expect. We also do some preprocessing in order to enhance computational efficiency.

## Setting up genotype files 
The script 0_setup_genotypes.R takes an input plink formatted genotype files and preprocesses them for EGRET. This script does two main things: prune SNPs by LD and removes rare variants (MAF < 0.05).

```
Rscript 0_setup_genotypes.R \
  --bfile {unprocessed plink files} \
  --LD_r2 {LD threshold used for pruning SNPs} \
  --LD_window {window size} \
  --LD_chunk {chunk size} \
  --out {output path} \
  --plink_path {path to plink}
```

## Setting up expression files
The script 1_setup_expression.R takes as input expression (already normalized), tissue, covariates (for example sex, age, genotype PCs, gene expression PCs), gene information (such as gene type, location), and individuals to create a processed expression output. This script selects expression from only the individual ids specified. Next, we remove genes which are not protein coding, anti-sense RNA, or linc RNAs. Next the script regresses the provided covariates. The script outputs both processed expression (no covariate regression) and expression after covariate regression.

```
Rscript 1_setup_expression.R \
  --expression {expression file} \
  --tissue {tissue} \
  --individuals {individual ids to keep} \
  --out {output file name} \
  --covariates {covariate file}
```

## Setting up folds
The script 2_setup_folds.R divides individuals among folds prior to any model training. This essential for cross-validation R^2. Since many of the feature selection methods are lengthy, we have opted to store individuals for each fold and their respective gene expression and results. As a default we recommend 5-fold cross-validation, however, we have also provided the option to create more folds.

| Option          | Type        | Default    | Description                                                                                                |
| :-------------- | :---------- | :--------- | :--------------------------------------------------------------------------------------------------------- |
| `--expression`  | `character` | *Required* | Path to gene expression matrix (genes × individuals).                                                      |
| `--individuals` | `character` | *Required* | File containing list of individuals to include. Used to ensure consistent sample ordering across datasets. |
| `--tissue`      | `character` | *Required* | Tissue name used for labeling or directory structure.                                                      |
| `--covariates`  | `character` | *Optional* | Path to covariate file (e.g., PEER factors, genotype PCs). Must correspond to the same individuals.        |
| `--folds`       | `numeric`   | `5`        | Number of cross-validation folds to create.                                                                |
| `--plink_path`  | `character` | *Optional* | Path to the PLINK executable (used for genotype subsetting).                                               |
| `--bfile`       | `character` | *Optional* | Path prefix to PLINK-format genotype files (.bed/.bim/.fam).                                               |


```
Rscript 2_setup_folds.R \
  --expression path/to/expression_matrix.txt \
  --individuals path/to/individual_list.txt \
  --tissue "Whole_Blood" \
  --covariates path/to/covariates.txt \
  --folds 5 \
  --plink_path /path/to/plink \
  --bfile path/to/genotypes
```

## Running Matrix eQTL, GBAT, and trans-PCO
In EGRET, we implement Matrix eQTL, GBAT, and _trans_-PCO methods to identify _trans_-variants with potential to enhance gene expression models. Each one is run separatly and the outputs of each is used in the final training of the genome-wide model.

### Running Matrix eQTL
We utilize the R package, Matrix eQTL, to conduct pairwise association between all genetic variants and gene expression. As input this script takes the tissue you wish to carry this analysis on. It is assumed that the setup scripts were run so it looks for their respective outputs.

```
Rscript 7_MatrixeQTL.R --tissue {tissue}
```

### Running GBAT





