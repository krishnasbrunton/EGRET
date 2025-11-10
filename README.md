# EGRET (Estimating Genome-wide Regulatory Effects on the Transcriptome)

EGRET is a multivariate linear model designed to identify genome-wide loci that are predictive of gene expression levels. EGRET integrates predictions from existing trans-eQTL mapping approaches (Matrix eQTL, GBAT, and trans-PCO) and determines the optimal weighted combination of regulatory variants that best explains gene expression. Finally, EGRET utilizes a genome-wide summary statistics-based TWAS to identify novel gene-disease associations.

## Setting up the environment
Install dependences
- create yml file

R packages
- data.table
- optparse
## Setting up genotype files 
The script 0_setup_genotypes.R takes an input plink formatted genotype files and preprocesses them for EGRET. This script does two main things: prune SNPs by LD and removes rare variants (MAF < 0.05).

```
Rscript 0_setup_genotypes.R --bfile {unprocessed plink files} --LD_r2 {LD threshold used for pruning SNPs} --LD_window {window size} --LD_chunk {chunk size} --out {output path} --plink_path {path to plink}
```
