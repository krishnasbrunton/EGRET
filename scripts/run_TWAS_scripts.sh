#!/bin/bash
tissue=$1
trait=$2
gwas_sumstat_dir=$3
$output_dir=$3


home=/expanse/lustre/projects/ddp412/kbrunton

Rscript 16.5_create_EGRET_pos.R \
    -tissue $tissue

ld=${home}/fusion_twas-master/LDREF/1000G.EUR.merged

# run TWAS
Rscript ${home}/fusion_twas-master/FUSION.assoc_test_trans.R \
    --chr $chr \
    --ref_ld_chr $ld \
    --sumstats $gwas \
    --weights $wgt \
    --weights_dir $wgtdir \
    --out ${out}/${output}.dat


    