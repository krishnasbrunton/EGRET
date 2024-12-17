#!/bin/bash
#SBATCH --job-name="Running TWAS"
#SBATCH --output="batch_submissions/TWAS_output.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 12:00:00

home=/expanse/lustre/projects/ddp412/kbrunton

#gwas=${home}/TCSC/sumstats/PASS_ADHD_Demontis2018.sumstats.gz

wgt=${home}/full_workflow/pos_files/Whole_Blood/cis.pos
wgtdir=${home}/full_workflow/xtune_fusion_models/Whole_Blood/cis/
#wgtdir=${home}/transeQTLs/weights/Whole_Blood_gtex_tf_chip_abc_nocis_concatenated
#wgtdir=${home}/transeQTLs/weights/Whole_Blood_gtex_allsnps_basil_concatenated

out=${home}/full_workflow/TWAS_assocation_results/Whole_Blood/cis
mkdir -p $out

ld=${home}/fusion_twas-master/LDREF/1000G.EUR.merged

immune_sumstats=("PASS_IBD_deLange2017.sumstats.gz" "PASS_Rheumatoid_Arthritis.sumstats.gz" "UKB_460K.blood_PLATELET_COUNT.sumstats.gz" "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz" "UKB_460K.blood_RED_COUNT.sumstats.gz" "UKB_460K.blood_WHITE_COUNT.sumstats.gz" "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz" "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz")

#cvd_sumstats=("PASS_AtrialFibrillation_Nielsen2018.sumstats.gz" "PASS_HDL.sumstats.gz" "PASS_IschemicStroke_Malik2018.sumstats.gz" "PASS_LDL.sumstats.gz" "UKB_460K.biochemistry_Cholesterol.sumstats.gz" "UKB_460K.body_BMIz.sumstats.gz" "UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz")

chr=1  #doesn't matter if using trans script
#for gwas in $(ls ../TCSC/sumstats/ | grep gz)
#do

for gwas in "${immune_sumstats[@]}"
do
echo $gwas
trait="${gwas::-12}"	
output=$trait
gwas=${home}/TCSC/sumstats/$gwas
#gwas=${home}/TCSC/sumstats/PASS_ADHD_Demontis2018.sumstats.gz
Rscript ${home}/fusion_twas-master/FUSION.assoc_test_trans.R --chr $chr --ref_ld_chr $ld --sumstats $gwas --weights $wgt --weights_dir $wgtdir --out ${out}/${output}.dat
echo "done with "$trait
done
