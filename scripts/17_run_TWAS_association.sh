#!/bin/bash
#SBATCH --job-name="Running TWAS"
#SBATCH --output="batch_submissions/TWAS_output.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
home=/expanse/lustre/projects/ddp412/kbrunton

# Directory where pos files are located
#wgt=${home}/EGRET/pos_files/${tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1.pos
wgt=${home}/EGRET/pos_files/${tissue}/cis.pos

# Directory where weight files are located
wgtdir=${home}/EGRET/FUSION/${tissue}/cis/
#wgtdir=${home}/EGRET/xtune_fusion_models/${tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1/

# Directory to output TWAS results
#out=${home}/EGRET/TWAS_association_results/${tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1
out=${home}/EGRET/TWAS_association_results/${tissue}/cis

mkdir -p $out

ld=${home}/fusion_twas-master/LDREF/1000G.EUR.merged

#immune_sumstats=("PASS_IBD_deLange2017.sumstats.gz" "PASS_Rheumatoid_Arthritis.sumstats.gz" "UKB_460K.blood_PLATELET_COUNT.sumstats.gz" "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz" "UKB_460K.blood_RED_COUNT.sumstats.gz" "UKB_460K.blood_WHITE_COUNT.sumstats.gz" "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz" "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz")

#cvd_sumstats=("PASS_AtrialFibrillation_Nielsen2018.sumstats.gz" "PASS_HDL.sumstats.gz" "PASS_IschemicStroke_Malik2018.sumstats.gz" "PASS_LDL.sumstats.gz" "UKB_460K.biochemistry_Cholesterol.sumstats.gz" "UKB_460K.body_BMIz.sumstats.gz" "UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz")

#neuro_sumstats=("PASS_ADHD_Demontis2018.sumstats.gz" "PASS_Alzheimers_deRojas2021.sumstats.gz" "PASS_Autism.sumstats.gz" "PASS_Schizophrenia_Pardinas2018.sumstats.gz" "PASS_AnorexiaNervosa_Watson2019.sumstats.gz")

general_sumstats=("PASS_IBD_deLange2017.sumstats.gz" "PASS_Rheumatoid_Arthritis.sumstats.gz" "UKB_460K.blood_PLATELET_COUNT.sumstats.gz" "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz" "UKB_460K.blood_RED_COUNT.sumstats.gz" "UKB_460K.blood_WHITE_COUNT.sumstats.gz" "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz" "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz" "PASS_AtrialFibrillation_Nielsen2018.sumstats.gz" "PASS_HDL.sumstats.gz" "PASS_IschemicStroke_Malik2018.sumstats.gz" "PASS_LDL.sumstats.gz" "UKB_460K.biochemistry_Cholesterol.sumstats.gz" "UKB_460K.body_BMIz.sumstats.gz" "UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz" "PASS_ADHD_Demontis2018.sumstats.gz" "PASS_Alzheimers_deRojas2021.sumstats.gz" "PASS_Autism.sumstats.gz" "PASS_Schizophrenia_Pardinas2018.sumstats.gz" "PASS_AnorexiaNervosa_Watson2019.sumstats.gz")

#general_sumstats=("UKB_460K.blood_WHITE_COUNT.sumstats.gz")

chr=1  #doesn't matter if using trans script
#for gwas in $(ls ../TCSC/sumstats/ | grep gz)
#do

for gwas in "${general_sumstats[@]}"
do
	
	echo $gwas
	trait="${gwas::-12}"	
	output=$trait
	gwas=${home}/TCSC/sumstats/$gwas
	if [ -f "${out}/${output}.dat" ]; then
		    continue
	fi
	Rscript ${home}/fusion_twas-master/FUSION.assoc_test_trans.R --chr $chr --ref_ld_chr $ld --sumstats $gwas --weights $wgt --weights_dir $wgtdir --out ${out}/${output}.dat
	echo "done with "$trait
done
