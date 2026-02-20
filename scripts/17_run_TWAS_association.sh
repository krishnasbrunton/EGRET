#!/bin/bash
#SBATCH --job-name="Running TWAS"
#SBATCH --output="/expanse/lustre/scratch/kbrunton/temp_project/batch_submission_output/TWAS_output/%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 48:00:00

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

general_sumstats=("UKB_460K.cov_EDU_YEARS.sumstats.gz" "UKB_460K.lung_FVCzSMOKE.sumstats.gz" "PASS_Intelligence_SavageJansen2018.sumstats.gz" "UKB_460K.mental_NEUROTICISM.sumstats.gz" "UKB_460K.body_HEIGHTz.sumstats.gz" "UKB_460K.other_MORNINGPERSON.sumstats.gz" "UKB_460K.repro_MENARCHE_AGE.sumstats.gz" "UKB_460K.biochemistry_TotalProtein.sumstats.gz" "UKB_460K.body_WHRadjBMIz.sumstats.gz" "UKB_460K.lung_FEV1FVCzSMOKE.sumstats.gz" "PASS_ReactionTime_Davies2018.sumstats.gz" "PASS_Myopia_Hysi2020.sumstats.gz" "UKB_460K.biochemistry_Creatinine.sumstats.gz" "UKB_460K.bmd_HEEL_TSCOREz.sumstats.gz" "PASS_SleepDuration_Dashti2019.sumstats.gz" "UKB_460K.biochemistry_IGF1.sumstats.gz" "PASS_GeneralRiskTolerance_KarlssonLinner2019.sumstats.gz" "PASS_Eosino_Vuckovic2020.sumstats.gz" "UKB_460K.biochemistry_AspartateAminotransferase.sumstats.gz" "PASS_Insomnia_Jansen2019.sumstats.gz" "PASS_Height1.sumstats.gz" "PASS_LymP_Vuckovic2020.sumstats.gz" "UKB_460K.repro_NumberChildrenEverBorn_Pooled.sumstats.gz" "UKB_460K.biochemistry_Testosterone_Male.sumstats.gz" "PASS_BMI1.sumstats.gz" "PASS_AgeFirstBirth.sumstats.gz" "PASS_Glaucoma_Craig2020.sumstats.gz" "PASS_RTC_Vuckovic2020.sumstats.gz" "UKB_460K.body_BALDING1.sumstats.gz" "UKB_460K.biochemistry_AlkalinePhosphatase.sumstats.gz" "PASS_DrinksPerWeek_Liu2019.sumstats.gz" "UKB_460K.biochemistry_Phosphate.sumstats.gz" "PASS_MedicationUse_Wu2019.sumstats.gz" "UKB_460K.biochemistry_VitaminD.sumstats.gz" "PASS_MDD_Wray2018.sumstats.gz" "PASS_RD_Zhao2021.sumstats.gz" "PASS_MonoP_Vuckovic2020.sumstats.gz" "UKB_460K.biochemistry_TotalBilirubin.sumstats.gz" "UKB_460K.repro_MENOPAUSE_AGE.sumstats.gz" "PASS_SA_Grasby2020.sumstats.gz" "PASS_HipOA_Tachmazidou2019.sumstats.gz" "UKB_460K.pigment_SUNBURN.sumstats.gz" "PASS_NumberChildrenEverBorn.sumstats.gz" "PASS_BipolarDisorder_Ruderfer2018.sumstats.gz" "PASS_Years_of_Education1.sumstats.gz" "PASS_PancreasVol_Liu2021.sumstats.gz" "PASS_CaudateVol_Satizabal2019.sumstats.gz" "PASS_BrainstemVol_Satizabal2019.sumstats.gz" "PASS_TH_Grasby2020.sumstats.gz" "PASS_BasoP_Vuckovic2020.sumstats.gz" "PASS_Anorexia.sumstats.gz" "PASS_Ever_Smoked.sumstats.gz" "PASS_MO_Zhao2021.sumstats.gz" "PASS_ProstateCancer.sumstats.gz" "UKB_460K.cancer_PROSTATE.sumstats.gz" "PASS_AccumbensVol_Satizabal2019.sumstats.gz" "UKB_460K.cancer_BREAST.sumstats.gz" "PASS_Type_2_Diabetes.sumstats.gz" "PASS_BipolarDisorder_Ruderfer2018.sumstats.gz")


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
