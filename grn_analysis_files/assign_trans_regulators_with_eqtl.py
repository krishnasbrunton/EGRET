import pandas as pd 
from convert_variant_id_to_rsid import convert_variant_id_to_rsid

tissues = [
    "Adipose_Visceral_Omentum", "Artery_Tibial", "Brain_Cortex", "Colon_Transverse",
    "Esophagus_Mucosa", "Heart_Left_Ventricle", "Lung", "Muscle_Skeletal",
    "Skin_Not_Sun_Exposed_Suprapubic", "Whole_Blood"
]

tissues = [
    "Adipose_Subcutaneous","Adipose_Visceral_Omentum","Artery_Tibial","Brain_Cortex","Colon_Transverse",
    "Esophagus_Mucosa","Heart_Left_Ventricle","Lung","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic",
    "Whole_Blood","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue",
    "Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid",
    "Esophagus_Gastroesophageal_Junction","Esophagus_Muscularis","Heart_Atrial_Appendage",
    "Kidney_Cortex","Liver","Minor_Salivary_Gland","Nerve_Tibial","Ovary","Pancreas","Pituitary",
    "Prostate","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen",
    "Stomach","Testis","Thyroid","Uterus","Vagina"
]
tissues = ['Whole_Blood']


for tissue in tissues:
    print(tissue)
    convert_variant_id_to_rsid(tissue)
    # Load data
    trans_connections = pd.read_csv(
        f"tissue_trans_connections_updated/{tissue}_trans_connections.txt",
        sep='\t'
    )
    eqtl_sumstats = pd.read_csv(
        f"GTEx_Analysis_v8_eQTL/{tissue}_with_rsid.txt",
        sep='\t'
    )

    # Merge on 'rsid'
    merged = trans_connections.merge(
        eqtl_sumstats[['rsid', 'gene_id']],
        left_on='snp',
        right_on = 'rsid',
        how='left'
    )

    # Add gene_id as new column 'egene'
    merged.rename(columns={'gene_id': 'egene'}, inplace=True)

    merged.loc[merged['egene'].isna(),"egene"] = merged['nearest_gene_id']
    merged.to_csv(f"trans_connections_by_eqtl/{tissue}_trans_connections_by_eqtl.txt", sep="\t", index=False)

    print(len(merged['target_gene'].unique()))
    print((len(merged['egene'].unique())))
    print(len(merged['snp'].unique()))
