import pandas as pd
import numpy as np
import pyreadr
from scipy.spatial import KDTree
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tissue', help='tissue for analysis')
args = parser.parse_args()
tissue = args.tissue

sumstats_path = f"results_sumstats/{tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1.txt"
sumstats = pd.read_csv(sumstats_path, header = 0, sep = '\t')


column_mapping = {
    'r2_gw': 'gw',
    'r2_cis_part': 'cis',
    'r2_trans_part': 'trans'
}

sumstats['r2_gw'] = sumstats['r2_gw'].fillna(0)
sumstats['r2_cis_part'] = sumstats['r2_cis_part'].fillna(0)
sumstats['r2_trans_part'] = sumstats['r2_trans_part'].fillna(0)

sumstats['best_model'] = sumstats[['r2_gw', 'r2_cis_part', 'r2_trans_part']].idxmax(axis=1).map(column_mapping)

sumstats['trans_component'] = (
    (sumstats['p_trans_part'] < 0.01) & 
    ((sumstats['best_model'] == 'gw') | (sumstats['best_model'] == 'trans'))
)

tissue_genes_path = f"expression_files/{tissue}_expression.txt.gz"
tissue_genes = pd.read_csv(tissue_genes_path, header = 0, sep = '\t')

# Load gene annotation (GTEx format assumed)
gene_annot = pd.read_csv("../data/GTEx_V8.txt.gz", sep='\t')
gene_annot = gene_annot[['geneId', 'name', '#chrom', 'chromStart', 'chromEnd']]
gene_annot['#chrom'] = gene_annot['#chrom'].str.replace('chr', '').astype(str)
gene_annot = gene_annot[gene_annot['geneId'].isin(tissue_genes['gene_id'])]

# Build a KDTree for gene positions by chromosome
gene_kdtrees = {}
for chrom in gene_annot['#chrom'].unique():
    positions = gene_annot[gene_annot['#chrom'] == chrom][['chromStart']].values
    gene_kdtrees[chrom] = KDTree(positions)

results = []

for _, row in sumstats[sumstats['trans_component']].iterrows():
    gene = row['gene']
    best_model = row['best_model']
    gene_chr = gene_annot.loc[gene_annot['geneId'] == gene, '#chrom'].values[0]
    
    # Load RDat file
    rdat_path = f"xtune_fusion_models/{tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1/{gene}.wgt.RDat"
    print(rdat_path)
    try:
        rdat = pyreadr.read_r(rdat_path)
    except Exception as e:
        print(f"Failed to read {rdat_path}: {e}")
        continue

    print(rdat['cv.performance'])
    model = rdat['cv.performance'].idxmax(axis=1).iloc[0]
    model_weights  = rdat['wgt.matrix'][model]
    model_weights = model_weights[model_weights != 0]

    # Get the SNP metadata, aligning by SNP ID
    snp_info = rdat['snps']
    snp_info.index = rdat['wgt.matrix'].index  # Align index with weights

    # Get the chromosome per SNP
    snp_chr = snp_info.loc[model_weights.index, 'V1'].astype(str)
    trans_weights = model_weights[snp_chr != gene_chr]
    

    for snp_id, weight in trans_weights.items():
        snp_row = snp_info.loc[snp_id]
        snp_chr = str(snp_row['V1'])
        snp_pos = snp_row['V4'] 

        if snp_chr in gene_kdtrees:
            dist, idx = gene_kdtrees[snp_chr].query([[snp_pos]])
            nearest_gene = gene_annot[(gene_annot['#chrom'] == snp_chr)].iloc[idx[0]]

            results.append({
                'target_gene': gene,
                'snp': snp_id,
                'snp_chr': snp_chr,
                'snp_pos': snp_pos,
                'nearest_gene': nearest_gene['name'],
                'nearest_gene_id': nearest_gene['geneId'],
                'distance': dist[0],
                'weight': weight
            })

trans_sources = pd.DataFrame(results)
print(trans_sources.head())
output_path = f"grn_analysis_files/tissue_trans_connections_updated/{tissue}_trans_connections.txt"
trans_sources.to_csv(output_path, sep="\t", index=False)
