import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import os

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

for tissue in tissues:
    trans_connections = pd.read_csv(
        f"trans_connections_by_eqtl/{tissue}_trans_connections_by_eqtl.txt",
        sep="\t", header=0
    )

    trans_connections.loc[trans_connections['egene'].isna(),"egene"] = trans_connections['nearest_gene_id']
    trans_connections = trans_connections[trans_connections['egene'].notna()]

    # Extract unique genes for regulators and targets
    regulators = trans_connections["egene"].unique()
    targets = trans_connections["target_gene"].unique()

    # Initialize adjacency matrix (rows=regulators, cols=targets)
    trans_connection_matrix = pd.DataFrame(
        0, index=regulators, columns=targets, dtype=int
    )

    # Fill in matrix with 1 where connection exists
    for _, row in trans_connections.iterrows():
        trans_connection_matrix.loc[row["egene"], row["target_gene"]] = 1

    # Transpose so rows = target genes, cols = regulators
    G = trans_connection_matrix.T @ trans_connection_matrix

    regulator_to_coreg_genes = trans_connection_matrix @ G

    norm_sim = G.div(np.sqrt(np.outer(np.diag(G), np.diag(G))), axis=0)
    dist_matrix = 1 - norm_sim.fillna(0)

    # Hierarchical clustering
    condensed_dist = squareform(dist_matrix, checks=False)
    Z = linkage(condensed_dist, method="average")

    plt.figure(figsize=(6,4))
    dendrogram(Z, labels=G.index, leaf_rotation=90)
    plt.title("Hierarchical clustering of genes (based on shared regulators)")
    plt.savefig("grn_dendrogram.png",format = 'png')

    clusters = fcluster(Z, t=0.9, criterion='distance')
    gene_clusters = pd.DataFrame({"Gene": G.index, "Cluster": clusters})
    cluster_sizes = gene_clusters.groupby('Cluster').size()

    # Keep clusters with more than 1 gene
    multi_gene_clusters = cluster_sizes[cluster_sizes > 3].index

    output_dir = f"coregulation_grns/{tissue}"
    os.makedirs(output_dir, exist_ok=True)

    regulator_output_dir = f"coregulation_grn_regulators/{tissue}"
    os.makedirs(regulator_output_dir, exist_ok=True)

    # Get clusters with more than 3 genes
    cluster_sizes = gene_clusters.groupby('Cluster').size()
    multi_gene_clusters = cluster_sizes[cluster_sizes > 3].index

    # Filter and save each cluster separately
    for cluster_id in multi_gene_clusters:
        filtered_genes = gene_clusters[gene_clusters['Cluster'] == cluster_id]['Gene']
        
        submatrix = trans_connection_matrix.loc[:, filtered_genes]


        regulator_counts = submatrix.sum(axis=1).sort_values(ascending=False)
        top_regulators = regulator_counts[regulator_counts >= 2]

        out_path = os.path.join(output_dir, f"coregulation_cluster_{cluster_id}.txt")
        filtered_genes.to_csv(out_path, sep='\t', index=False, header=False)

        # Save the corresponding top regulators for this cluster
        reg_path = os.path.join(regulator_output_dir, f"coregulation_cluster_{cluster_id}_regulators.txt")
        top_regulators.to_csv(reg_path, sep='\t', header=['num_targets'])

        print(f"{tissue} | Cluster {cluster_id}: {len(filtered_genes)} genes, {len(top_regulators)} key regulators")


exit()
# Hierarchical clustering
#sns.clustermap(G, metric="cosine", method="average", cmap="Reds")

plt.savefig("grn_matrix_clustered_eqtl.png",format = 'png')

exit()
gene_info = pd.read_csv("../../data/GTEx_V8.txt.gz", header=0, sep="\t")


trans_connections = trans_connections.merge(
    gene_info[["geneId", "#chrom", "chromStart"]], 
    left_on="target_gene", right_on="geneId", 
    how="left"
)


chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
chrom_map = {chrom: i+1 for i, chrom in enumerate(chromosomes)}

# Assign numerical chromosome positions
trans_connections["target_gene_chr_num"] = trans_connections["#chrom"].map(chrom_map)


chromosome_lengths = {
1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
21: 46709983, 22: 50818468
}
chrom_offsets = {}
cumulative_offset = 0

for chrom, length in chromosome_lengths.items():
    chrom_offsets[chrom] = cumulative_offset
    cumulative_offset += length  # Move the offset for the next chromosome


trans_connections["target_gene_start_adjusted"] = (
    trans_connections["chromStart"] + trans_connections["target_gene_chr_num"].map(chrom_offsets)
)
trans_connections["trans_gene_start_adjusted"] = (
    trans_connections["snp_pos"] + trans_connections["snp_chr"].map(chrom_offsets)
)
x_min,x_max, y_min, y_max = 0,0,0,0
x_min, x_max = min(trans_connections["trans_gene_start_adjusted"].min(),x_min), max(trans_connections["trans_gene_start_adjusted"].max(),x_max)
y_min, y_max = min(y_min,trans_connections["target_gene_start_adjusted"].min()), max(y_max,trans_connections["target_gene_start_adjusted"].max())

fig, ax = plt.subplots(figsize=(10, 10))
scatter = ax.scatter(
    trans_connections["trans_gene_start_adjusted"], 
    trans_connections["target_gene_start_adjusted"], 
    alpha=0.4
)

legend_marker_size = 15


chrom_borders = list(chrom_offsets.values()) + [cumulative_offset]  # Add last chromosome end
chrom_centers = [(chrom_borders[i] + chrom_borders[i + 1]) / 2 for i in range(len(chrom_borders) - 1)]

# Set X and Y ticks at chromosome boundaries
ax.set_xticks(chrom_borders)
ax.set_xticklabels([""] * len(chrom_borders))  # Hide tick labels at borders

ax.set_yticks(chrom_borders)
ax.set_yticklabels([""] * len(chrom_borders))

# Place chromosome labels at their centers
for i, chrom in enumerate(chromosome_lengths.keys()):
    ax.text(chrom_centers[i], -cumulative_offset * 0.015, str(chrom), ha="center", fontsize=10)  # X-axis labels (bottom)
    ax.text(-cumulative_offset * 0.015, chrom_centers[i], str(chrom), va="center", fontsize=10, rotation=90)  # Y-axis labels (left)

# Set the axis limits to only cover the data range
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Grid at chromosome boundaries
ax.grid(True, linestyle="--", alpha=0.6)

ax.set_xlabel("Upstream Trans Gene Position")
ax.set_ylabel("Target Gene Position")
ax.set_title("Trans Regulators Plot")


plt.tight_layout()
plt.savefig("EGRET_trans_regulators.pdf", format="pdf", bbox_inches="tight")