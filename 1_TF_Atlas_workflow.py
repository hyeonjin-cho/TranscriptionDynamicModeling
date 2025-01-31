import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np

# Load an H5AD file
adata = sc.read_h5ad("GSE217460_210322_TFAtlas_raw.h5ad")

# Step 1: filter out all the genes that have low expression
raw_counts = adata.X
cell_names = adata.obs_names
gene_names = adata.var_names

raw_counts_df = pd.DataFrame(
    raw_counts,  # Already a NumPy array
    index=cell_names,
    columns=gene_names
).T

gene_sums = raw_counts_df.sum(axis=1)
raw_counts_df = raw_counts_df[gene_sums > 0]

gene_sums = raw_counts_df.sum(axis=1)
threshold = np.percentile(gene_sums, 10)
filtered_matrix = raw_counts_df[gene_sums > threshold]
filtered_genes = gene_sums[gene_sums > threshold].index
adata_subset = adata[:, adata.var_names.isin(filtered_genes)]

with open("gene_pool_Pawel_sorted.txt", "r") as f:
	decay_rate_genes_list = [line.strip() for line in f]

adata_subset = adata_subset[:, adata_subset.var_names.isin(decay_rate_genes_list)]


# Step 2: Find the indices where the 'TF' column contains TF of interest
'''
GFP: TFORF3549
FOXA1: TFORF3365
GATA2: TFORF2997|TFORF2998
SOX9: TFORF2531
GR: TFORF3432|TFORF1432|TFORF1433

TBP: TFORF1771|TFORF1771
TFIIA: TFORF1840|TFORF3076|TFORF3238
TFIIB: TFORF3401 
TFIIE: TFORF1086
TFIIF: TFORF0832|TFORF0833 
TFIIH: TFORF1445|TFORF1446|TFORF3186|TFORF0779|TFORF0780 
'''

indices = adata_subset.obs[adata_subset.obs['TF'].str.contains('TFORF1445|TFORF1446|TFORF3186|TFORF0779|TFORF0780')].index

# Step 3: Subset the AnnData object using these indices
adata_tf_subset = adata_subset[indices, :]

filtered_df = adata_tf_subset.X
cell_names = adata_tf_subset.obs_names
gene_names = adata_tf_subset.var_names

filtered_df = pd.DataFrame(
    filtered_df,  # Already a NumPy array
    index=cell_names,
    columns=gene_names
).T

filtered_df.to_csv("TFIIH/adata_subset_transposed_rawCounts_filtered_TFIIH.csv")
