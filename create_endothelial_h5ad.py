import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import loompy as lp
sc.settings.verbosity = 3

# Get args
args = sys.argv
IN_DIR = args[1] # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/endothelial_perturb_seq"
OUT_DIR = args[2] # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/endothelial_perturb_seq/preprocess/subset"

# Define files
MTX_FILE = os.path.join(IN_DIR, "sparse.aggregated.2kG.library.cell_x_gene.mtx.gz")
VAR_FILE = os.path.join(IN_DIR, "features.tsv")
OBS_FILE = os.path.join(IN_DIR, "barcodes.expanded.df.txt")

# Load in matrix take ~30m on my machine
adata = sc.read_mtx(MTX_FILE)

# Read in the barcode information and what gene was targeted
obs_df = pd.read_csv(OBS_FILE, delimiter="\t", index_col=0)
obs_df.head()

# Add feature data
var_df = pd.read_csv(VAR_FILE, delimiter="\t", header=None)
var_df[["gene_name", "gene_id"]] = [row for row in var_df[0].str.split(":")]
var_df = var_df.drop(0, axis=1).set_index("gene_name")
var_df.head()

# Add these to anndata
adata.obs = obs_df
adata.var = var_df

# Make var names unique
adata.var_names_make_unique()

# Create a subset, comment out when ready
#np.random.seed(13)
#adata = sc.pp.subsample(adata, n_obs=500, copy=True, random_state=13)
#adata = adata[:, adata.var.index.isin(np.random.choice(adata.var_names, 10000))]

# Write raw h5ad
adata.write(os.path.join(OUT_DIR, "raw.h5ad"))

# Prelim filtters, takes a bout a minute
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# How many mitochondrial genes do we have in our feature set
adata.var_names.str.startswith('MT-').sum()

# Another way of doing the exact same thing
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Perform the actual filtering
mito_filter = 15
n_counts_filter = 6000
adata = adata[adata.obs.n_genes_by_counts < n_counts_filter, :]
adata = adata[adata.obs.pct_counts_mt < mito_filter, :]

# Save the raw for downstream transforms
adata.raw = adata

# Write raw h5ad
adata.write(os.path.join(OUT_DIR, "filtered.h5ad"))

# Save as loom
row_attrs = dict(zip(adata.var.reset_index().columns, adata.var.reset_index().values.T))
col_attrs = dict(zip(adata.obs.reset_index().columns, adata.obs.reset_index().values.T))
row_attrs["Gene"] = np.array(adata.var_names)
col_attrs["CellID"] = np.array(adata.obs_names)
col_attrs["nGene"] = np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten()
col_attrs["nUMI"] = np.array(np.sum(adata.X.transpose(), axis=0)).flatten()
lp.create(os.path.join(OUT_DIR, "filtered.loom"), adata.X.transpose(), row_attrs, col_attrs)
