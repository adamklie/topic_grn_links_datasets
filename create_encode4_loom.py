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
H5AD_FILE = args[1] # "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/subset/filtered.h5ad"
OUT_DIR = args[2] # "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/subset"

# Load in anndata
print("Reading...")
adata = sc.read_h5ad(H5AD_FILE)

# Create a subset, comment out when ready
#np.random.seed(13)
#adata = sc.pp.subsample(adata, n_obs=500, copy=True, random_state=13)
#adata = adata[:, adata.var.index.isin(np.random.choice(adata.var_names, 10000))]

# Prelim filtters, takes a bout a minute
print("Prelim filters...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# How many mitochondrial genes do we have in our feature set
print("Calculate metrics...")
adata.var_names.str.startswith('mt-').sum()

# Another way of doing the exact same thing
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Save the raw for downstream transforms
adata.raw = adata

# Save as loom
print("Saving as loom...")
row_attrs = dict(zip(adata.var.reset_index().columns, adata.var.reset_index().values.T))
col_attrs = dict(zip(adata.obs.reset_index().columns, adata.obs.reset_index().values.T))
row_attrs["Gene"] = np.array(adata.var_names)
col_attrs["CellID"] = np.array(adata.obs_names)
col_attrs["nGene"] = np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten()
col_attrs["nUMI"] = np.array(np.sum(adata.X.transpose(), axis=0)).flatten()
lp.create(os.path.join(OUT_DIR, "filtered.loom"), adata.X.transpose(), row_attrs, col_attrs)
