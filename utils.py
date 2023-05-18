import os

import loompy as lp
import numpy as np
import scanpy as sc


def make_dirs(path):
    if not os.path.exists(path):
        os.makedirs(path)

def save_h5ad(adata, out_dir, filename="All"):
    # Function to save as h5ad
    print(f"Saving as h5ad to {out_dir}/h5ad/{filename}.h5ad...")
    adata.write(os.path.join(out_dir, "h5ad", f"{filename}.h5ad"))

def save_loom(adata, out_dir, filename="All"):
    # Function to save as loom
    print(f"Saving as loom to {out_dir}/loom/{filename}.loom...")
    row_attrs = dict(zip(adata.var.reset_index().columns, adata.var.reset_index().values.T))
    col_attrs = dict(zip(adata.obs.reset_index().columns, adata.obs.reset_index().values.T))
    row_attrs["Gene"] = np.array(adata.var_names)
    col_attrs["CellID"] = np.array(adata.obs_names)
    col_attrs["nGene"] = np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten()
    col_attrs["nUMI"] = np.array(np.sum(adata.X.transpose(), axis=0)).flatten()
    lp.create(os.path.join(out_dir, "loom", f"{filename}.loom"), adata.X.transpose(), row_attrs, col_attrs)


def subset_genes(adata, use_variable_genes=True, n_genes=None, cell_frac=None, name="All"):
    # Function to subset genes using either variable genes or cell fraction
    print(f"Subsetting genes for {name} with use_variable_genes={use_variable_genes}, n_genes={n_genes}, cell_frac={cell_frac}...")
    if use_variable_genes:
        assert n_genes is not None
        adata_copy = sc.pp.normalize_total(adata, target_sum=1e4, copy=True)
        sc.pp.log1p(adata_copy)
        sc.pp.highly_variable_genes(adata_copy, n_top_genes=n_genes, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata_copy.var.highly_variable]
    else:
        assert cell_frac is not None
        sc.pp.filter_genes(adata, min_cells=adata.shape[0]*cell_frac)
    return adata