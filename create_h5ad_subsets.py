# Script to subset an h5ad expression dataset and save different file formats
# Usage: python create_h5ad_subsets.py <h5ad_file> <out_dir> <gene_arg> <subset_columns>
# h5ad_file example: /cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/subset/filtered.h5ad
# out_dir example: /cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna
# gene_arg example: "0.05"
# subset_columns example: "celltypes, timepoint"

import argparse
import os
import sys

import loompy as lp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from utils import save_h5ad, save_loom

sc.settings.verbosity = 3
np.random.seed(13)

def make_dirs(path):
    if not os.path.exists(path):
        os.makedirs(path)

def qc(adata, out_dir, mito_filter=15, n_counts_filter=6000):
    # Funciton to plot QC metrics
    print(f"Calculating metrics and saving plot to {out_dir}/qc/qc.png...")

    # Another way of doing the exact same thing
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Total counts and mito filters, #draw horizontal red lines indicating thresholds.
    fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax = axs[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax = axs[1], show = False)
    axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
    axs[1].hlines(y = n_counts_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed')
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, "qc", "qc.png"))

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

def save_tsv(adata, out_dir, filename, normalize=False):
    # Function to save as tsv 
    if normalize:
        adata_copy = sc.pp.normalize_total(adata, target_sum=1e4, copy=True)
        sc.pp.log1p(adata_copy)
        print(f"Saving normalized tsv to {out_dir}/tsv/{filename}.normalized.tsv...")
        expr = pd.DataFrame(adata_copy.X.todense(), index=adata_copy.obs_names, columns=adata_copy.var_names).T
        expr.index.name = "gene"
        expr.to_csv(os.path.join(out_dir, "tsv", f"{filename}.normalized.tsv"), sep="\t")
    else:
        print(f"Saving counts tsv to {out_dir}/tsv/{filename}.counts.tsv...")
        expr = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names).T
        expr.index.name = "gene"
        expr.to_csv(os.path.join(out_dir, "tsv", f"{filename}.counts.tsv"), sep="\t")

def main(args):
    # Main script function

    # Parse args
    print("Parsing command line args...")
    h5ad_file = args.h5ad_file
    if "." in args.gene_arg:
        cell_frac = float(args.gene_arg)
        n_genes = None
        use_variable_genes = False
    else:
        cell_frac = None
        n_genes = int(args.gene_arg)
        use_variable_genes = True
    out_dir = os.path.join(args.out_dir, args.gene_arg)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    make_dirs(os.path.join(args.out_dir, "qc"))
    make_dirs(os.path.join(out_dir, "h5ad"))
    make_dirs(os.path.join(out_dir, "loom"))
    make_dirs(os.path.join(out_dir, "tsv"))
    save_tsv_flag = args.save_tsv
    subset_cols = [column.strip(",") for column in args.subset_columns]
    column_values = [value.strip(",") for value in args.column_values]
    if len(column_values) == 0:
        column_values = None
    all_cells_flag = args.skip_all_cells
    normalize_flag = args.normalize_tsv
    if column_values is not None:
        print(f"Will only keep cells with {column_values} in {subset_cols}...")
    else:
        print(f"Will keep all cells in {subset_cols}...")
    print(f"Will save tsv: {save_tsv_flag} with normalization: {normalize_flag}...")
    print(f"Will save all cells: {all_cells_flag}...")
    
    # Load in anndata
    print(f"Loading {h5ad_file}...")
    adata = sc.read_h5ad(h5ad_file)

    # Prelim filtters, takes a bout a minute
    print("Prelim filters using ScanPy min_genes=200 and min_cells=3...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.raw = adata

    # Plot QC
    qc(adata, args.out_dir)
    
    # Save object for All cells
    if all_cells_flag:
        curr_ad = subset_genes(adata, use_variable_genes=use_variable_genes, n_genes=n_genes, cell_frac=cell_frac, name="All")
    
        save_loom(curr_ad, out_dir, filename="All")
        save_h5ad(curr_ad, out_dir, filename="All")
        if save_tsv_flag:
            save_tsv(curr_ad, out_dir, filename="All", normalize=normalize_flag) 

    # Subset -- 500 cells, 2000 genes
    print("Subsetting genes for Subset to 500 cells, 2000 genes...")
    curr_ad = adata[adata.obs.index.isin(np.random.choice(adata.obs_names, 500)), adata.var.index.isin(np.random.choice(adata.var_names, 2000))]
    save_loom(curr_ad, out_dir, filename="Subset")
    save_h5ad(curr_ad, out_dir, filename="Subset")
    if save_tsv_flag:
        save_tsv(curr_ad, out_dir, filename="Subset", normalize=normalize_flag)

    # Column subsets
    print(f"Found {subset_cols} columns to subset on")
    for column in subset_cols:
        options = adata.obs[column].unique()
        print(f"Found {options} values for {column}")
        for value in options:
            print(f"Processing {column}={value}...")
            if column_values is not None and value not in column_values:
                continue
            curr_ad = adata[adata.obs[column] == value]
            curr_ad = subset_genes(curr_ad, use_variable_genes=use_variable_genes, n_genes=n_genes, cell_frac=cell_frac, name=value)
            save_loom(curr_ad, out_dir=out_dir, filename=value)
            save_h5ad(curr_ad, out_dir=out_dir, filename=value)
            if save_tsv_flag:
                save_tsv(curr_ad, out_dir, filename=value, normalize=normalize_flag)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad_file", help="h5ad file to subset", required=True)
    parser.add_argument("--out_dir", help="output directory to store all the prepared files, files will be saved in a subdirectory named after the gene_arg", required=True)
    parser.add_argument("--gene_arg", help="if an int, use that many variable genes, if a float, genes that are greater than that fraction of cells", required=True)
    parser.add_argument("--save_tsv", help="save tsv format of expression matrix, the reason for not doing this by default is that it takes a long time", action="store_true")
    parser.add_argument("--normalize_tsv", help="normalize the tsv expression matrix", action="store_true")
    parser.add_argument("--subset_columns", help="columns to slice the object by and save as separate files, if multiple they should be comma separated", nargs="+", default=[])
    parser.add_argument("--column_values", help="values to subset the columns by, if multiple they should be comma separated", nargs="+", default=[])
    parser.add_argument("--skip_all_cells", help="skip the All cells subset", action="store_false")
    args = parser.parse_args()
    main(args)