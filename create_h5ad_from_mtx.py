# Script to prepare inputs for encode4 mouse data
# Usage: python prepare_inputs_encode4_mouse.py <h5ad_file> <out_dir> <gene_arg> <subset_columns>
# h5ad_file example: "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/subset/filtered.h5ad"
# out_dir example: "/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna"
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

def main(args):

    # Define files
    MTX_FILE = args.mtx_file
    OBS_FILE = args.obs_file
    VAR_FILE = args.var_file
    OUT_DIR = args.out_dir

    # Load in matrix
    print(f"Reading in {MTX_FILE} using scanpy read_mtx...")
    adata = sc.read_mtx(MTX_FILE)

    # Read in the barcode information
    print(f"Reading in {OBS_FILE}...")
    obs_df = pd.read_csv(OBS_FILE, delimiter="\t", index_col=0)
    obs_df.head()
    adata.obs = obs_df

    # Add feature data
    print(f"Reading in {VAR_FILE}...")
    var_df = pd.read_csv(VAR_FILE, delimiter="\t", header=None)
    var_df[["gene_name", "gene_id"]] = [row for row in var_df[0].str.split(":")]
    var_df = var_df.drop(0, axis=1).set_index("gene_name")
    var_df.head()
    adata.var = var_df
    adata.var_names_make_unique()

    # Write raw h5ad
    save_h5ad(adata, OUT_DIR, filename="All.raw")

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
    save_h5ad(adata, OUT_DIR, filename="All.filtered")

    # Save as loom
    save_loom(adata, OUT_DIR, filename="All.filtered")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create ScanPy compatible h5ad and loom files from mtx, obs, and var files")
    parser.add_argument("--mtx_file", help="matrix market file to read in", required=True)
    parser.add_argument("--obs_file", help="tsv file with cell barcodes and metadata", required=True)
    parser.add_argument("--var_file", help="tsv file with gene names and metadata", required=True)
    parser.add_argument("--out_dir", help="directory to save output files", required=True)
    args = parser.parse_args()
    main(args)