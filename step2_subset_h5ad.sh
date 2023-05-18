#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o ./out/%x.out
#SBATCH -e ./err/%x.err
##############################################
# USAGE:
# sbatch \
#   --job-name=subset_h5ad \
#   --cpus-per-task=1 \
#   --mem=128G \
#   $script \
#   $H5AD_FILE $OUT_DIR $GENE_ARG $SUBSET_COLUMNS
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scenicplus

# Set-up dirs
H5AD_FILE=$1 # /cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_adrenal/prepare_inputs/snrna/Converted.h5ad
OUT_DIR=$2 # /cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_adrenal/prepare_inputs/snrna/
GENE_ARG=$3 # 0.05
SUBSET_COLUMNS=$4 # "celltypes, timepoint"
SCRIPT=/cellar/users/aklie/projects/igvf/topic_grn_links/data/scripts/prepare_inputs_encode4_mouse.py

# Print messages
echo -e "Preparing inputs from ENCODE4 h5ad file: $H5AD_FILE"
echo -e "Saving inputs to folders in $OUT_DIR"
echo -e "Using gene argument: $GENE_ARG"
echo -e "No saving csv files"
echo -e "Subsetting columns: $SUBSET_COLUMNS"

# Create h5ad and loom files in Python
CMD="python $SCRIPT \
    --h5ad_file $H5AD_FILE \
    --out_dir $OUT_DIR \
    --gene_arg $GENE_ARG
    --subset_columns $SUBSET_COLUMNS"

echo "Running:\n $CMD\n"
$CMD

date
