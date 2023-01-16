#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o ./out/%x.out
#SBATCH -e ./err/%x.err
##############################################
# USAGE: sbatch --job-name=preprocess_encode_python --cpus-per-task=1 --mem-per-cpu=64 step1_preprocess_encode4_python.sh $IN_FILE $OUT_DIR
# Date 01/12/2023
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-py38

# Set-up dirs
IN_FILE=$1 # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_adrenal/preprocess/filtered.h5ad"
OUT_DIR=$2 # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_adrenal/preprocess/subset"

# Print messages
echo -e "Preprocessing ENCODE4 data object: $IN_FILE"
echo -e "Saving preprocessed data in $OUT_DIR\n"

# Create h5ad and loom files in Python
CMD="python create_encode4_loom.py \
    $IN_FILE \
    $OUT_DIR"
echo -e "Running:\n $CMD\n"
$CMD

date
