#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o ./out/%x.out
#SBATCH -e ./err/%x.err
##############################################
# USAGE: sbatch --job-name=preprocess_endothelial_python --cpus-per-task=1 --mem-per-cpu=256G step4_preprocess_endothelial_python.sh
# Date 01/12/2023
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-py38

# Set-up dirs
IN_DIR=$1 # /cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq
OUT_DIR=$2 # /cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq/preprocess/subset

# Print messages
echo -e "Preprocessing endothelial cell data in $IN_DIR"
echo -e "Saving preprocessed data in $OUT_DIR\n"

# Create h5ad and loom files in Python
CMD="python create_endothelial_h5ad.py \
    $IN_DIR \
    $OUT_DIR"
echo -e "Running:\n $CMD\n"
$CMD

date
