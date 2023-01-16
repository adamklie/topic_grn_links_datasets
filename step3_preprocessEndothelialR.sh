#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o ./out/%x.out
#SBATCH -e ./err/%x.err
##############################################
# USAGE: sbatch --job-name=preprocessEndothelialR --cpus-per-task=1 --mem-per-cpu=256G step3_preprocessEndothelialR.sh
# Date 01/12/2023
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-R413

# Set-up dirs
IN_DIR=$1 # /cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq
OUT_DIR=$2 # /cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq/preprocess/subset

# Print messages
echo -e "Preprocessing endothelial cell data in $IN_DIR"
echo -e "Saving preprocessed data in $OUT_DIR\n"

# Convert RDS matrix to Seurat object and preprocess
CMD="Rscript --vanilla createEndothelialRDS.R \
    $IN_DIR \
    $OUT_DIR"
echo -e "Running:\n $CMD\n"
$CMD

date
