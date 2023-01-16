#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH -o ./out/%x.out
#SBATCH -e ./err/%x.err
##############################################
# USAGE: sbatch --job-name=preprocessEncode4R --cpus-per-task=1 --mem-per-cpu=64 step1_preprocessEncode4R.sh $IN_FILE $OUT_DIR
# Date 01/12/2023
##############################################

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-R413

# Set-up dirs
IN_FILE=$1 # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/auxiliary_data/snrna/heart_Parse_10x_integrated.rds"
OUT_DIR=$2 # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/preprocess/snrna/subset"

# Print messages
echo -e "Preprocessing ENCODE4 data object: $IN_FILE"
echo -e "Saving preprocessed data in $OUT_DIR\n"

# Convert RDS matrix to Seurat object and preprocess
CMD="Rscript --vanilla createEncode4H5AD.R \
    $IN_FILE \
    $OUT_DIR"
echo -e "Running:\n $CMD\n"
$CMD

date
