#########################################################################################
### Step 1 -- Create RNA count h5ad objects from Seurat object for all ENCODE 4 datasets

# Mouse adrenal tissue
IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/auxiliary_data/snrna/adrenal_Parse_10x_integrated.rds
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/prepare_inputs/snrna
sbatch \
    --job-name=prepare_inputsAdrenal4R \
    --cpus-per-task=1 \
    --mem-per-cpu=64G \
    step1_prepare_inputsEncode4R.sh \
    $IN_FILE $OUT_DIR

# Mouse heart tissue
IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/auxiliary_data/snrna/heart_Parse_10x_integrated.rds
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/prepare_inputs/snrna
sbatch \
    --job-name=prepare_inputsHeart4R \
    --cpus-per-task=1 \
    --mem-per-cpu=64G \
    step1_prepare_inputsEncode4R.sh \
    $IN_FILE $OUT_DIR

#########################################################################################

#########################################################################################
### Step 2 - Create RNA count loom objects from h5ad for all ENCODE 4 datasets
GENE_ARG=0.05
SAVE_CSV=False

# Mouse adrenal tissue
IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/prepare_inputs/snrna/filtered.h5ad
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/prepare_inputs/snrna
sbatch \
    --job-name=prepare_inputs_adrenal_python \
    --cpus-per-task=1 \
    --mem-per-cpu=128G \
    step2_prepare_inputs_encode4_python.sh \
    $IN_FILE $OUT_DIR

# Mouse heart tissue
IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/prepare_inputs/snrna/filtered.h5ad
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/prepare_inputs/snrna
sbatch --job-name=prepare_inputs_heart_python --cpus-per-task=1 --mem-per-cpu=128G step2_prepare_inputs_encode4_python.sh $IN_FILE $OUT_DIR
#########################################################################################
# Step 3 - Create a prepare_inputsed RNA Seurat object for the endothelial perturb seq dataset
IN_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seqprepare_inputs/
sbatch --job-name=prepare_inputsEndothelialR --cpus-per-task=1 --mem-per-cpu=256G step3_prepare_inputsEndothelialR.sh $IN_DIR $OUT_DIR

# Step 4 - Create a filtered loom file for the endothelial perturb seq dataset
IN_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seqprepare_inputs/
sbatch --job-name=prepare_inputs_endothelial_python --cpus-per-task=1 --mem-per-cpu=256G step4_prepare_inputs_endothelial_python.sh $IN_DIR $OUT_DIR
