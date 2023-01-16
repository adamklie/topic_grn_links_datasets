# Step 1 -- Create RNA count h5ad objects from Seurat object for all ENCODE 4 datasets
IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/auxiliary_data/snrna/adrenal_Parse_10x_integrated.rds
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna
sbatch --job-name=preprocessAdrenal4R --cpus-per-task=1 --mem-per-cpu=64G step1_preprocessEncode4R.sh $IN_FILE $OUT_DIR

IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/auxiliary_data/snrna/heart_Parse_10x_integrated.rds
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/preprocess/snrna
sbatch --job-name=preprocessHeart4R --cpus-per-task=1 --mem-per-cpu=64G step1_preprocessEncode4R.sh $IN_FILE $OUT_DIR

# Step 2 - Create RNA count loom objects from h5ad for all ENCODE 4 datasets
IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna/filtered.h5ad
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_adrenal/preprocess/snrna
sbatch --job-name=preprocess_adrenal_python --cpus-per-task=1 --mem-per-cpu=128G step2_preprocess_encode4_python.sh $IN_FILE $OUT_DIR

IN_FILE=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/preprocess/snrna/filtered.h5ad
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/mouse_heart/preprocess/snrna
sbatch --job-name=preprocess_heart_python --cpus-per-task=1 --mem-per-cpu=128G step2_preprocess_encode4_python.sh $IN_FILE $OUT_DIR

# Step 3 - Create a preprocessed RNA Seurat object for the endothelial perturb seq dataset
IN_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq/preprocess
sbatch --job-name=preprocessEndothelialR --cpus-per-task=1 --mem-per-cpu=256G step3_preprocessEndothelialR.sh $IN_DIR $OUT_DIR

# Step 4 - Create a filtered loom file for the endothelial perturb seq dataset
IN_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq
OUT_DIR=/cellar/users/aklie/data/igvf/topic_grn_links/endothelial_perturb_seq/preprocess
sbatch --job-name=preprocess_endothelial_python --cpus-per-task=1 --mem-per-cpu=256G step4_preprocess_endothelial_python.sh $IN_DIR $OUT_DIR
