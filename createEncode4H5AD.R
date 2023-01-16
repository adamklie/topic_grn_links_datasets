suppressMessages(library(Seurat))
suppressMessages(library(loomR))
suppressMessages(library(PISCES))
suppressMessages(library(SeuratDisk))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
theme_set(theme_cowplot())

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)
IN_FILE <- args[1] #  "/cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/auxiliary_data/snrna/heart_Parse_10x_integrated.rds"
OUT_DIR <- args[2] # "/cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/preprocess/snrna/subset"

# Read in the R object
print("Reading...")
dat <- readRDS(IN_FILE)

# Grab only the RNA counts
print("Putting the object on a diet...")
DefaultAssay(dat) <- "RNA"
rna_dat <- DietSeurat(dat, assays=c("RNA"))

# Subset, comment out when ready
#print("Subsampling to 500 cells and saving...")
#n.cells = 500
#dat <- dat[, sample(colnames(dat), size=n.cells, replace=F)]
#rna_dat <- rna_dat[, sample(colnames(rna_dat), size=n.cells, replace=F)]
#saveRDS(dat, file.path(OUT_DIR, "normalized.rds"))
saveRDS(rna_dat, file.path(OUT_DIR, "filtered.rds"))

# Save as h5seuat
print("Saving h5seurat...")
SaveH5Seurat(rna_dat, filename = file.path(OUT_DIR, "filtered.h5Seurat"))

# Convert to h5ad
print("Saving h5ad...")
Convert(file.path(OUT_DIR, "filtered.h5Seurat"), dest = "h5ad")

# Grab the normalized data
norm.mat <- dat@assays$SCT@scale.data

# Save as tsv
print("Saving ARACNe tsv...")
mat.name <- file.path(OUT_DIR, 'normalized.tsv')
ARACNeTable(norm.mat, mat.name, subset = FALSE)
