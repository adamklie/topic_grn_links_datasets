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
rds_file <- args[1] #  /cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/auxiliary_data/snrna/heart_Parse_10x_integrated.rds
out_dir <- args[2] # /cellar/users/aklie/projects/igvf/topic_grn_links/data/mouse_heart/prepare_inputs/snrna/

# Read in the R object
print(sprintf("Reading in %s...", rds_file))
dat <- readRDS(rds_file)

# Grab only the RNA counts
print("Putting the object on a diet...")
DefaultAssay(dat) <- "RNA"
rna_dat <- DietSeurat(dat, assays=c("RNA"))

# Save as h5seuat
print("Saving h5seurat...")
SaveH5Seurat(rna_dat, filename = file.path(out_dir, "tmp", "Converted.h5Seurat"))

# Convert to h5ad
print("Saving h5ad...")
Convert(file.path(out_dir, "tmp", "Converted.h5Seurat"), dest = "h5ad")

# Move this out of tmp
print("Moving h5ad...")
file.rename(file.path(out_dir, "tmp", "Converted.h5ad"), file.path(out_dir, "Converted.h5ad"))

# Subset, comment out when ready
print("Subsampling to 500 cells and saving...")
n.cells = 500
n.genes = 2000
subset_dat <- dat[sample(rownames(dat), size=n.genes, replace=F), sample(colnames(dat), size=n.cells, replace=F)]
subset_rna_dat <- rna_dat[sample(rownames(rna_dat), size=n.genes, replace=F), sample(colnames(rna_dat), size=n.cells, replace=F)]
saveRDS(dat, file.path(out_dir, "normalized.rds"))
saveRDS(rna_dat, file.path(out_dir, "filtered.rds"))

# Save as h5seuat
print("Saving h5seurat...")
SaveH5Seurat(rna_dat, filename = file.path(out_dir, "tmp", "Converted.h5Seurat"))

# Convert to h5ad
print("Saving h5ad...")
Convert(file.path(out_dir, "tmp", "Converted.h5Seurat"), dest = "h5ad")

# Grab the normalized data
norm.mat <- dat@assays$SCT@scale.data

# Save as tsv
print("Saving ARACNe tsv...")
mat.name <- file.path(out_dir, 'normalized.tsv')
ARACNeTable(norm.mat, mat.name, subset = FALSE)
