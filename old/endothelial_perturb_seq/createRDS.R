suppressMessages(library(Seurat))
suppressMessages(library(WGCNA))
suppressMessages(library(PISCES))
suppressMessages(library(loomR))
suppressMessages(library(SeuratDisk))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
theme_set(theme_cowplot())

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)
IN_DIR <- args[1]
OUT_DIR <- args[2]

# Read in the R object
dat <- readRDS(file.path(IN_DIR, "aggregated.2kG.library.mtx.cell_x_gene.RDS"))

# Transpose matrix using WGCNA
t_dat <- transposeBigData(x = dat)
rm(dat)

# Cell metadata
cell_dat <- read.table(file.path(IN_DIR, "barcodes.expanded.df.txt"), sep = "\t", header = TRUE, row.names = 1)

# Feature metadata
var_dat0 <- read.table(file.path(IN_DIR, "features.tsv"), sep = "\t", header = FALSE)
var_dat <- head(data.frame(do.call('rbind', strsplit(as.character(var_dat0$V1), ':', fixed=TRUE))))
colnames(var_dat) <- c("gene_name", "gene_id")

# Create object
seurat_obj <- CreateSeuratObject(counts=t_dat, project="endothelial_perturb_seq", assay="RNA", meta.data = cell_dat, min.cells = 3, min.features = 200)

# Subset, comment out when ready
n.cells = 500
seurat_obj <- seurat_obj[, sample(colnames(seurat_obj), size=n.cells, replace=F)]

# Save the raw object
saveRDS(seurat_obj, file.path(OUT_DIR, "raw.rds"))

# Filter cells
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < 6000 & percent.mt < 15)

# Save the filtered object
saveRDS(seurat_obj, file.path(OUT_DIR, "filtered.rds"))

# Preprocess
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Save the raw object
saveRDS(seurat_obj, file.path(OUT_DIR, "normalized.rds"))

# Grab the normalized data
norm.mat <- subset_obj@assays$RNA@data

# Save as tsv
mat.name <-  file.path(OUT_DIR, 'normalized.tsv')
ARACNeTable(data.frame(norm.mat), mat.name, subset = FALSE)

