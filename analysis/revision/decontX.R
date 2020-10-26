library(Seurat);
library(celda);
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
decontxModel <- decontX(counts = seurat.all@raw.data, z = seurat.all@meta.data$cluster);
