library(BoutrosLab.plotting.general);
library(Seurat);
library(scran);
library(scater);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/seurat/");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
newColDat <- readRDS(conf$colData);
seurat.all@meta.data$type <- newColDat$type;
####
cd8 <- seurat.all@data['CD8A', ];
seurat.cd8 <- SubsetData(object = seurat.all, cells.use = rownames(seurat.all@meta.data[seurat.all@meta.data$type%in%c('T')|cd8>0, ]));
