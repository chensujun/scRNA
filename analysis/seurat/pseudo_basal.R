library(BoutrosLab.plotting.general);
library(Seurat);
library(scran);
library(scater);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/seurat/");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
newColDat <- readRDS(conf$colData);
seurat.all@meta.data$type <- newColDat$type;
markers <- readRDS(conf$marker_format);
#### find variable genes
seurat.154 <- SubsetData(object = seurat.all, cells.use = rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident == 'JD1800154SL'&seurat.all@meta.data$type%in%c('Basal/intermediate', 'Luminal'), ]));
seurat.154 <- FindVariableGenes(seurat.154)
seurat.154 <- ScaleData(seurat.154);
### chose genes with higher than median dispersiom and expression
#seurat.154 <- RunPCA(seurat.154, pc.genes = rownames((seurat.154@hvg.info[seurat.154@hvg.info$gene.dispersion>1.0769&seurat.154@hvg.info$gene.mean>0.11373, ])));
seurat.154 <- RunPCA(seurat.154);
###
pdf(generate.filename('all154', 'pca', 'pdf'), width = 8);
PCAPlot(seurat.154, dims = 1:2, reduction = "pca");
dev.off();
seurat.154 <- FindClusters(seurat.154, resolution = 0.8);
pdf(generate.filename('all154', 'cluster', 'pdf'), width = 8);
DimPlot(seurat.154, dims = 1:2, reduction = "pca");
dev.off();
pdf(generate.filename('all154', 'markers', 'pdf'), width = 8);
FeaturePlot(seurat.154, features = markers[['basal_intermediate']], reduction = "pca");
DimPlot(seurat.154, group.by = 'type', reduction = "pca");
dev.off();
