library(BoutrosLab.plotting.general);
library(Seurat);
library(pheatmap);
library(gtools);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/ghwei/objects');
source('/cluster/home/sujunc/chensj/scRNA/script/myfunctions/plot_functions.R');
source('/cluster/projects/hansengroup/sujunc/scRNA/script/myfunctions/plotcnv_functions.R');
name <- 'bph4';
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/ghwei/objects/p1234.combined.susu.rds');
seurat.all@meta.data$type <- 'Epi';
seurat.all@meta.data[seurat.all@ident%in%c(0, 9), ]$type <- 'LE';
seurat.all@meta.data[seurat.all@ident%in%c(1, 10, 11), ]$type <- 'SM'
seurat.all@meta.data[seurat.all@ident%in%c(2), ]$type <- 'Neu'
seurat.all@meta.data[seurat.all@ident%in%c(3, 14, 17), ]$type <- 'Fib'
seurat.all@meta.data[seurat.all@ident%in%c(6), ]$type <- 'T'
seurat.all@meta.data[seurat.all@ident%in%c(4, 12), ]$type <- 'Endo'
seurat.all@meta.data[seurat.all@ident%in%c(8), ]$type <- 'Macro'
seurat.all@meta.data[seurat.all@ident%in%c(16), ]$type <- 'B'
saveRDS(seurat.all, file = generate.filename('seurat', 'bph4', 'rds'));
pdf(generate.filename('plottsne', paste0(name, '_type'), 'pdf'));
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'type', do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = FALSE) + labs(x = 'Dimension 1', y = 'Dimension 2')
dev.off();
####
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/ghwei/objects/Pca_susu.rds');
name <- 'pca1';
seurat.all@meta.data$type <- 'Epi';
seurat.all@meta.data[seurat.all@ident%in%c(0, 1, 6), ]$type <- 'T';
seurat.all@meta.data[seurat.all@ident%in%c(2, 9), ]$type <- 'Endo'
seurat.all@meta.data[seurat.all@ident%in%c(7, 8, 13), ]$type <- 'CAF'
seurat.all@meta.data[seurat.all@ident%in%c(5), ]$type <- 'NK'
seurat.all@meta.data[seurat.all@ident%in%c(4, 10), ]$type <- 'B'
seurat.all@meta.data[seurat.all@ident%in%c(12), ]$type <- 'Mono'
seurat.all@meta.data[seurat.all@ident%in%c(11), ]$type <- 'Mast'
pdf(generate.filename('plottsne', paste0(name, '_type'), 'pdf'));
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'type', do.label = TRUE, label.size = 5, no.legend = TRUE, no.axes = FALSE) + labs(x = 'Dimension 1', y = 'Dimension 2')
dev.off();
saveRDS(seurat.all, file = generate.filename('seurat', 'pca1', 'rds'));
