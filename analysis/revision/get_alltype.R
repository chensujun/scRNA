library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
seurat.e <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/endo/objects/Manual_seurat_endo.rds');
seurat.f <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/fibroblast/objects/13_sample_fibroblast_20190718/GraphClust.seuset.rds');
seurat.m <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/macrophage/objects/Manual_seurat_macro.rds')
pdf(generate.filename('plottsne', paste0('mono', '_cluster'), 'pdf'), width = 6, height = 6);
DimPlot(seurat.m, reduction.use = 'tsne', pt.size = 1,  do.label = TRUE, vector.friendly = TRUE, label.size = 24, no.legend = TRUE, no.axes = FALSE) +
 labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();

mytype <- seurat.all@meta.data[, 'type', drop = FALSE];
seurat.t@meta.data$type <- 'CD8_effector';
seurat.t@meta.data[seurat.t@meta.data$res.0.8%in%c(1,4), ]$type <- 'CD8_naive';
seurat.t@meta.data[seurat.t@meta.data$res.0.8%in%c(0), ]$type <- 'CD4_conv'; 
seurat.t@meta.data[seurat.t@meta.data$res.0.8%in%c(6), ]$type <- 'CD4_Treg';
###
seurat.e@meta.data$type <- 'aEC';
seurat.e@meta.data[seurat.e@meta.data$res.0.8%in%c(0, 1), ]$type <- 'endothelial';
###
seurat.f@meta.data$type <- 'fib_S1';
seurat.f@meta.data[seurat.f@meta.data$res.0.8%in%c(2), ]$type <- 'fib_S2';
seurat.f@meta.data[seurat.f@meta.data$res.0.8%in%c(4), ]$type <- 'fib_S3';
###
seurat.m@meta.data$type <- 'TAM';
seurat.m@meta.data[seurat.m@meta.data$res.0.8%in%c(3), ]$type <- 'monocyte';
seurat.m@meta.data[seurat.m@meta.data$res.0.8%in%c(1), ]$type <- 'DC';
###
mytype[rownames(seurat.t@meta.data), ] <- seurat.t@meta.data$type;
mytype[rownames(seurat.e@meta.data), ] <- seurat.e@meta.data$type;
mytype[rownames(seurat.f@meta.data), ] <- seurat.f@meta.data$type;
mytype[rownames(seurat.m@meta.data), ] <- seurat.m@meta.data$type;
mytype$type_ori <- seurat.all@meta.data$type;
mytype <- mytype[!mytype$type%in%c('T', 'Myeloid', 'Endothelial', 'Fibroblast'), ];
saveRDS(mytype, generate.filename('celltype', 'all', 'rds'));
