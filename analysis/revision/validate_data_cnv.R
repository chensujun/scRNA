library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision');
name <- 'all';
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
mytype <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_celltype_scnorm.txt');
pdf(generate.filename('plottsne', paste0(name, '_type2'), 'pdf'), width = 6, height = 6)
DimPlot(seurat.all, reduction.use = 'tsne', group.by = 'type_scnorm', pt.size = 0.5, do.label = TRUE, vector.friendly = FALSE, no.axes = TRUE, no.legend = TRUE) ;
dev.off();
####
mywidth <- 16;
myheight <- 12;
myres <- 300;
#### F1A
tiff(generate.filename('plottsne', paste0(name, '_type'), 'tiff'), width = mywidth, height = mywidth, res = myres, units = 'in')
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'type_scnorm', do.label = TRUE, label.size = 24, no.legend = FALSE, no.axes = FALSE) 
dev.off();
