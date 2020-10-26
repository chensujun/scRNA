library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
name <- 'all'
genes.marker <- c('ACOX2')
pdf(generate.filename('plotgene', paste0(name, '_ACOX2'), 'pdf'), width = 5, height = 5);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 1, nCol = 1, min.cutoff = 0,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();

genes.marker <- c('PECAM1', 'PDGFRB', 'ITGAV', 'ITGB3');
pdf(generate.filename('plotgene', paste0(name, '_markers'), 'pdf'), width = 5, height = 5);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 1, nCol = 2, min.cutoff = 0,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
