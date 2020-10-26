library(BoutrosLab.plotting.general);
library(Seurat);
library(readxl);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/normal');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/fibroblast/objects/13_sample_fibroblast_20190718/GraphClust.seuset.rds');
name <- 'caf';
genes.marker <- c('PTGS2', 'NCOA7', 'IGKC', 'NR4A2', 'A2M', 'BGN')
pdf(generate.filename('plotgene', paste0(name, '_1'), 'pdf'), width = 15, height = 9);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 4, nCol = 3, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();

genes.marker <- c('CD74', 'CFH', 'HSPE1', 'TPM4', 'TSPYL2')
pdf(generate.filename('plotgene', paste0(name, '_2'), 'pdf'), width = 15, height = 9);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 4, nCol = 3, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
###
genes <- data.frame(read_excel('Cluster4_receptor.xlsx'));
genes <- c(intersect(genes$allmarkers.gene, rownames(seurat.all@data)), 'IGHA1','IGHG3','IGHG1','IGHG4','IGLC2','IGKC','IGFBP2');
genes.exp <- rowMeans(seurat.all@data[genes, ])
genes.marker <- names(genes.exp[genes.exp<1]) 
pdf(generate.filename('plotviolin', paste0(name, '_marker_normal_1'), 'pdf'), width = 16, height = 16);
VlnPlot(seurat.all, features.plot = genes[1:40], point.size.use = 0, nCol = 5);
dev.off();
pdf(generate.filename('plotviolin', paste0(name, '_marker_normal_2'), 'pdf'), width = 16, height = 16);
VlnPlot(seurat.all, features.plot = genes[41:80], point.size.use = 0, nCol = 5);
dev.off();
pdf(generate.filename('plotviolin', paste0(name, '_marker_normal_3'), 'pdf'), width = 16, height = 16);
VlnPlot(seurat.all, features.plot = genes[81:104], point.size.use = 0, nCol = 5);
dev.off();
###
pdf(generate.filename('plotgene', paste0(name, '_marker_normal_1'), 'pdf'), width = 16, height = 16);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes[1:30], pt.size = 4, nCol = 5, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
pdf(generate.filename('plotgene', paste0(name, '_marker_normal_2'), 'pdf'), width = 16, height = 16);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes[31:60], pt.size = 4, nCol = 5, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
pdf(generate.filename('plotgene', paste0(name, '_marker_normal_3'), 'pdf'), width = 16, height = 16);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes[31:60], pt.size = 4, nCol = 5, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
pdf(generate.filename('plotgene', paste0(name, '_marker_normal_4'), 'pdf'), width = 16, height = 16);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes[61:90], pt.size = 4, nCol = 5, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
pdf(generate.filename('plotgene', paste0(name, '_marker_normal_5'), 'pdf'), width = 16, height = 16);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes[91:104], pt.size = 4, nCol = 5, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
####
genes.marker <- c('PTGS2', 'TFRC', 'SLC3A2', 'FSP1', 'GPX4', 'SLC7A1', 'IL1A', 'IL1B');
genes.marker <- intersect(genes.marker, rownames(seurat.all@data));
pdf(generate.filename('plotgene', paste0(name, '_iron'), 'pdf'), width = 15, height = 9);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 2, nCol = 4, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = FALSE, no.axes = TRUE);
dev.off();
###
genes.marker <- c('GPX4', 'SLC7A11', 'SLC3A2', 'FSP1', 'ACSL4', 'TFRC');
genes.marker <- intersect(genes.marker, rownames(seurat.all@data));
pdf(generate.filename('plotgene', paste0(name, '_iron2'), 'pdf'), width = 15, height = 9);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 2, nCol = 4, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = FALSE, no.axes = TRUE);
dev.off();
