library(BoutrosLab.plotting.general);
library(Seurat);
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-11_seurat_reduce_all.rds');
all_dr <- data.frame(cbind(seurat.all@dr$tsne@cell.embeddings, cluster = seurat.all@meta.data$cluster));
saveRDS(all_dr, file = generate.filename('embedding_cluster', 'all', 'rds'));
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
name <- 'tcell'
pdf(generate.filename('plottsne', paste0(name, '_cluster'), 'pdf'), width = 6, height = 6)
DimPlot(seurat.t, reduction.use = 'tsne', pt.size = 3, vector.friendly = FALSE);
dev.off();
t_dr <- data.frame(cbind(seurat.t@dr$tsne@cell.embeddings, cluster = seurat.t@meta.data$cluster));
saveRDS(t_dr, file = generate.filename('embedding_cluster', 't', 'rds'));
genes.marker <- 'KLK3'
pdf(generate.filename('plotgene', name, 'pdf'), width = 9, height = 9);
FeaturePlot(seurat.t, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 1, nCol = 1, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = FALSE, no.axes = TRUE);
dev.off();
pdf(generate.filename('plotgene', 'all', 'pdf'), width = 9, height = 9);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 1, nCol = 1, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = FALSE, no.axes = TRUE);
dev.off();
####
seurat.e <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/endo/objects/Manual_seurat_endo.rds');
name <- 'endo';
pdf(generate.filename('plottsne', paste0(name, '_cluster'), 'pdf'), width = 6, height = 6)
DimPlot(seurat.e, reduction.use = 'tsne', pt.size = 3, vector.friendly = FALSE);
dev.off();
endo_dr <- data.frame(cbind(seurat.e@dr$tsne@cell.embeddings, cluster = seurat.e@meta.data$res.0.8));
saveRDS(endo_dr, file = generate.filename('embedding_cluster', 'endo', 'rds'));
