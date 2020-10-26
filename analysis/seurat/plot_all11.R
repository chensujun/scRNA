library(BoutrosLab.plotting.general);
library(Seurat);
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$seurat_all11);
name <- 'all11';
pdf(generate.filename('plottsne', paste0(name, '_type'), 'pdf'), width = 8, height = 6)
DimPlot(seurat.all, reduction.use = 'TSNE', group.by = 'type', pt.size = 1, vector.friendly = TRUE);
dev.off();
genes.major <- c('CD3D', 'CD79A', 'C1QA', 'CD163', 'EPCAM', 'PECAM1', 'ACTA2', 'DCN', 'S100A9');
pdf(generate.filename('marker_major', name, 'pdf'));
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.major, pt.size = 2,
        cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE)
dev.off();
