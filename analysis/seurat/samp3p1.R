library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(readxl);
library(plyr);
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/samp3p1');
source('/u/schen/svn/singleCell/myfunctions/grep_term.R');
source('~/svn/singleCell/myfunctions/plot_enrich_list.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
iseurat <- readRDS(conf$samp3p1);
name <- 'samp3p1';
mytype <- read.table('objects/CellType_Rename.txt', header = TRUE);
iseurat@meta.data$type <- mytype[match(rownames(iseurat@meta.data), mytype$Cell), ]$Cluster;
pdf(generate.filename('plotgene_samp', name, 'pdf'), width = 7, height = 6);
DimPlot(iseurat, reduction.use = 'tsne', group.by = 'orig.ident', pt.size = 2, no.axes = TRUE, vector.friendly = TRUE, cols.use = c('darkblue', 'orange', 'green3'));
dev.off();

pdf(generate.filename('plotgene_type', name, 'pdf'), width = 6, height = 6);
DimPlot(iseurat, reduction.use = 'tsne', group.by = 'type', pt.size = 2, no.axes = TRUE, vector.friendly = TRUE, do.label = TRUE, no.legend = TRUE);
dev.off();
####
pdf(generate.filename('plotviolin', paste0(name, '_klk3'), 'pdf'), width = 10);
VlnPlot(iseurat, features.plot = 'KLK3', x.lab.rot = TRUE, point.size.use = 0, group.by = 'samp_type')  + ggtitle('') +
	labs(x = '', y = expression('log(nUMI)')) +
	stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95) + 
	geom_hline(yintercept = 0.05, linetype = 'dashed', size = 2);
dev.off();

pdf(generate.filename('plotviolin', paste0('all', '_klk3'), 'pdf'), width = 10);
VlnPlot(seurat.all, features.plot = 'KLK3', x.lab.rot = TRUE, point.size.use = 0, group.by = 'type')  + ggtitle('') +
	labs(x = '', y = expression('log(nUMI)')) +
	stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95) + 
	geom_hline(yintercept = 0.05, linetype = 'dashed', size = 2);
dev.off();

pdf(generate.filename('dotmap', 'markers', 'pdf'), width = 4, height = 12);
DotPlot(iseurat, genes.plot = 'KLK3', cols = c('blue', 'red'), dot.scale = 6, group.by = 'samp_type',
	x.lab.rot = TRUE, plot.legend = TRUE);
dev.off();

