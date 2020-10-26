library(BoutrosLab.plotting.general);
library(Seurat);
library(matrixStats);
library(reshape);
library(plyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/GS')
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
# remove TMEM16G as it's low GS associated
# TEGT to TMBIM6
gs <- c("ANXA5", "APLP2", "COL1A2", "CRIM1", "ECHS1", "HMGB1", "MAP3K5", "NGFRAP1", "NPM1", "NUB1", "PDIA3", "PLA2G2A", "ROCK1", "TMBIM6", "TMEM69", "TRAF4", "VCAN", "VCP", "VDAC1");
seurat.all@meta.data$gs <- colMeans(seurat.all@data[rownames(seurat.all@data)%in%gs, ]);
myexp <- data.frame(t(scale(t(seurat.all@data[gs, ]))));
seurat.all@meta.data$gs.s <- colMeans(myexp)
pdf(generate.filename('plotvln', paste0('gs_sig'), 'pdf'), width = 5, height = 5);
VlnPlot(seurat.all, features.plot = 'gs', x.lab.rot = TRUE, point.size.use = 0, group.by = 'type') + 
	stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95) + 
	geom_hline(yintercept = 0.63, linetype = 'dashed', size = 1);
dev.off();
####
pdf(generate.filename('plotvln', paste0('gs_sig_scale'), 'pdf'), width = 5, height = 5);
VlnPlot(seurat.all, features.plot = 'gs.s', x.lab.rot = TRUE, point.size.use = 0, group.by = 'type') + 
	stat_summary(fun.y = median, geom = 'point', size = 10, colour = 'black', shape = 95) + 
	geom_hline(yintercept = median(seurat.all@meta.data[seurat.all@meta.data$type=='Luminal', ]$gs.s), linetype = 'dashed', size = 1);
dev.off();
####
myexp <- data.frame(t(seurat.all@data[gs, ]));
myexp$type <- seurat.all@meta.data$type;
myexp <- ddply(myexp, 'type', numcolwise('mean'));
rownames(myexp) <- myexp$type;
myexp <- myexp[, -1];
myfc <- sapply(seq(ncol(myexp)), function(y) max(sapply(c(seq(1,3), seq(5, 7)), function(x) myexp[x, y]/myexp[4, y])));
names(myfc) <- colnames(myexp);
#genes.marker <- colnames(myexp[, myexp['Luminal', ]<0.1]);
genes.marker <- names(myfc[myfc>2]);
pdf(generate.filename('plotvln', paste0('gs_genes'), 'pdf'), width = 5, height = 15);
VlnPlot(seurat.all, features.plot = genes.marker, x.lab.rot = TRUE, point.size.use = 0, group.by = 'type', nCol = 1) + 
	stat_summary(fun.y = median, geom = 'point', size = 10, colour = 'black', shape = 95);
dev.off();
###
myexp <- data.frame(t(seurat.all@data[gs, ]));
myexp$type <- seurat.all@meta.data$type;
myexp <- ddply(myexp, 'type', numcolwise('median'));
rownames(myexp) <- myexp$type;
myexp <- myexp[, -1];
myexp <- myexp + 0.01
myfc <- sapply(seq(ncol(myexp)), function(y) max(sapply(c(seq(1,3), seq(5, 7)), function(x) myexp[x, y]/myexp[4, y])));
names(myfc) <- colnames(myexp);
#genes.marker <- colnames(myexp[, myexp['Luminal', ]<0.1]);
genes.marker <- names(myfc[myfc>2]);
