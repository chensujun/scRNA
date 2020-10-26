library(BoutrosLab.plotting.general);
library(Seurat);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/seurat/");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
name <- 'all'
mycount <- seurat.all@data;
##### compare pam50 signature
mypam <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/pam50/2019-04-05_pam_score_c17.rds');
seurat.all@meta.data[, (ncol(seurat.all@meta.data) + 1):(ncol(seurat.all@meta.data)+3)] <- data.frame(mypam[rownames(seurat.all@meta.data), ]);
mywidth = 8;
myheight = 4;
iseurat <- SubsetData(seurat.all, cells.use = rownames(seurat.all@meta.data[!seurat.all@meta.data$type%in%c('Unknown'), ]))
fontsize2 <- theme(axis.text=element_text(size=12), axis.title.y=element_text(size=12, colour = 'black'), axis.title.x = element_text(size = 0),
plot.title = element_text(size = 12, face = 'plain'),
plot.margin=unit(c(.5,.5,1,2.5),"cm"), axis.text.x = element_text(angle = 45, hjust = 1, size = 12));
fontsize1 <- theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, colour = 'black'), axis.title.x = element_text(size = 0),
    plot.title = element_text(size = 12, face = 'plain'), 
    plot.margin=unit(c(.5,.5,1,2.5),"cm"), axis.text.x = element_text(size = 0));

pdf(generate.filename('plotviolin', paste0(name, '_pam_luminalA'), 'pdf'), width = mywidth, height = myheight-1);
VlnPlot(iseurat, features.plot = 'lum.a.score', x.lab.rot = TRUE, point.size.use = 0) + ggtitle('Luminal A') +
        labs(x = '', y = expression(paste('Signature score'))) + fontsize1 + 
        stat_summary(fun.y = median, geom = 'point', size = 12, colour = 'black', shape = 95) + 
        geom_hline(yintercept = 0, linetype = 'dashed', size = 1);
dev.off();

pdf(generate.filename('plotviolin', paste0(name, '_pam_luminalB'), 'pdf'), width = mywidth, height = myheight-1);
VlnPlot(iseurat, features.plot = 'lum.b.score', x.lab.rot = TRUE, point.size.use = 0) + ggtitle('Luminal B') +
        labs(x = '', y = expression(paste('Signature score'))) + fontsize1 + 
        stat_summary(fun.y = median, geom = 'point', size = 12, colour = 'black', shape = 95) + 
        geom_hline(yintercept = 0, linetype = 'dashed', size = 1);
dev.off();

pdf(generate.filename('plotviolin', paste0(name, '_pam_basal'), 'pdf'), width = mywidth, height = myheight);
VlnPlot(iseurat, features.plot = 'basal.score', x.lab.rot = TRUE, point.size.use = 0) + ggtitle('Basal') +
        labs(x = '', y = expression(paste('Signature score'))) + fontsize2 + 
        stat_summary(fun.y = median, geom = 'point', size = 12, colour = 'black', shape = 95) + 
        geom_hline(yintercept = 0, linetype = 'dashed', size = 1);
dev.off();
####
iseurat <- readRDS(conf$seurat_endo);
iseurat@meta.data <- seurat.all@meta.data[rownames(iseurat@meta.data), ];
iseurat@meta.data$cluster_new <- iseurat@ident;
####
seurat.all@meta.data$cluster_new <- as.vector(iseurat@meta.data[match(rownames(seurat.all@meta.data), rownames(iseurat@meta.data)), ]$cluster_new);
seurat.all@meta.data$cluster_new[is.na(seurat.all@meta.data$cluster_new)] <- '';
jseurat <- SubsetData(seurat.all, cells.use = rownames(seurat.all@meta.data[seurat.all@meta.data$type%in%c('Endothelial', 'Basal/intermediate'), ]));
jseurat@meta.data$group <- paste0(jseurat@meta.data$type, jseurat@meta.data$cluster_new);
name <- 'basal_endo';
pdf(generate.filename('plotviolin', paste0(name, '_pam_basal'), 'pdf'), width = mywidth, height = myheight);
VlnPlot(jseurat, features.plot = 'basal.score', x.lab.rot = TRUE, point.size.use = 0, group.by = 'group') + ggtitle('Basal') +
        labs(x = '', y = expression(paste('Signature score'))) + fontsize2 + 
        stat_summary(fun.y = median, geom = 'point', size = 12, colour = 'black', shape = 95) + 
        geom_hline(yintercept = 0, linetype = 'dashed', size = 1);
dev.off();
####
