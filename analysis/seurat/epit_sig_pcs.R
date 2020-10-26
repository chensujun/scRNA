library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(Seurat);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/epitumour");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
name <- 'epi';
iseurat <- readRDS(conf$seurat_epi13);
iseurat <- SetIdent(iseurat, ident.use = iseurat@meta.data$fig.patient);
pcs1 <- c('STMN1', 'MCM4', 'CCNB1', 'CDC6', 'CDKN3', 'EZH2', 'TPX2', 'FOXM1', 'HMMR', 'MKI67', 'KNTC1');
iseurat@meta.data$pcs1 <- colMeans(iseurat@data[pcs1, ]);
fontsize2 <- theme(axis.text=element_text(size=12), axis.title.y=element_text(size=12, colour = 'black'), axis.title.x = element_text(size = 0),
plot.title = element_text(size = 12, face = 'plain'),plot.margin=unit(c(.5,.5,1,2.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1, size = 12));

mywidth = 10
myheight = 4
pdf(generate.filename('plotviolin', paste0(name, '_sig_pcs1'), 'pdf'), width = mywidth, height = myheight);
VlnPlot(iseurat, features.plot = 'pcs1', x.lab.rot = TRUE, point.size.use = 0, group.by = 'fig.cluster') +
        labs(x = 'Cluster', y = expression(paste('Signature score'))) + fontsize2 +
        stat_summary(fun.y = median, geom = 'point', size = 12, colour = 'black', shape = 95) +
        geom_hline(yintercept = 0, linetype = 'dashed', size = 1);
dev.off();

