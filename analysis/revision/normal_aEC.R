library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/EC_normal');
seurat.e <- readRDS('GraphClust.seuset-cellreport-prostate-EC.rds');
#seurat.f <- readRDS('GraphClust.seuset-cellreport-prostate-Fibroblast.rds');
mywidth <- 12;
myheight <- 12;
myres <- 300;
name <- 'endo';
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));
fontsize1 <- theme(axis.text=element_text(size=24), axis.title=element_text(size=24), legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1));

genes.marker <- c('PECAM1', 'THY1');
pdf(generate.filename('plotgene', name, 'pdf'), width = 15, height = 6);
FeaturePlot(seurat.e, reduction.use = 'umap', features.plot = genes.marker, pt.size = 4, nCol = 2, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
