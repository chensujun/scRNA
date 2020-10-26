library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/SF1');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
name <- 'all';
mywidth <- 12;
myheight <- 9;
myres <- 300;
#### SF2C
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));
tiff(generate.filename('plottsne', paste0(name, '_type2'), 'tiff'), width = mywidth, height = myheight, res = myres, units = 'in')
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = .5, group.by = 'type_scnorm') + fontsize + labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
