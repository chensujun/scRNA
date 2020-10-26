library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
name <- 'tcell';
doublet <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/data/Doublet_Result/AllSample.doublet.txt', header = TRUE);
seurat.t@meta.data$doublet <- doublet[match(rownames(seurat.t@meta.data), doublet$Cell), ]$isDoublet;

mywidth <- 12;
myheight <- 9;
myres <- 300;
#### SF2C
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));

tiff(generate.filename('plottsne', paste0(name, '_doublet'), 'tiff'), width = mywidth, height = myheight, res = myres, units = 'in')
DimPlot(seurat.t, reduction.use = 'tsne', pt.size = 2, group.by = 'doublet') + fontsize + labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
###
myexp <- data.frame(exp = seurat.t@data['KLK3', ]);
myexp$doublet <- seurat.t@meta.data$doublet;
create.violinplot(
    formula = exp~doublet,
    data = myexp,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance'),
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_doublet', 'pdf'),
    #add.text = TRUE,
    #text.x = 1.5,
    #text.y = 7,
    #text.label = pval
    );
###
myexp <- data.frame(exp = log(seurat.t@raw.data['KLK3', colnames(seurat.t@data)]+1));
myexp$doublet <- seurat.t@meta.data$doublet;
create.violinplot(
    formula = exp~doublet,
    data = myexp,
    xlab.label = 'Group',
    ylab.label = 'log(UMI+1)',
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_doublet', 'pdf'),
    #add.text = TRUE,
    #text.x = 1.5,
    #text.y = 7,
    #text.label = pval
    );
