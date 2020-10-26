library(BoutrosLab.plotting.general);
library(Seurat);
library(viridis);
library(reshape2);
library(readxl);
library(canprot);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/F1');
source('~/chensj/scRNA/script/myfunctions/get_dotplot_data.R');
source('~/chensj/scRNA/script/myfunctions/DotPlot_flip.R');
source('~/chensj/scRNA/script/myfunctions/utilities_internal.R');
source('~/chensj/scRNA/script/myfunctions/plot_distribution_new.R');
p_theme<- theme(
  text=element_text(family="Arial"),
  axis.title.x =element_text(color="black",size=18,family="Arial") ,
  axis.title.y =element_text(color="black",size=18,family="Arial") ,
  axis.text.x =element_text(color="black",size=16,family="Arial") ,
  axis.text.y =element_text(color="black",size=16,family="Arial") ,
  legend.text =element_text(color="black",size=16,family="Arial"),
  legend.title=element_text(color="black",size=18,family="Arial")
)
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
name <- 'all';
seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monocytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
seurat.all <- SetAllIdent(seurat.all, id = 'type');
####
mywidth <- 16;
myheight <- 12;
myres <- 300;
#### F1A
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));
tiff(generate.filename('plottsne', paste0(name, '_type'), 'tiff'), width = mywidth, height = mywidth, res = myres, units = 'in')
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'type', do.label = TRUE, label.size = 19, no.legend = TRUE, no.axes = FALSE) + 
fontsize + labs(x = 'Dimension 1', y = 'Dimension 2')
dev.off();

tiff(generate.filename('plottsne', paste0(name, '_type_nolabel'), 'tiff'), width = mywidth, height = mywidth, res = myres, units = 'in')
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'type', do.label = FALSE, label.size = 19, no.legend = TRUE, no.axes = FALSE) + 
fontsize + labs(x = 'Dimension 1', y = 'Dimension 2')
dev.off();

pdf(generate.filename('plottsne', paste0(name, '_type'), 'pdf'), width = mywidth, height = myheight);
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'type', do.label = TRUE, vector.friendly = TRUE, label.size = 24, no.legend = TRUE, no.axes = FALSE) + 
fontsize + labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
####F1B
markers <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/no_hk_filter/annotate_type/2019-03-01_annotate_type_markers_27.rds');
cellTypeMarks_csv= read.csv('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/collectedBiomarkers_3_celltype.csv', stringsAsFactors=FALSE);
temp_geneSymbols<- strsplit(cellTypeMarks_csv$geneSymbol, ", ")
# create a list object to store the markers
markers=list()
for(i in 1:nrow(cellTypeMarks_csv))
{
  temp_celltype=gsub(' |\\/', '_', cellTypeMarks_csv$cellName[i])
  temp_markers=temp_geneSymbols[[i]]
  if(temp_celltype %in% names(markers))
  {
    markers[[temp_celltype]]=union(temp_markers,markers[[temp_celltype]])
  }else 
  {
    markers[[temp_celltype]]=c(temp_markers)
  }
  markers[[temp_celltype]] <- na.omit(markers[[temp_celltype]])
  markers[[temp_celltype]] <- sapply(markers[[temp_celltype]], function(x) gsub('\\?', '', x))
};
markers[['basal_intermediate']] <- c('TP63', 'KRT14', 'KRT5');

mymarker <- vector();
for(i in c('T_cell', 'Macrophage', 'basal_intermediate', 'Luminal_cell', 'Mast_cells_cell201801', 'Endothelial_cells_cell201801', 'Myofibroblast')){
	mymarker <- c(mymarker, markers[[i]])
};
mymarker <- mymarker[!duplicated(mymarker)]
mymarker <- mymarker[mymarker%in%rownames(seurat.all@data)];

iseurat <- seurat.all;
iseurat@ident <- factor(iseurat@ident, levels = (c('T', 'Monolytic', 'Basal/intermediate', 'Luminal', 'Mast', 'Endothelial', 'Fibroblast')));
mymarker <- rev(mymarker);

tiff(generate.filename('dotmap', 'markers', 'tiff'), width = 5, height = 10, res = 300, units = 'in')
DotPlot_flip(iseurat, genes.plot = mymarker, cols = c('white', 'blue'), dot.scale = 6, 
	x.lab.rot = TRUE, plot.legend = TRUE) + 
theme(axis.text=element_text(size=24, colour = 'black'), legend.justification = 'top') + p_theme;
dev.off();

pdf(generate.filename('dotmap', 'markers', 'pdf'), width = 6, height = 10)
DotPlot_flip(iseurat, genes.plot = mymarker, cols = c('white', 'blue'), dot.scale = 6, 
	x.lab.rot = TRUE, plot.legend = TRUE) + 
theme(axis.text=element_text(size=24, colour = 'black'), legend.justification = 'top') + p_theme;
dev.off();
#####
#### F1C
pdf(generate.filename('plottsne', paste0(name, '_nGene'), 'pdf'), width = 4, height = 3);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = c('nGene'), pt.size = 0.5, vector.friendly = TRUE, 
  nCol = , no.axes = FALSE, no.legend = FALSE);
dev.off();
#####
name <- 'all';
to.plot <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/multi//2020-05-26_cnv_score_pri_0518.rds');
rownames(to.plot) <- gsub('\\.', '-', rownames(to.plot));
to.plot$type <- seurat.all@meta.data[rownames(to.plot), ]$type;
seurat.all@meta.data$cnv <- to.plot[rownames(seurat.all@meta.data), ]$score;
fontsize2 <- theme(axis.text=element_text(size=48), axis.title.y=element_text(size=48, colour = 'black'), axis.title.x = element_text(size = 0),
  plot.title=element_text(size = 0),
    plot.margin=unit(c(.5,.5,1,2.5),"cm"), axis.text.x = element_text(angle = 45, hjust = 1, size = 48));

wilcox.test(to.plot[to.plot$type%in%c('Luminal', 'Basal/intermediate'), ]$score, to.plot[!to.plot$type%in%c('Luminal', 'Basal/intermediate'), ]$score)$p.value;
CLES(to.plot[to.plot$type%in%c('Luminal', 'Basal/intermediate'), ]$score, to.plot[!to.plot$type%in%c('Luminal', 'Basal/intermediate'), ]$score);

seurat.all@meta.data$type <- gsub('Basal/intermediate', 'Basal/int.', seurat.all@meta.data$type);
tiff(generate.filename('plotviolin', paste0(name, '_cnv'), 'tiff'), width = mywidth , height = myheight - 4, res = myres, units = 'in')
VlnPlot(seurat.all, features.plot = 'cnv', x.lab.rot = TRUE, point.size.use = 0, group.by = 'type') + fontsize2 + ggtitle('') +
  #labs(x = '', y = expression(paste('(CNV-1)'^'2'))) +
  labs(x = '', y = expression('Ave. deviation')) +
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95) + 
  geom_hline(yintercept = 0.05, linetype = 'dashed', size = 2);
dev.off();


pdf(generate.filename('plottsne', paste0(name, '_cnv'), 'pdf'), width = myheight, height = myheight);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = 'cnv', pt.size = 1, nCol = 1, min.cutoff = 0.01,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);

dev.off();


seurat.all@meta.data$cnv[seurat.all@meta.data$cnv>0.2] <- 0.2;
tiff(generate.filename('plotviolin', paste0(name, '_cnv0.2'), 'tiff'), width = mywidth , height = myheight - 4, res = myres, units = 'in')
VlnPlot(seurat.all, features.plot = 'cnv', x.lab.rot = TRUE, point.size.use = 0, group.by = 'type') + fontsize2 + ggtitle('') +
  #labs(x = '', y = expression(paste('(CNV-1)'^'2'))) +
  labs(x = '', y = expression('Ave. deviation')) +
  stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95) + 
  geom_hline(yintercept = 0.05, linetype = 'dashed', size = 2);
dev.off();


#####
iseurat <- seurat.all;
annot <- data.frame(read_excel('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F1/Table_S1_SampleInfo.xlsx'));
annot$pGS <- annot$pGleason_score;
annot$IDCP <- annot$IDC.CA*100;
annot$pre.treatment.psa <- annot$PSA;
annot$pathological_t <- annot$pT;
annot$sample <- paste0(annot$sample, 'SL');
annot$idc <- ifelse(annot$IDCP>0, TRUE, FALSE);
percent.tab <- as.data.frame.matrix(table(iseurat@meta.data$orig.ident, iseurat@meta.data$type));
mycol <- readRDS('~/chensj/scRNA/primary/scran/data/2019-07-25_mycol_type.rds');
plot_distribution(percent.tab, 'type_all', annot, ncell.max = max(percent.tab), mytype = 'tiff', width = 10, height = 2.8, mycol = mycol);
names(mycol) <- gsub('Myeloid', 'Monocytic', names(mycol))
plot_distribution(percent.tab, 'type_all', annot, ncell.max = max(percent.tab), mytype = 'tiff', width = 12.5, height = 3.5, mycol = mycol);
