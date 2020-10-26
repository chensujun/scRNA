library(BoutrosLab.plotting.general);
library(clusterProfiler);
library(Seurat);
library(qusage);
library(viridis);
library(plyr);
library(tidyr);
library(dplyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/SF1/');
source('~/chensj/scRNA/script/myfunctions/run_qusage_seurat.R');
source('~/chensj/scRNA/script/myfunctions/get_dotplot_data.R');
source('~/chensj/scRNA/script/myfunctions/DotPlot_flip.R');
source('~/chensj/scRNA/script/myfunctions/utilities_internal.R');
p_theme<- theme(
  text=element_text(family="Arial"),
  axis.title.x =element_text(color="black",size=18,family="Arial") ,
  axis.title.y =element_text(color="black",size=18,family="Arial") ,
  axis.text.x =element_text(color="black",size=16,family="Arial") ,
  axis.text.y =element_text(color="black",size=16,family="Arial") ,
  legend.text =element_text(color="black",size=16,family="Arial"),
  legend.title=element_text(color="black",size=18,family="Arial")
);
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monolytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
seurat.all <- SetAllIdent(seurat.all, id = 'type');
####
seurat.f <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/fibroblast/objects/13_sample_fibroblast_20190718/GraphClust.seuset.rds');
annot.keg <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-05_kegg_human_20190613.rds');
annot.go <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
####
seurat.all@meta.data$cluster <- as.numeric(factor(seurat.all@meta.data$type));
gs <- list(emt = annot.go[annot.go$name=='epithelial to mesenchymal transition', ]$gene);
myresults <- run_qusage(seurat.all, 'all_emt', gs);
saveRDS(myresults, generate.filename('qusage', 'all_emt', 'rds'));
####
h1 <- clusterProfiler::read.gmt('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/data/h.all.v7.0.symbols.gmt');
gs[['h1_emt']] <- h1[h1$ont=='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', ]$gene
seurat.all@meta.data$emt <- colMeans(seurat.all@data[rownames(seurat.all@data)%in%gs[[1]], ]);
seurat.all@meta.data$emt.h1 <- colMeans(seurat.all@data[rownames(seurat.all@data)%in%gs[[2]], ]);
seurat.all@meta.data$acta2 <- seurat.all@data['ACTA2', ];


mymarker <- c('S100A4', 'SPARC', 'ACTA2', 'PDGFRB', 'CAV1', 'VIM');
tiff(generate.filename('dotmap', 'markers', 'tiff'), width = 5, height = 5, res = 300, units = 'in')
DotPlot_flip(seurat.all, genes.plot = mymarker, cols = c('white', 'blue'), dot.scale = 6, 
	x.lab.rot = TRUE, plot.legend = TRUE) + 
theme(axis.text=element_text(size=24, colour = 'black'), legend.justification = 'top') + p_theme;
dev.off();

to.plot <- seurat.all@meta.data;
to.plot <- to.plot[to.plot$type%in%c('Luminal', 'Basal/intermediate'), ]
to.plot <- ddply(to.plot, 'orig.ident', numcolwise(mean));
to.plot$nCAF <- table(seurat.all@meta.data$orig.ident, seurat.all@meta.data$type)[, 'Fibroblast'];
to.plot$nCell <- table(seurat.all@meta.data$orig.ident);
to.plot$pCAF <- 100*to.plot$nCAF/to.plot$nCell;
to.plot$nCAF.2 <- table(seurat.all@meta.data[seurat.all@meta.data$acta2>0&seurat.all@meta.data$type=='Fibroblast', ]$orig.ident)[to.plot$orig.ident]
to.plot[is.na(to.plot)] <- 0;
to.plot$pCAF.2 <- 100*to.plot$nCAF.2/to.plot$nCell;
####
create.scatterplot(
	formula = pCAF.2~emt.h1,
	data = to.plot,
	xlab.label = 'EMT score in epithelial cells',
	ylab.label = 'ACTA2+ CAF percentage',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	#xlimits = c(0, 2.5),
	filename = generate.filename('fib', 'emt', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$pCAF.2,
                y = to.plot$emt.h1,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);

create.scatterplot(
	formula = pCAF.2~emt.h1,
	data = to.plot,
	xlab.label = 'EMT score in epithelial cells',
	ylab.label = 'ACTA2+ CAF percentage',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	#xlimits = c(0, 2.5),
	filename = generate.filename('fib', 'emt', 'tiff'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$pCAF.2,
                y = to.plot$emt.h1,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);
###
write.csv(to.plot, generate.filename('plotdata', 'F_S1_reviewer', 'csv'));