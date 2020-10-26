library(BoutrosLab.plotting.general);
#library(clusterProfiler);
library(pheatmap);
library(Seurat);
library(plyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
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
seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monocytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
seurat.all <- SetAllIdent(seurat.all, id = 'type');
annot.go <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');

#seurat.all <- ScaleData(seurat.all);
#saveRDS(seurat.all@scale.data, file= generate.filename('seurat_scale', 'all', 'rds'));
mygenes <- c(intersect(annot.go[annot.go$name=='T cell costimulation', ]$gene, rownames(seurat.all@data)), 'KIRREL')
pdf(generate.filename('heatmap', 'Tcostimulation', 'pdf'), width = 10);
print(DoHeatmap(seurat.all, genes.use = mygenes, 
        slim.col.label = TRUE, remove.key = FALSE));
dev.off();

pdf(generate.filename('plotgene', paste0('all', '_tcellCostimul'), 'pdf'), width = 16, height = 16);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = mygenes, pt.size = 2, nCol = 9, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE
                );
dev.off();
####
to.plot <- seurat.all@data[mygenes, ];
to.plot <- to.plot[rowMeans(to.plot)>log(0.1+1), ];
####
group <- ifelse(colnames(to.plot)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]), 1, 2);
pval <- data.frame(pval = apply(to.plot, 1, function(x) wilcox.test(x~group)$p.value));
pval$lfc <- apply(to.plot+0.1, 1, function(x) log2(mean(x[group==1])/mean(x[group==2])));

to.plot <- data.frame(t(scale(t(to.plot))));
ann_col <- data.frame(type = seurat.all@meta.data$type);
rownames(ann_col) <- gsub('-', '.', rownames(seurat.all@meta.data));
ann_col <- ann_col[order(ann_col$type), ,drop = FALSE];
mycol.type <- readRDS('~/chensj/scRNA/primary/scran/data/2019-07-25_mycol_type.rds');
names(mycol.type) <- gsub('Myeloid', 'Monocytic', names(mycol.type))
to.plot <- to.plot[, match(rownames(ann_col), colnames(to.plot))];
mybreaks <- unique(c(seq(-13, -1.5, 10), seq(-1.5, 1.5, 0.01), seq(1.5, 13, 10)));
mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
ann_row <- data.frame(meanUMI = rowMeans(exp(seurat.all@data[rownames(to.plot), ]-1)));
rownames(ann_row) <- rownames(to.plot)
ann_row$meanUMI <- exp(ann_row$meanUMI)-1
ann_colors <- list(
  type = mycol.type
  )
pdf(generate.filename('pheatmap', 'Tcostimulation', 'pdf'), width = 10, height = 5);
pheatmap((to.plot), cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
	breaks = mybreaks, color = mycol, annotation_col = ann_col, annotation_colors = ann_colors);
dev.off();

png(generate.filename('pheatmap', 'Tcostimulation', 'png'), width = 8, height = 4, unit= 'in', res = 300);
pheatmap((to.plot), cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
  breaks = mybreaks, color = mycol, annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_colors);
dev.off();

#set.seed(100);
set.seed(1);
hc <- hclust(dist(t(to.plot)));
save(hc, file = generate.filename('Tcostimulation_clusterCol', 'hc', 'rds'));
to.plot <- to.plot[, hc$order];
ann_col <- ann_col[colnames(to.plot), ,drop =FALSE];
png(generate.filename('pheatmap', 'Tcostimulation_clusterCol', 'png'), width = 8, height = 4, unit= 'in', res = 300);
pheatmap((to.plot), cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
  breaks = mybreaks, color = mycol, annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_colors);
dev.off();

####
to.plot <- data.frame(t(seurat.all@data[mygenes, ]));
to.plot$type <- seurat.all@meta.data$type;
to.plot <- ddply(to.plot, 'type', numcolwise(mean));
rownames(to.plot) <- to.plot$type;
to.plot <- data.frame(t(to.plot[, -1]));

mymarkers <- rownames(to.plot[to.plot$Luminal>=to.plot['DPP4', ]$Luminal, ]);
pdf(generate.filename('plotviolin', 'higher_tcellCostimul', 'pdf'), width = 6, height = 4.5);
VlnPlot(seurat.all, features.plot = mymarkers, point.size.use = 0, nCol = 3);
dev.off();
