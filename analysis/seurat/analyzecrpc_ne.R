library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(reshape2);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm");
source('~/svn/singleCell/myfunctions/plot_clust.R');
source('~/svn/singleCell/myfunctions/enrich_gprofiler.R');
source('~/svn/singleCell/myfunctions/go_simplify.R');
source('~/svn/singleCell/myfunctions/plot_enrich.R');
source('~/svn/singleCell/myfunctions/plot_distribution.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.cr <- readRDS(conf$seurat_crpc);
name <- 'crpc'
mycol <- readRDS(conf$col_type);
###

mycount <- seurat.cr@data;
genes.marker <- readRDS(conf$Epi_subtype);
genes.ne <- colnames(genes.marker)[1:7]
iseurat <- SubsetData(seurat.cr, cells.use = rownames(seurat.cr@meta.data[seurat.cr@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]));
genes.marker <- data.frame(t(mycount[c(genes.ne, "KLK3", "AR", "ACPP", "KRT8"), rownames(iseurat@meta.data)]));
genes.marker$patient <- iseurat@meta.data[rownames(genes.marker), ]$orig.ident;
genes.marker$group <- 'Others';
genes.marker[rowSums(genes.marker[, genes.ne])>0&genes.marker$AR==0, ]$group <- 'NE';
genes.marker[rowSums(genes.marker[, genes.ne])==0&genes.marker$AR==0, ]$group <- 'DN';
genes.marker[rowSums(genes.marker[, genes.ne])==0&genes.marker$KLK3>3.7&genes.marker$AR>0, ]$group <- 'Lum';

mywidth <- 7;
myheight <- 6;
myres <- 300;

pdf(generate.filename('plottsne', paste0(name, '_type'), 'pdf'), width = mywidth, height = myheight)
DimPlot(seurat.cr, reduction.use = 'TSNE', pt.size = 1, group.by = 'type', no.axes = TRUE, 
	cols.use = mycol[levels(seurat.cr@ident)]) 
dev.off();

pdf(generate.filename('plottsne', paste0(name, '_sample'), 'pdf'), width = mywidth, height = myheight)
DimPlot(seurat.cr, reduction.use = 'TSNE', pt.size = 1, group.by = 'orig.ident', no.axes = TRUE) 
dev.off();
#####
annot <- readRDS(conf$annot);
iseurat <- seurat.cr
percent.tab <- as.data.frame.matrix(table(iseurat@meta.data$orig.ident, iseurat@meta.data$type));
to.plot <- apply(percent.tab, 1, function(x) 100*x/sum(x));
colnames(to.plot) <- gsub('JD1800|SL', '', colnames(to.plot));
to.plot <- melt(to.plot);
#to.plot$Var1 <- gsub('Normal epithelial', 'Epithelial', to.plot$Var1);
to.plot$Var2 <- factor(to.plot$Var2);
to.plot$Var1 <- factor(to.plot$Var1);
bar.plot <- create.barplot(
	file = generate.filename('distribution', name, 'pdf'),
	formula = Var2~value,
	data = to.plot,
	groups = Var1,
	col = mycol[levels(to.plot$Var1)],
	stack = TRUE,
	plot.horizontal = TRUE,
	ylab.label = 'Patient',
	xlab.label = 'Percentage',
	xlimits = c(0, 100),
	xat = seq(0, 100, 20),
	yaxis.rot = 90,
	style = 'Nature',
	width = 6,
	height = 4,
	y.spacing = -1,
	border.col = 'white',
	box.ratio = 8
	);


#for(gene in c('C1QA', 'CD3D', 'S100A9', 'CCR7')){
#	for(gene in c('CD14', 'CD163', 'CSF1R', 'FCGR3A', 'LYZ', 'CSF1R', 'FCGR2A', 'CD163')){
#	pdf(generate.filename('marker', paste0(name, '_', gene, '_tsne'), 'pdf'), width = 5, height = 5);
#	FeaturePlot(seurat.cr, features.plot = gene, reduction.use = 'TSNE', pt.size = .3, cols.use = c('grey', 'red'));
#	dev.off();
#};

####
iseurat <- SubsetData(seurat.cr, cells.use = rownames(seurat.cr@meta.data[seurat.cr@meta.data$orig.ident=='QJZ', ]));
name <- 'QJZ';
pdf(generate.filename('plottsne', paste0(name, '_type'), 'pdf'), width = mywidth, height = myheight)
DimPlot(iseurat, reduction.use = 'TSNE', pt.size = 1, group.by = 'type', cols.use = mycol[levels(iseurat@ident)],
	no.axes = TRUE) 
dev.off();
####
iseurat <- SubsetData(seurat.cr, cells.use = rownames(seurat.cr@meta.data[seurat.cr@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]));
iseurat <- SetIdent(iseurat, ident.use = ifelse(iseurat@meta.data$orig.ident=='QJZ', 1, 2));
gene.go <- readRDS('~/circRNA/star-circ/circRNA_landscape/rnaseq_landscape/scRNA/2019-01-05_gp_hsapiens.GO.Name.rds');
find_diff <- function(iseurat, name, mycompare, gene.go){
	jm.out <- FindMarkers(iseurat, ident.1 = 1, ident.2 = 2, logfc.threshold = 0);
	iseurat <- ScaleData(iseurat);
	idiff <- jm.out;
	pdf(generate.filename('heatmap', paste0(name, '_diff50_NEvsCR'), 'pdf'), width = mywidth, height = myheight);
	print(DoHeatmap(iseurat, genes.use = c(rownames(idiff[order(idiff$avg_logFC), ])[1:25], rownames(idiff[order(-idiff$avg_logFC), ])[1:25]),
		slim.col.label = TRUE, remove.key = FALSE));
	dev.off();

	jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
	ego.up <- enrich_gprofiler(jup, rownames(seurat.cr@data));
	ego.dn <- enrich_gprofiler(jdn, rownames(seurat.cr@data));
	red.ego.up <- run_simplify(ego.up);
	red.ego.dn <- run_simplify(ego.dn);
	save(jm.out, ego.up, ego.dn, red.ego.up, red.ego.dn, file = generate.filename('diff', paste0(name, '_', mycompare), 'rda'));
	ego.all <- enrich_gprofiler(rownames(jm.out[jm.out$p_val_adj<0.05, ]), rownames(seurat.cr@data));
	red.ego.all <- run_simplify(ego.all);
	save(ego.all, red.ego.all, file = generate.filename('diff_together', paste0(name, '_', mycompare), 'rda'));
	myego <- list(up = red.ego.up, down = red.ego.dn);
	plot_enrich(myego, gene.go, 'BP', 'NEvsCR_alldiff', mybg = rownames(seurat.cr@data), width = 10, 
		spot.size.function = function(x) {abs(x)});
	return(jm.out);
};
jm.out <- find_diff(iseurat, 'QJZ', 'NEvsCR', gene.go);
#### compare NE to AR/KLK3+ cells
iseurat <- SubsetData(seurat.cr, cells.use = rownames(genes.marker[genes.marker$group%in%c('NE', 'Others')&genes.marker$patient=='QJZ', ]));
iseurat@meta.data$group <- genes.marker[rownames(iseurat@meta.data), ]$group;
iseurat <- SetIdent(iseurat, ident.use = ifelse(iseurat@meta.data$group=='NE', 1, 2));

jm.out <- find_diff(iseurat, 'QJZ', 'NEvsLum', gene.go);

#### metabolism 
myfile <- 'raw_data/Metabolism_20190418.xlsx';
meta.markers <- list{};
for(i in seq(4)){
	ifile <- read.xls(myfile, sheet = i, header = TRUE);
};
names(meta.markers) <- c('Glycolysis', 'TCA', 'OxyPhospho', 'Mitochondria');
#####
d