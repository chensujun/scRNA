library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(readxl);
library(plyr);
library(viridis);
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/Figures/F4');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
iseurat <- readRDS(conf$seurat_crpc_endo);
name <- 'crpc_endo';
pdf(generate.filename('plottsne', paste0(name, '_cluster'), 'pdf'));
DimPlot(iseurat, reduction.use = 'tsne', pt.size = 1);
dev.off();
flist <- list.files('objects/rsc_190822/rscEndothelial--20190822/3.PathWayAnalysis/', pattern = 'xlsx', full.names = TRUE);
for(i in seq(length(flist))){
	ifile <- data.frame(read_excel(flist[i]));
	assign(paste0('mykegg.', gsub('\\..*', '', gsub('.*/', '', flist[i]))), ifile)
};
genes.marker <- c('PECAM1', 'CAV1', 'S100A4', 'ACTA2', 'THY1', 'EDNRB');
pdf(generate.filename('plotgene', name, 'pdf'), width = 9, height = 12);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 4, nCol = 2, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();

####
iseurat@meta.data$group <- ifelse(grepl('waizhou', iseurat@meta.data$orig.ident), 'n', 'c');
ptime <- read.table('objects/rsc_190822/rscEndothelial--20190822/4.Pseudotime/GraphClust.Pseudotime.Summary_Cell.txt', header = TRUE);
to.plot <- ddply(ptime, 'FinalCluster', numcolwise(mean));
to.plot$ncell <- as.numeric(table(iseurat@meta.data$res.0.8));
to.plot$nCRPC <- as.numeric(table(iseurat@ident[iseurat@meta.data$group=='c']));
to.plot$perc <- 100*(to.plot$nCRPC)/to.plot$ncell;
create.heatmap(
	to.plot[, 'Pseudotime', drop = FALSE],
	cluster.dimensions = 'none',
    same.as.matrix = TRUE,
    colour.scheme = c('gold', 'red'),
    xaxis.rot = 0,
    xaxis.fontface = 'plain',
    yaxis.fontface = 'plain',
    yaxis.lab = NULL,
    width = 5,
    height = 1.5,
    print.colour.key = TRUE,
    colourkey.cex = 1.5,
    file = generate.filename('Pseudotime', name, 'pdf')
	);

create.heatmap(
	to.plot[, 'perc', drop = FALSE],
	cluster.dimensions = 'none',
    same.as.matrix = TRUE,
    colour.scheme = c('#B2FFFF', '#00308F'),
    xaxis.rot = 0,
    xaxis.fontface = 'plain',
    yaxis.fontface = 'plain',
    yaxis.lab = NULL,
    width = 5,
    height = 1.5,
    print.colour.key = TRUE,
    colourkey.cex = 1.5,
    file = generate.filename('percent', name, 'pdf')
	);
####pie chart schematics
#to.schem <- data.frame(prop = matrix(0.5, nrow = 2, ncol = 1), group = c('ptime', 'percent'));
#ggplot(to.schem, aes(x="", y=prop, fill = c('blue', 'red'))) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + 
#    geom_text(aes(label = unlist(round(to.plot[1, c('Pseudotime', 'perc')], 2))), position = position_stack(vjust = 0.5)) + 
#    labs(x = NULL, y = NULL, fill = NULL, title = "Cluster 0")
mycol2 <- colorRampPalette(c('gold', 'red'))(7);
mycol1 <- colorRampPalette(c('#B2FFFF', '#00308F'))(7);
pdf(generate.filename('Pseudotime_percent', name, 'pdf'), width = 8, height = 3);
layout(matrix(seq(7), nrow = 1, ncol = 7));
for(i in seq(7)){
    par(mar = rep(0, 4));
    pie(c(1, 1), labels = round(to.plot[i, c('Pseudotime', 'perc')], 2), main = paste0('Cluster ', i - 1),
        col = c(mycol1[rank(to.plot$Pseudotime)[i]], mycol2[rank(to.plot$perc)[i]]), border = NA);
};
dev.off();

#### run qusage and get plotting data
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
gs.ecm <- list(
	ecm = annot.keg[annot.keg$name=='ECM-receptor interaction', ]$gene,
	focal = annot.keg[annot.keg$name=='Focal adhesion', ]$gene
	);
iseurat@meta.data$cluster <- as.numeric(as.vector(iseurat@meta.data$res.0.8)) + 1;
myresult <- run_qusage(iseurat, paste0(name, '_ecm'), gs.ecm);
load('2019-08-22_crpc_endo_ecm_qusage_results_2_100.rda');
mycol <- readRDS(conf$col_crpc_endo);
#mycol <- viridis(number.clusters);
max.x.rg <- 0;
min.x.rg <- 0;
max.y.rg <- 0;
for (i in 1:number.clusters ){
qs <- get(paste0("results.",i))
if (max(qs$path.mean)>max.x.rg){
  max.x.rg <- max(qs$path.mean)
}
if (min(qs$path.mean)<min.x.rg){
  min.x.rg <- min(qs$path.mean)
}
if (max(qs$path.PDF)>max.y.rg){
  max.y.rg <- max(qs$path.PDF)
}}
#Plot correlation plots by geneset
for (i in 1:length(gs)){
pdf(generate.filename(paste0('qusage_', names(gs)[i]), name, 'pdf'), width = 8, height = 4);
for (j in 1:number.clusters){
  qs <- get(paste0("results.",j))
  if (j==1){
    plotDensityCurves(qs,path.index=i,col=mycol[j],main=names(gs)[i],
    	xlim=c(min.x.rg-0.01,max.x.rg+0.01),ylim=c(0,5*ceiling(max.y.rg/5)),
    	xlab="Gene Set Activation",lwd=5,cex.main=2.5,cex.axis=1.5,cex.lab=2)
  } else {
    plotDensityCurves(qs,path.index=i,add=TRUE,col=mycol[j],lwd=5)
}}
#leg <- paste0("Cluster ",1:number.clusters)
#legend("topright",legend=leg,lty=1,col=mycol,lwd=5,cex=2, bty="n",
#	pt.cex=2, ncol = ceiling(number.clusters/10))
#box(lwd=5)
dev.off()
};
