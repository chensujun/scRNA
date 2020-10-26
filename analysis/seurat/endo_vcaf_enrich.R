library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(readxl);
library('viridis');
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/vCAF');
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
comm <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pvalues_0.05_communication.txt', header = TRUE, sep = '\t', as.is = TRUE);
mytype <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data/objects/forCellPhone-Novel-final-version.txt', header = TRUE);
seurat.all <- readRDS(conf$sseurat_all);
iseurat <- readRDS(conf$seurat_endo);
iseurat@meta.data <- seurat.all@meta.data[rownames(iseurat@meta.data), ];
iseurat@meta.data$cluster_new <- iseurat@ident;
name <- 'vCAF';
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
####
flist <- list.files('Endothelial/GOAnalysis_result/', pattern = 'xlsx', full.names = TRUE);
for(i in grep('BP', flist)){
	ifile.b <- data.frame(read_excel(flist[i]));
	ifile.c <- data.frame(read_excel(flist[(i+1)]));
	ifile.m <- data.frame(read_excel(flist[(i+2)]));
	ifile <- rbind(cbind(ifile.b, type = 'BP'),
		cbind(ifile.c, type = 'CC'), cbind(ifile.m, type = 'MF'))	
	assign(paste0('mygo.', gsub('\\..*', '', gsub('.*/', '', flist[i]))), ifile)
};
####
flist <- list.files('Endothelial/PathWayAnalysis_result/', pattern = 'xlsx', full.names = TRUE);
for(i in seq(length(flist))){
	ifile <- data.frame(read_excel(flist[i]));
	assign(paste0('mykegg.', gsub('\\..*', '', gsub('.*/', '', flist[i]))), ifile)
};
####
myego <- list(mykegg.0, mykegg.1, mykegg.2, mykegg.3, mykegg.4, mykegg.5);
names(myego) <- paste0('C', seq(0, 5));
grep_data <- function(myego){
	my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(myego)));
	colnames(my.bg) <- names(myego);
	my.or <- my.bg;
	for(i in seq(length(myego))){
		ego <- myego[[i]];
		ego <- ego[ego$FDR<0.05, ];
		if(nrow(ego)>0){
			ego <- ego[order(-ego$Enrichment), ];
			bg.dat <- -log10(ego[1:min(5, nrow(ego)), 'FDR', drop = FALSE]);
			or.dat <- ego[1:min(5, nrow(ego)), 'Enrichment', drop = FALSE];
			rownames(bg.dat) <- ego[1:min(5, nrow(ego)), ]$PathwayTerm;
			rownames(or.dat) <- rownames(bg.dat);
			my.bg[rownames(bg.dat), names(myego)[i]] <- bg.dat[, 1];
			my.or[rownames(or.dat), names(myego)[i]] <- or.dat[, 1];

		};
	};
      for(i in seq(length(myego))){
            ego <- myego[[i]];
            my.bg[, names(myego)[i]] <- -log10(ego[match(rownames(my.bg), ego$PathwayTerm), 'FDR'])
            my.or[, names(myego)[i]] <- ego[match(rownames(my.or), ego$PathwayTerm), 'Enrichment'];
      }
	return(list(bg = my.bg, or = my.or))
};

to.plot <- grep_data(myego);
saveRDS(to.plot, generate.filename('enrich_cluster', name, 'rds'));

my.or <- to.plot[['or']];
my.bg <- to.plot[['bg']];
width <- 6;
height <- 8;
xrot <- 0
create.dotmap(
      file = generate.filename(paste0('dotmap_kegg'), name, 'pdf'),
      x = my.or,
      xaxis.cex = 1,
      yaxis.cex = 1.2,
      left.padding = 0,
      bottom.padding = 4,
      # use specified spot size and colour functions
      spot.size.function = spot.size.function,
      spot.colour.function = spot.colour.function,
      # create a legend matching the dot sizes
      #key = dot.key,
      key.top = 1,
      xaxis.lab = gsub('JD1800|SL', '', colnames(my.bg)),
      yaxis.lab = rownames(my.or),
      xaxis.rot = xrot,
      pch = 21,
      pch.border.col = 'transparent',
      # add the background
      bg.data = my.bg,
      # add a colourkey
      colourkey = FALSE,
      colour.scheme = c("white", "black"),
      total.colour = 5,
      bg.alpha = 1,
      at = c(0, -log10(0.1), -log10(0.01), 5),
      colourkey.labels.at = c(0, -log10(0.1), 2, 5),
      colourkey.labels = c(1, expression(0.1), expression(10^-2), expression(''<=10^-5)),
      width = width,
      height = height,
      na.spot.size = 3,
      add.grid = TRUE,
      col.lwd = 1,
      style = 'Nature',
      col.colour = 'black', 
      row.colour = 'black', 
      );
####
jm.out <- readRDS('2019-08-15_diff_all_vCAF.rds');
jm.out <- jm.out[abs(jm.out$avg_logFC)>0.25, ]; 
plot.gene <- intersect(annot.keg[annot.keg$name%in%c('ECM-receptor interaction', 'Focal adhesion'), ]$gene, jm.out$gene);
plot.gene <- unique(jm.out[jm.out$gene%in%plot.gene&jm.out$avg_logFC>log(2), ]$gene);
to.plot <- data.frame(t(scale(t(iseurat@data[plot.gene, ]))));
to.plot <- to.plot[, order(iseurat@ident)];
to.plot[to.plot>1] <- 1;
to.plot[to.plot<(-1)] <- -1;
col.col <- rep('transparent', ncol(to.plot));
col.col[sapply(seq(5), function(x) sum(table(iseurat@ident)[1:x]))] <- 'black';
create.heatmap(
        to.plot,
        cluster.dimensions = 'row',
        same.as.matrix = TRUE,
        xaxis.lab = NULL,
        yaxis.lab = NULL,
        width = 5,
        col.lines = sapply(seq(5), function(x) sum(table(iseurat@ident)[1:x])),
        filename = generate.filename('path_ECM_focal', name, 'pdf')
        );
###
plot.dat <- data.frame(C0 = rowMeans(iseurat@data[, iseurat@ident==0]), C1 = rowMeans(iseurat@data[, iseurat@ident==1]),
        C2 = rowMeans(iseurat@data[, iseurat@ident==2]), C3 = rowMeans(iseurat@data[, iseurat@ident==3]),
        C4 = rowMeans(iseurat@data[, iseurat@ident==4]), C5 = rowMeans(iseurat@data[, iseurat@ident==5])
        );
to.plot <- t(scale(t(plot.dat[plot.gene, ])))
create.heatmap(
        to.plot,
        cluster.dimensions = 'row',
        same.as.matrix = TRUE,
        xaxis.rot = 0,
        xaxis.fontface = 'plain',
        yaxis.fontface = 'plain',
        xaxis.lab = NA,
        yaxis.lab = NA,
        width = 5,
        #filename = generate.filename('path_ECM_focal', paste0(name, '_ave'), 'pdf')
        filename = generate.filename('path_sel', paste0(name, '_ave'), 'pdf')
        );
###
ipath <- 'NF-kappa B signaling pathway';
plot.gene <- intersect(annot.keg[annot.keg$name==ipath, ]$gene, rownames(iseurat@data));
to.plot <- na.omit(t(scale(t(plot.dat[plot.gene, ]))));
create.heatmap(
        to.plot,
        cluster.dimensions = 'row',
        same.as.matrix = TRUE,
        xaxis.rot = 0,
        xaxis.fontface = 'plain',
        yaxis.fontface = 'plain',
        xaxis.lab = NA,
        yaxis.lab = NA,
        width = 5,
        height = 10,
        #filename = generate.filename('path_ECM_focal', paste0(name, '_ave'), 'pdf')
        filename = generate.filename('path_sel', paste0('NFKB', '_ave'), 'pdf')
        );
#### NF-KB pathway
path.nf <- read.xls('objects/pathway_NFKB.xlsx', header = FALSE);
gs.nf <- list(
	nonc = as.vector(path.nf$V1),
	cano = annot.keg[annot.keg$name==ipath, ]$gene
	);
iseurat@meta.data$cluster <- as.numeric(as.vector(iseurat@meta.data$cluster_new)) + 1;
myresult <- run_qusage(iseurat, paste0(name, '_nf'), gs.nf);
gs.ecm <- list(
	ecm = annot.keg[annot.keg$name=='ECM-receptor interaction', ]$gene,
	focal = annot.keg[annot.keg$name=='Focal adhesion', ]$gene
	);
myresult <- run_qusage(iseurat, paste0(name, '_ecm'), gs.ecm);

flist <- list();
for(i in 1:number.clusters){
	flist[[i]], paste0('mygo.'), i)
}