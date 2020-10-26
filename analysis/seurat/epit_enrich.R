library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(readxl);
library(reshape);
library(viridis);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/epitumour");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');

name <- 'epi';
jm.out <- readRDS("epiTumor.basal.markers.3types.batchEffectRemovedMatirx.seuratClustered.rds");
jm.out <- jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0&jm.out$pct.1>0.5, ];
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
ego.b <- test_enrich(annot.keg, jm.out[jm.out$cluster=='basal/intermediate', ]$gene, mybg = NULL);
ego.c <- test_enrich(annot.keg, jm.out[jm.out$cluster=='cell-cycle-enhanced luminal', ]$gene, mybg = NULL);
ego.l <- test_enrich(annot.keg, jm.out[jm.out$cluster=='luminal', ]$gene, mybg = NULL);
grep_data <- function(myego, npath = 5){
        my.bg <- data.frame(matrix(NA, nrow = 0, ncol = length(myego)));
        colnames(my.bg) <- names(myego);
        my.or <- my.bg;
        for(i in seq(length(myego))){
                ego <- myego[[i]];
                #ego <- ego[ego$FDR<0.05, ];
                if(nrow(ego)>0){
                        ego <- ego[order(-ego$Enrichment), ];
                        bg.dat <- -log10(ego[1:min(npath, nrow(ego)), 'FDR', drop = FALSE]);
                        or.dat <- ego[1:min(npath, nrow(ego)), 'Enrichment', drop = FALSE];
                        rownames(bg.dat) <- ego[1:min(npath, nrow(ego)), ]$TermName;
                        rownames(or.dat) <- rownames(bg.dat);
                        my.bg[rownames(bg.dat), names(myego)[i]] <- bg.dat[, 1];
                        my.or[rownames(or.dat), names(myego)[i]] <- or.dat[, 1];

                };
        };
      for(i in seq(length(myego))){
            ego <- myego[[i]];
            my.bg[, names(myego)[i]] <- -log10(ego[match(rownames(my.bg), ego$TermName), 'FDR'])
            my.or[, names(myego)[i]] <- ego[match(rownames(my.or), ego$TermName), 'Enrichment'];
      }
        return(list(bg = my.bg, or = my.or))
};
myego <- list(basal = ego.b, cycle = ego.c, luminal = ego.l);
to.plot <- grep_data(myego);
saveRDS(to.plot, file = generate.filename('enrich_type', name, 'rds'))

width <- 6;
height <- 8;
xrot <- 0

####
myego.m <- list();
for(i in names(myego)){
	iego <- myego[[i]];
	iego <- iego[grep('PATH:00|PATH:01', rownames(iego)), ];
	iego <- iego[iego$P.val<0.05&iego$FDR<0.25, ];
	myego.m[[i]] <- iego;
	print(nrow(iego))
};
to.plot <- grep_data(myego.m, 10);
plot_enrich(to.plot, paste0(name, '_metab'), 6, 6, xrot = 45, cut1 = 0.25, dot.key = TRUE, colourkey = TRUE);

annot_pathgene <- function(pathway, genelist, annot.keg = annot.keg){
	myannot <- data.frame(matrix(NA, ncol = 2, nrow = length(genelist)));
	colnames(myannot) <- c('gene', 'pathway');
	myannot$gene <- genelist;
	for(i in pathway){
		pathgene <- annot.keg[annot.keg$name==i, ]$gene;
		myannot[myannot$gene%in%pathgene, ]$pathway <- paste0(myannot[myannot$gene%in%pathgene, ]$pathway, ',', i);
	};
	myannot <- na.omit(myannot);
	myannot$pathway <- gsub('^NA,', '', myannot$pathway);
	return(myannot);
};

annot.gene <- annot_pathgene(rownames(to.plot[[1]]), jm.out$gene, annot.keg);
saveRDS(annot.gene, file = generate.filename('metabolism_gene', name, 'rds'));
#### calculate qusage score 
iseurat <- readRDS(conf$seurat_epi13);
iseurat@meta.data$cluster <- as.numeric(as.factor(iseurat@meta.data$cell.type));
iseurat.13 <- iseurat
gs <- list();
for(i in unique(annot.keg$name)){
	gs[[i]] <- annot.keg[annot.keg$name==i, ]$gene;
};
myresult <- run_qusage(iseurat, paste0(name, '_keg'), gs);
saveRDS(myresult, file = generate.filename('qusage_keg', name, 'rds'));
###
mypath <- myresult[[1]];
mypath.fc <- cast(mypath, pathway.name~Cluster, mean, value = 'log.fold.change');
mypath.p <- cast(mypath, pathway.name~Cluster, mean, value = 'FDR');
mypath.fc <- mypath.fc[mypath.fc$pathway.name%in%mypath.p[rowSums(mypath.p[, -1]<0.05)>0, ]$pathway.name, ];
create.heatmap(
	x = mypath.fc[, -1], 
    cluster.dimensions = 'row',
    same.as.matrix = TRUE,
    xaxis.lab = NULL,
    yaxis.lab = NULL,
    width = 5,
    filename = generate.filename('qusage_keg', name, 'pdf')
	);
#to.plot <- mypath.fc[mypath.fc$pathway.name%in%annot.keg[grep('PATH:00|PATH:01', annot.keg$ont), ]$name, ];
to.plot <- mypath.fc
to.plot <- to.plot[rowSums(to.plot[, -1]>0.1)>0, ];
create.heatmap(
	x = to.plot[, -1], 
    cluster.dimensions = 'row',
    same.as.matrix = TRUE,
    xaxis.lab = c('basal', 'cycle', 'luminal'),
    xat = seq(3),
    yaxis.lab = to.plot$pathway.name,
    yat = seq(nrow(to.plot)),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    xaxis.rot = 0,
    width = 8,
    height = 10,
    colourkey.cex = 1.5,
    filename = generate.filename('qusage_keg_0.1', name, 'pdf')
	);
####
iseurat <- readRDS(conf$seurat_epi4);
name <- 'epi4crpc';
iseurat@meta.data$cluster <- as.numeric(as.factor(iseurat@meta.data$cell.type));
myresult <- run_qusage(iseurat, paste0(name, '_keg'), gs);
saveRDS(myresult, file = generate.filename('qusage_keg', name, 'rds'));
mypath <- myresult[[1]];
mypath.fc <- cast(mypath, pathway.name~Cluster, mean, value = 'log.fold.change');
mypath.p <- cast(mypath, pathway.name~Cluster, mean, value = 'FDR');
mypath.fc <- mypath.fc[mypath.fc$pathway.name%in%mypath.p[rowSums(mypath.p[, -1]<0.05)>0, ]$pathway.name, ];
to.plot <- mypath.fc
to.plot <- to.plot[rowSums(to.plot[, -1]>0.1)>0, ];
create.heatmap(
	x = to.plot[, -1], 
    cluster.dimensions = 'row',
    same.as.matrix = TRUE,
    xaxis.lab = c('basal', 'luminal'),
    xat = seq(3),
    yaxis.lab = to.plot$pathway.name,
    yat = seq(nrow(to.plot)),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    xaxis.rot = 0,
    width = 8,
    height = 10,
    colourkey.cex = 1.5,
    filename = generate.filename('qusage_keg_0.1', 'epi4crpc', 'pdf')
	);
#### cell communication
comm <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pvalues_0.05_communication.txt', header = TRUE, sep = '\t', as.is = TRUE);
#comm.e <- comm[grepl('lum|CellCycle|Basal', comm$a)&grepl('^CD|TAM|DC|Mono', comm$b), ];
comm.e <- comm[grepl('lum|CellCycle|Basal', comm$a)&grepl('^CD|TAM', comm$b), ];
pair.l <- unique(comm.e[grepl('luminal', comm.e$a), ]$interacting_pair);
pair.b <- unique(comm.e[grepl('Basal', comm.e$a), ]$interacting_pair);
pair.c <- unique(comm.e[grepl('CellCycle', comm.e$a), ]$interacting_pair);

pair.ll <- unique(gsub(' .*', '', pair.l));
pair.lb <- unique(gsub(' .*', '', pair.b));
pair.lc <- unique(gsub(' .*', '', pair.c));

comm.e$group <- 'neut';
comm.e[grepl('Cluster2|Cluster3|Cluster5|TAM-M1', comm.e$b), ]$group <- 'pro';
comm.e[grepl('Cluster6|TAM-M2', comm.e$b), ]$group <- 'anti';
