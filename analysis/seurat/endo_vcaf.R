library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(UpSetR);
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
comm.e <- comm[grepl('lum|CellCycle|Basal', comm$a)&grepl('Endo|Fibro|vCAF', comm$b)&!comm$b%in%c('Endothelial6', 'Endothelial7'), ];
pair.v <- unique(comm.e[grepl('vCAF', comm.e$b), ]$interacting_pair);
pair.e <- unique(comm.e[grepl('Endo', comm.e$b), ]$interacting_pair);
pair.f <- unique(comm.e[grepl('Fibro', comm.e$b), ]$interacting_pair);
venn.diagram(main='',
        list('Endothelial' = pair.e,
        'Fibroblast' = pair.f,
        'vCAF' = pair.v),
        fill=c("orange","chartreuse4","darkorchid4"),
        filename = generate.filename('venn_comm', 'vcaf_epithelia_R', 'tiff'),
        alpha=c(0.5,0.5, 0.5),
        col=c("orange","chartreuse4","darkorchid4"),
        main.cex=2,
        imagetype = "tiff",
        cex = 2
        );
pair.uv <- pair.v[!pair.v%in%c(pair.f, pair.e)];
pair.uv <- data.frame(l = gsub(' >.*', '', pair.uv), r = gsub('.*> ', '', pair.uv));
links <- data.frame(from = pair.uv$l, to = pair.uv$r);
net <- graph_from_data_frame(d = links, directed = TRUE);
l <- do.call(layout_in_circle, list(net));
#plot(net, edge.arrow.size = .4, layout = l);
#dev.off();
V(net)$label.color <- ifelse(!names(V(net))%in%pair.uv$l, '#117F2D', 'tomato');
pdf(generate.filename('network', 'vCAF_uniq_in', 'pdf'));
plot(net, edge.arrow.size = .4, layout = l, vertex.color = 'white');
legend('bottomleft', c('vCAF', 'Epithelial'), pch = 21, pt.bg = c('#117F2D', 'tomato'), 
        pt.cex = 2, bty = 'n', ncol = 3);
dev.off();
#### out communication
comm.e <- comm[grepl('lum|CellCycle|Basal', comm$b)&grepl('Endo|Fibro|vCAF', comm$a)&!comm$b%in%c('Endothelial6', 'Endothelial7'), ];
pair.v <- unique(comm.e[grepl('vCAF', comm.e$a), ]$interacting_pair);
pair.e <- unique(comm.e[grepl('Endo', comm.e$a), ]$interacting_pair);
pair.f <- unique(comm.e[grepl('Fibro', comm.e$a), ]$interacting_pair);
venn.diagram(main='',
        list('Endothelial' = pair.e,
        'Fibroblast' = pair.f,
        'vCAF' = pair.v),
        fill=c("orange","chartreuse4","darkorchid4"),
        filename = generate.filename('venn_comm', 'vcaf_epithelia_L', 'tiff'),
        alpha=c(0.5,0.5, 0.5),
        col=c("orange","chartreuse4","darkorchid4"),
        main.cex=2,
        imagetype = "tiff",
        cex = 2
        );
####
pair.uv <- pair.v[!pair.v%in%c(pair.f, pair.e)];
pair.uv <- data.frame(l = gsub(' >.*', '', pair.uv), r = gsub('.*> ', '', pair.uv));
links <- data.frame(from = pair.uv$l, to = pair.uv$r);
net <- graph_from_data_frame(d = links, directed = TRUE);
l <- do.call(layout_in_circle, list(net));
plot(net, edge.arrow.size = .4, layout = l);
dev.off();
V(net)$label.color <- ifelse(names(V(net))%in%pair.uv$l, '#117F2D', 'tomato');
pdf(generate.filename('network', 'vCAF_uniq', 'pdf'));
plot(net, edge.arrow.size = .4, layout = l, vertex.color = 'white');
legend('bottomleft', c('vCAF', 'Epithelial'), pch = 21, pt.bg = c('#117F2D', 'tomato'), 
        pt.cex = 2, bty = 'n', ncol = 3);
dev.off();
####
pair.v0 <- unique(comm.e[grepl('vCAF0', comm.e$a), ]$interacting_pair);
pair.v1 <- unique(comm.e[grepl('vCAF1', comm.e$a), ]$interacting_pair);
pair.v2 <- unique(comm.e[grepl('vCAF2', comm.e$a), ]$interacting_pair);
pair.v3 <- unique(comm.e[grepl('vCAF3', comm.e$a), ]$interacting_pair);
pair.v4 <- unique(comm.e[grepl('vCAF4', comm.e$a), ]$interacting_pair);
pair.v5 <- unique(comm.e[grepl('vCAF5', comm.e$a), ]$interacting_pair);
###

to.plot <- data.frame(matrix(0, nrow = length(pair.v), ncol = 4));
rownames(to.plot) <- pair.v;
colnames(to.plot) <- paste0('C', seq(2, 5));
to.plot[pair.v2, 1] <- 1;
to.plot[pair.v3, 2] <- 1;
to.plot[pair.v4, 3] <- 1;
to.plot[pair.v5, 4] <- 1;
to.plot[is.na(to.plot)] <- 0;
pdf(generate.filename('overlap_comm', name, 'pdf'), width = 10);
upset(to.plot, nsets = ncol(to.plot), nintersects = NA, keep.order = TRUE, 
        sets = rev(colnames(to.plot)));
dev.off();


venn.diagram(main='',
        list('C2' = pair.v2,
        'C3' = pair.v3,
        'C4' = pair.v4,
        'C5' = pair.v5),
        fill=default.colours(4),
        filename = generate.filename('venn_comm', 'vcaf_L', 'tiff'),
        alpha=c(0.5,0.5, 0.5, 0.5),
        col=default.colours(4),
        main.cex=2,
        imagetype = "tiff",
        cex = 2
        );

####
plot.dat <- data.frame(C0 = rowMeans(iseurat@data[, iseurat@ident==0]), C1 = rowMeans(iseurat@data[, iseurat@ident==1]),
        C2 = rowMeans(iseurat@data[, iseurat@ident==2]), C3 = rowMeans(iseurat@data[, iseurat@ident==3]),
        C4 = rowMeans(iseurat@data[, iseurat@ident==4]), C5 = rowMeans(iseurat@data[, iseurat@ident==5])
        );
plot.gene <- c('S100A4', 'PDGFRA', 'TNC', 'IL6', 'VIM', 'PDPN', 'TGFB3', 'MMP11', 'ACTA2', 'MYL9', 'PDGFA', 'PDGFRB', 
        'MYLK', 'SPARC', 'THY1', 'MCAM', 'FAP', 'MMP2', 'EDNRB', 'PECAM1', 'POSTN', 'SPP1', 'CSPG4', 'ENG', 'TAGLN');
to.plot <- t(scale(t(plot.dat[plot.gene, ])))
create.heatmap(
        to.plot,
        same.as.matrix = TRUE,
        xaxis.lab = NA,
        yaxis.lab = NA,
        width = 5,
        filename = generate.filename('markers', name, 'pdf')
        );
pdf(generate.filename('markers_dotmap', name, 'pdf'), width = 12, height = 5);
DotPlot(iseurat, genes.plot = plot.gene, cols = c('lightgrey', 'blue'), dot.scale = 6,
        x.lab.rot = TRUE, plot.legend = TRUE);
dev.off();
####
jm.out <- FindAllMarkers(iseurat, logfc.threshold = 0);
#saveRDS(jm.out, file = generate.filename('diff', name, 'rds'));
saveRDS(jm.out, file = generate.filename('diff_all', name, 'rds'));
jm.out <- jm.out[abs(jm.out$avg_logFC)>0.25, ]; 
to.plot <- t(scale(t(plot.dat[rownames(plot.dat)%in%rownames(jm.out), ])));
create.heatmap(
        to.plot,
        same.as.matrix = TRUE,
        xaxis.lab = NA,
        yaxis.lab = NULL,
        xaxis.rot = 0,
        width = 4,
        filename = generate.filename('diff', name, 'pdf')
        );
#iseurat <- ScaleData(iseurat);
#to.plot <- iseurat@scale.data[rownames(iseurat@scale.data)%in%rownames(jm.out[jm.out$p_val_adj<0.05, ]), ]
to.plot <- plot.dat[unique(jm.out$gene), ];
myclust <- as.dendrogram(hclust(dist(t(to.plot))));
myclust <- color_labels(myclust, k = 4);
myclust <- color_branches(myclust, k = 4);
pdf(generate.filename('average_profile', 'celltypes', 'pdf'), width = 4);
par(mar = c(3,1,1,10));
plot(myclust, horiz = TRUE);
dev.off();
####
#plot CAF gene
gene.marker <- c('TNC','PDGFRA','PDPN','TGFB3','MMP11','MYL9','PDGFA','MYLK','SPARC','MCAM','MMP2','FAP');

####
genes.lr <- unique(c(gsub(' .*', '', comm$interacting_pair), gsub('.* ', '', comm$interacting_pair)));
jup <- jm.out[jm.out$avg_logFC>0&jm.out$p_val_adj<0.05&jm.out$gene%in%genes.lr, ];
write.csv(jup, generate.filename('diff_genes', 'LR', 'csv'));
jup <- jup[jup$pct.1>jup$pct.2, ];
jup <- jup[jup$avg_logFC>1|jup$pct.2<0.1, ];
g0 <- jup[jup$cluster==0, ]$gene;
g1 <- jup[jup$cluster==1, ]$gene;
g2 <- jup[jup$cluster==2, ]$gene;
g3 <- jup[jup$cluster==3, ]$gene;
g4 <- jup[jup$cluster==4, ]$gene;
g5 <- jup[jup$cluster==5, ]$gene;

g.e <- intersect(g0, g1);
g.v1 <- g2;
g.v2 <- intersect(g3, g4);
g.v3 <- g5;

g.e <- g.e[!g.e%in%c(g.v3, g.v1, g.v2)];
g.v1 <- g.v1[!g.v1%in%c(g.e, g.v3, g.v2)];
g.v2 <- g.v2[!g.v2%in%c(g.e, g.v1, g.v3)];
g.v3 <- g.v3[!g.v3%in%c(g.e, g.v1, g.v2)];

g.all <- c(g.e, g.v1, g.v2, g.v3)
pdf('tmp.pdf', height = 40);
        Plot(iseurat, features.plot = g.all, nCol = 2);
dev.off();
####
intersect(g.e, (mgene[mgene$module%in%gsub('Module_', '', rownames(module.e)), ]$id));
intersect(g.v1, (mgene[mgene$module%in%gsub('Module_', '', rownames(module.v1)), ]$id));
intersect(g.v2, (mgene[mgene$module%in%gsub('Module_', '', rownames(module.v2)), ]$id));
intersect(g.v3, (mgene[mgene$module%in%gsub('Module_', '', rownames(module.v3)), ]$id));
####
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
#####
gene.c01 <- unique(jm.out[jm.out$cluster%in%c(0, 1)&jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]$gene);
go.c01 <- test_enrich(annot.go[grep('angiogenesis|vasculature', annot.go$name), ], gene.c01, mybg = NULL);

gene.c2 <- unique(jm.out[jm.out$cluster%in%c(2)&jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]$gene);
go.c2 <- test_enrich(annot.go[grep('angiogenesis|vasculature', annot.go$name), ], gene.c2, mybg = NULL);

gene.c34 <- intersect(jm.out[jm.out$cluster==3&jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]$gene, jm.out[jm.out$cluster==4&jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]$gene)
gene.c34 <- unique(jm.out[jm.out$cluster%in%c(3, 4)&jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]$gene);
go.c34 <- test_enrich(annot.go[grep('angiogenesis|vasculature', annot.go$name), ], gene.c34, mybg = NULL);

gene.c5 <- unique(jm.out[jm.out$cluster%in%c(5)&jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]$gene);
go.c5 <- test_enrich(annot.go[grep('angiogenesis|vasculature', annot.go$name), ], gene.c5, mybg = NULL);
###
gs <- list();
for(i in unique(annot.go[grep('angiogenesis|vasculature', annot.go$name), ]$ont)){
        gs[[i]] <- annot.go[annot.go$ont==i, ]$gene;
}
iseurat@meta.data$cluster <- as.numeric(as.vector(iseurat@meta.data$cluster_new)) + 1;
myresult <- run_qusage(iseurat, paste0(name, '_curated_brca'), gs);
saveRDS(myresult, file = generate.filename(paste0('qusage_', name), 'curated_brca', 'rds'));

enrich.keg <- read.table('Endothelial/QuSAGE_for_heatmap/KEGG_human_20190613_heatmap_input.txt');
enrich.41 <- read.table('Endothelial/QuSAGE_for_heatmap/41_human_GeneSets_heatmap_input.txt');
dim(enrich.41[rowSums(enrich.41[, 1:2]>0.1)==2&rowSums(enrich.41[, 3:6]>0)==0, ]);
####
comm.v <- comm[grepl('vCAF2', comm$a)&grepl('vCAF0|vCAF1|vCAF3|vCAF4|vCAF5', comm$b), ];
pair.v20 <- unique(comm.v[grepl('vCAF0', comm.v$b), ]$interacting_pair);
pair.v21 <- unique(comm.v[grepl('vCAF1', comm.v$b), ]$interacting_pair);
pair.v23 <- unique(comm.v[grepl('vCAF3', comm.v$b), ]$interacting_pair);
pair.v24 <- unique(comm.v[grepl('vCAF4', comm.v$b), ]$interacting_pair);
pair.v25 <- unique(comm.v[grepl('vCAF5', comm.v$b), ]$interacting_pair);

comm[comm$b=='vCAF2'&grepl('luminal', comm$a)&grepl('BAFFR|CD40|LTBR|FN14|RANK', comm$interacting_pair), ];
