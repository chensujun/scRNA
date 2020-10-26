library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
out <- readRDS(conf$results_mnn);
seurat.all <- readRDS(conf$sseurat_all);
mytype <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data/objects/forCellPhone-Novel-final-version.txt', header = TRUE);
iseurat <- SubsetData(seurat.all, cells.use = as.vector(mytype[grep('CD', mytype$type_epi), ]$nGene));
name <- 'tcell';
iseurat@meta.data$cluster_t <- droplevels(mytype[match(rownames(iseurat@meta.data), mytype$nGene), ]$type_epi);
mycount <- iseurat@data;
iseurat@meta.data$group <- gsub('Cluster', '', gsub('Cluster1|Cluster4', '1', iseurat@meta.data$cluster_t));
#iseurat <- run_tsne_seurat(iseurat, out, 'tcell');
#### KLK3 expression
iseurat@meta.data$cluster_t <- gsub('.*Cluster', '', iseurat@meta.data$cluster_t);
iseurat@meta.data$cluster_t <- factor(iseurat@meta.data$cluster_t, levels = c(6, 0, 1, 4, 2, 3, 5));
###
pdf(generate.filename('plotviolin_klk3', name, 'pdf'), height = 4);
VlnPlot(iseurat, features.plot = 'KLK3', group.by = 'cluster_t', point.size.use = 0) + 
 stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95);
dev.off();

pdf(generate.filename('plotviolin_IL10', name, 'pdf'), height = 4);
VlnPlot(iseurat, features.plot = 'IL10', group.by = 'cluster_t', point.size.use = 0) + 
 stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95);
dev.off();

pdf(generate.filename('plottsne_IL2RA', name, 'pdf'), height = 4);
FeaturePlot(iseurat, reduction.use = 'tsne', pt.size = .5, features.plot = 'IL2RA', cols.use = c('grey', 'red'),
        no.axes = TRUE, no.legend = FALSE, do.return = TRUE);
dev.off();

pdf(generate.filename('plottsne_PTPRC', name, 'pdf'), height = 4);
FeaturePlot(iseurat, reduction.use = 'tsne', pt.size = .5, features.plot = 'PTPRC', cols.use = c('grey', 'red'),
        no.axes = TRUE, no.legend = FALSE, do.return = TRUE);
dev.off();
pdf(generate.filename('plottsne_IL7R', name, 'pdf'), height = 4);
FeaturePlot(iseurat, reduction.use = 'tsne', pt.size = .5, features.plot = 'IL7R', cols.use = c('grey', 'red'),
        no.axes = TRUE, no.legend = FALSE, do.return = TRUE);
dev.off();
####

####
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
gs.keg <- list();
for(i in unique(annot.keg$name)){
	gs.keg[[i]] <- annot.keg[annot.keg$name==i, ]$gene;
}
iseurat@meta.data$cluster <- as.numeric(as.vector(iseurat@meta.data$cluster_t)) + 1;
myresult <- run_qusage(iseurat, paste0(name, '_kegg'), gs.keg);
####pathway curated from BRCA paper
flist <- list.files('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pathway_from_brca', full.name = TRUE);
gs <- list();
for(i in flist){
	ifile <- read.table(i, header = TRUE);
	gs[[gsub('.*/|\\..*', '', i)]] <- ifile[, 1]
};
myresult <- run_qusage(iseurat, paste0(name, '_curated_brca'), gs);
saveRDS(myresult, file = generate.filename(paste0('qusage_', name), 'curated_brca', 'rds'));
###
mypath <- myresult[[1]];
mypath <- mypath[grep('CD8|Pro|Treg|Cyto', mypath$pathway.name), ];
mypath.fc <- cast(mypath, pathway.name~Cluster, mean, value = 'log.fold.change');
rownames(mypath.fc) <- mypath.fc$pathway.name;
mypath.fc <- mypath.fc[, -1];

mypath.p <- cast(mypath, pathway.name~Cluster, mean, value = 'FDR');
rownames(mypath.p) <- mypath.p$pathway.name;
mypath.p <- mypath.p[, -1];
mypath.p <- -log10(mypath.p);
mypath.p[mypath.p==Inf] <- 22;
#
mypath.fc <- mypath.fc[, c(7, 1, 2, 5, 3, 4, 6)];
mypath.p <- mypath.p[, colnames(mypath.fc)];
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
        colours <- rep("white", length(x));
        colours[sign(x) == -1] <- default.colours(2, palette.type = "dotmap")[1]; 
        colours[sign(x) == 1] <- default.colours(2, palette.type = "dotmap")[2]; 
        return(colours);
}
dot.key <- list(
        # indicate which side of the plot the legend will appear
        space = "right",
                points = list(
                cex = spot.size.function(c(-0.5, 0, 0.5)),
                col = spot.colour.function(c(-0.5, 0, 0.5)),
                pch = 19
        ),
        # dot labels
        text = list(
                lab = as.character(c(-0.5, 0, 0.5)),
                cex = 1.5,
                adj = 1.0,
                fontface = "bold"
                )
        );

create.dotmap(
        x = mypath.fc,
        xaxis.cex = 1,
        yaxis.cex = 1,
        left.padding = 2,
        bottom.padding = 4,
        # use specified spot size and colour functions
        spot.size.function = spot.size.function,
        spot.colour.function = spot.colour.function,
        # create a legend matching the dot sizes
        key = dot.key,
        key.top = 1,
        xaxis.lab = NULL,
        xaxis.rot = 0,
        pch = 21,
        pch.border.col = 'transparent',
        # add the background
        bg.data = mypath.p,
        # add a colourkey
        colourkey = TRUE,
        colour.scheme = c("white", "black"),
        total.colour = 4,
        bg.alpha = .8,
        at = c(0, 1.3, 2, max(mypath.p)),
        colourkey.labels.at = c(1.3, 2),
        colourkey.labels = c(0.05, 0.01),
        filename = generate.filename('pathway_curated_brca', name, 'pdf'),
        width = 7,
        height = 3,
        na.spot.size = 3,
        style = 'Nature'
        );

### violnplot
###
####
p1 <- DimPlot(iseurat, reduction.use = 'tsne', group.by = 'cluster_t', pt.size = 1, no.axes = TRUE);
pdf(generate.filename('typec', name, 'pdf'), width = 7);
p1;
dev.off();
####
genes.major <- c('CD3D', 'CD79A', 'C1QA', 'CD163', 'EPCAM', 'PECAM1', 'VIM', 'DCN', 'S100A9');
pdf(generate.filename('marker_major', name, 'pdf'));
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.major, pt.size = 2,
	cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE)
dev.off();
####diff gene
iseurat <- SetIdent(iseurat, ident.use = ifelse(grepl('CD8', iseurat@meta.data$cluster_t), 'CD8', 'CD4'))
jm.out <- FindMarkers(iseurat, ident.1 = 'CD8', ident.2 = 'CD4', logfc.threshold = 0);
saveRDS(jm.out, file = generate.filename('diff', name, 'rds'));
jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
keg.up <- test_enrich(annot.keg, jup);
keg.dn <- test_enrich(annot.keg, jdn);
####
gene.k <- read.table('~/tmp/plot/KLF4_targets.human.tsv');
iseurat <- ScaleData(iseurat);
iseurat <- SetIdent(iseurat, ident.use = gsub('.*Cluster', '', iseurat@meta.data$cluster_t))
pdf(generate.filename('heatmap', 'KLF4', 'pdf'));
print(DoHeatmap(iseurat, genes.use = intersect(gene.k[gene.k$V3=='Activation', ]$V2, rownames(seurat.all@data)),
        slim.col.label = TRUE, remove.key = FALSE));
dev.off();
pdf(generate.filename('heatmap', 'target_gene', 'pdf'));
print(DoHeatmap(iseurat, genes.use = intersect(gene.reg[gene.reg$V1%in%c('NKX3-1', 'NR1H3', 'CUX2', 'HOXC6', 'KLF4'), ]$V2, rownames(seurat.all@data)),
        slim.col.label = TRUE, remove.key = FALSE));
dev.off();
### 
iseurat <- SubsetData(iseurat, cells.use = rownames(iseurat@meta.data[iseurat@meta.data$cluster_t%in%c(2,3,5), ]));
iseurat@meta.data$group <- ifelse(iseurat@data['KLK3', ]>0, 1, 0);
iseurat <- SetIdent(iseurat, ident.use = ifelse(iseurat@data['KLK3', ]>0, 1, 0));
jm.out <- FindMarkers(iseurat, ident.1 = '1', ident.2 = '0', logfc.threshold = 0);
jup <- rownames(jm.out[jm.out$p_val<0.01&jm.out$avg_logFC>0, ]);
jdn <- rownames(jm.out[jm.out$p_val<0.01&jm.out$avg_logFC<0, ]);
keg.up <- test_enrich(annot.keg, jup);
keg.dn <- test_enrich(annot.keg, jdn);
go.up <- test_enrich(annot.go, jup);
go.dn <- test_enrich(annot.go, jdn);
save(jm.out, keg.up, keg.dn, go.up, go.dn, file = generate.filename('diff_tcell', 'klk3', 'rda'));
###
iseurat <- SetIdent(iseurat, ident.use = ifelse(iseurat@meta.data$cluster_t==5, 1, 0));
jm.out <- FindMarkers(iseurat, ident.1 = '1', ident.2 = '0', logfc.threshold = 0);
jup <- rownames(jm.out[jm.out$p_val<0.01&jm.out$avg_logFC>0, ]);
jdn <- rownames(jm.out[jm.out$p_val<0.01&jm.out$avg_logFC<0, ]);
keg.up <- test_enrich(annot.keg, jup);
keg.dn <- test_enrich(annot.keg, jdn);
go.up <- test_enrich(annot.go, jup);
go.dn <- test_enrich(annot.go, jdn);
save(jm.out, keg.up, keg.dn, go.up, go.dn, file = generate.filename('diff_tcell', 'klk3', 'rda'));
####
fontsize <- theme(axis.title=element_text(size=48), axis.text.x = element_text(size = 0), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));

genes <- c('KLK2', 'KLK3', 'PMEPA1', 'ABCC4', 'NKX3-1', 'C1orf116', 'FKBP5', 'ACSL3', 'ZBTB10', 'HERC3', 'PTGER4', 'MPHOSPH9', 
        'EAF2', 'MED28', 'NNMT', 'MAF', 'GNMT', 'CENPN', 'ELL2', 'TMPRSS2');
pdf(generate.filename('AR_signature', name, 'pdf'), height = 8, width = 13);
VlnPlot(iseurat, features.plot = genes, nCol = 5, point.size.use = 0, size.x.use = 0);
dev.off();
