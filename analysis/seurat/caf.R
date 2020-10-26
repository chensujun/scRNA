library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(readxl);
library(dendextend);
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/fibroblast');
source('~/svn/singleCell/myfunctions/test_kegg.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
iseurat <- readRDS(conf$seurat_fibro);
iseurat@meta.data$orig.ident <- seurat.all@meta.data[rownames(iseurat@meta.data), ]$orig.ident;
name <- 'fibro';
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
pdf(generate.filename('plottsne_sample', name, 'pdf'), width = 8, height = 8);
DimPlot(iseurat, reduction.use = 'tsne', group.by = 'orig.ident', pt.size = 2, no.axes = TRUE, vector.friendly = TRUE);
dev.off();
pdf(generate.filename('plottsne_cluster', name, 'pdf'), width = 8, height = 8);
DimPlot(iseurat, reduction.use = 'tsne', group.by = 'ident', pt.size = 2, no.axes = TRUE, vector.friendly = TRUE);
dev.off();

mywidth <- 16;
myheight <- 12;
myres <- 300;
tiff(generate.filename('plottsne', paste0(name, '_cluster'), 'tiff'), width = myheight, height = myheight, res = myres, units = 'in')
DimPlot(iseurat, reduction.use = 'tsne', pt.size = 3, no.axes = TRUE, do.label = TRUE, no.legend = TRUE,  
	cols.use = viridis(length(unique(iseurat@meta.data$res.0.8))), label.size = 24, group.by = 'res.0.8');
dev.off();
####
genes.marker <- c('FAP', 'S100A4', 'SPARC', 'ACTA2', 'PDGFRA', 'PDGFRB', 'CAV1', 'VIM');
pdf(generate.filename('plottsne_genes', name, 'pdf'), height = 3, width = 5);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 5, vector.friendly = TRUE, 
	no.axes = TRUE, no.legend = TRUE, nCol = 4, max.cutoff = 5);
dev.off();
pdf(generate.filename('plottsne_genes_legend', name, 'pdf'), height = 3, width = 5);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 5, vector.friendly = TRUE, 
	no.axes = TRUE, no.legend = FALSE, nCol = 4, max.cutoff = 5);
dev.off();

#####
jm.out <- FindAllMarkers(iseurat, only.pos = TRUE);
saveRDS(jm.out, generate.filename('diff_cluster', name, 'rds'));
exp.avg <- sapply(seq(0, 6), function(x) rowMeans(iseurat@data[, iseurat@ident==x]));
colnames(exp.avg) <- paste0('C', seq(0, 6));
myclust <- hclust(dist(t(exp.avg[unique(jm.out$gene), ])));
mytype <- cutree(myclust, 3);
myclust <- as.dendrogram(myclust);
myclust <- color_labels(myclust, k = 3, col = c('#FFC000', '#00B050', '#00B0F0'));
myclust <- color_branches(myclust, k = 3, col = c('#FFC000', '#00B050', '#00B0F0'));
pdf(generate.filename('average_profile_subtype', name, 'pdf'), width = 4);
par(mar = c(3,1,1,10));
plot(myclust, horiz = TRUE);
dev.off();

iseurat@meta.data$subtype <- mytype[paste0('C', iseurat@ident)];
#####
iseurat <- SetIdent(iseurat, ident.use = iseurat@meta.data$subtype);
jm.out <- FindAllMarkers(iseurat, only.pos = TRUE);
saveRDS(jm.out, generate.filename('diff_subtype', name, 'rds'));
jm.out <- jm.out[jm.out$p_val_adj<0.05, ];
jm.out <- jm.out[order(-abs(jm.out$avg_logFC)), ];
for(i in seq(3)){
	idiff <- jm.out[jm.out$cluster==i, ]$gene[1:150];
	igo <- test_enrich(annot.go, idiff);
	ikeg <- test_enrich(annot.keg, idiff);
	assign(paste0('mygo', i), igo);
	assign(paste0('mykeg', i), ikeg);
};
####
filter_go <- function(igo){
	igo <- igo[igo$FDR<0.05&igo$nPath>10&igo$nPath<500, ];
	igo <- igo[order(-igo$OR), ];
	return(igo);
};
#####
myego <- list(S1 = mygo1, S2 = mygo2, S3 = mygo3);
to.plot <- grep_data(myego, cut.fdr = 0.05, Nmin = 10, Nmax = 500, nIntersect = 10);
plot_enrich(to.plot, name, width = 6, height = 6, spot.size.function = function(x) {ifelse(x>10, 2.5, abs(x)/4)}, 
	dot.key = TRUE, dotkeymax = 10, colourkey = TRUE);
#####
genes.marker <- intersect(annot.go[annot.go$name%in%rownames(to.plot[[1]])[1:3], ]$gene, jm.out[jm.out$pct.1>0.5, ]$gene);
pdf(generate.filename('plotviolin', paste0(name, '_marker_tam'), 'pdf'), width = 1.8, height = 12);
VlnPlot(iseurat, features.plot = genes.marker, point.size.use = 0, nCol = 1, cols.use = c('#00B0F0', '#00B050', '#FFC000'));
dev.off();


myego <- list(S1 = mykeg1, S2 = mykeg2, S3 = mykeg3);
to.plot <- grep_data(myego, cut.fdr = 0.05, Nmin = 0, Nmax = 500, nIntersect = 0);
###### 
###### TF
tf.m <- read.table('objects/reanalysis_mast_CAF/Fibroblast/6.SCENIC/GraphClust.aucell_output.tsv', header = TRUE, row.names = 1);
rownames(tf.m) <- gsub('\\(.*', '', rownames(tf.m));
rownames(tf.m) <- gsub('NKX3.1', 'NKX3-1', rownames(tf.m));
test.df <- data.frame(t(tf.m));
rownames(test.df) <- gsub('\\.', '-', rownames(test.df));
test.df$group <- iseurat@meta.data[rownames(test.df), ]$subtype;
test.pval <- data.frame(pval = apply(test.df[, 1:nrow(tf.m)], 2, function(x) summary(aov(x~test.df$group))[[1]][[5]][1]));
test.pval$fdr <- p.adjust(test.pval$pval);
rownames(test.pval) <- gsub('NKX3.1', 'NKX3-1', rownames(test.pval));

test.pval$S1 <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='1'], 1, median);
test.pval$S2 <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='2'], 1, median);
test.pval$S3 <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='3'], 1, median);
test.pval$S1.mean <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='1'], 1, mean);
test.pval$S2.mean <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='2'], 1, mean);
test.pval$S3.mean <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='3'], 1, mean);

test.pval$S3.pct <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='3'], 1, 
	function(x) length(x[x>0])/length(x));
test.pval$S1.pct <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='1'], 1, 
	function(x) length(x[x>0])/length(x));
test.pval$S2.pct <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='2'], 1, 
	function(x) length(x[x>0])/length(x));
test.pval[, 12:14] <- data.frame(
	auc.1 = rowMeans(tf.m[rownames(test.pval), colnames(tf.m)%in%gsub('-', '.', rownames(iseurat@meta.data[iseurat@meta.data$subtype=='1', ]))]),
	auc.2 = rowMeans(tf.m[rownames(test.pval), colnames(tf.m)%in%gsub('-', '.', rownames(iseurat@meta.data[iseurat@meta.data$subtype=='2', ]))]),
	auc.3 = rowMeans(tf.m[rownames(test.pval), colnames(tf.m)%in%gsub('-', '.', rownames(iseurat@meta.data[iseurat@meta.data$subtype=='3', ]))])	
	);

test.pval[, 15:17] <- data.frame(
	auc.med.1 = apply(tf.m[rownames(test.pval), colnames(tf.m)%in%gsub('-', '.', rownames(iseurat@meta.data[iseurat@meta.data$subtype=='1', ]))], 1, median),
	auc.med.2 = apply(tf.m[rownames(test.pval), colnames(tf.m)%in%gsub('-', '.', rownames(iseurat@meta.data[iseurat@meta.data$subtype=='2', ]))], 1, median),
	auc.med.3 = apply(tf.m[rownames(test.pval), colnames(tf.m)%in%gsub('-', '.', rownames(iseurat@meta.data[iseurat@meta.data$subtype=='3', ]))], 1, median)	
	);

test.diff <- test.pval[test.pval$fdr< 0.05, ];
test.comm <- test.pval[test.pval$fdr>0.5, ];
####
#gene.diff <- c(rownames(test.diff[order(-test.diff$auc.med.1), ])[1:5],
#	rownames(test.diff[order(-test.diff$auc.med.2), ])[1:5],
#	rownames(test.diff[order(-test.diff$auc.med.3), ])[1:5]);
#gene.diff <- gene.diff[!duplicated(gene.diff)];
gene.diff <- rownames(na.omit(test.diff[rowSums(test.diff[, 15:17]>0.1)>0, ]));
gene.diff <- c("BHLHE41", "CREB3L1", "MEOX2", "SALL3","TRPS1", "TWIST2", "DDIT3", "NFATC1");

gene.comm <- rownames(test.comm[order(-rowSums(test.comm[, 15:17])), ])[1:4];
#iseurat@meta.data[, 6:21] <- t(tf.m[c(gene.comm, gene.diff), gsub('-', '.', rownames(iseurat@meta.data))]);
#colnames(iseurat@meta.data)[6:21] <- paste0('tf_', c(gene.comm, gene.diff));
iseurat@meta.data[, 6:17] <- t(tf.m[c(gene.comm, gene.diff), gsub('-', '.', rownames(iseurat@meta.data))]);
colnames(iseurat@meta.data)[6:17] <- paste0('tf_', c(gene.comm, gene.diff));
genes.marker <- paste0('tf_', c(gene.comm, gene.diff));
pdf(generate.filename('plottsne_tf', name, 'pdf'), height = 4, width = 9);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 5, vector.friendly = TRUE, 
	no.axes = TRUE, no.legend = TRUE, nCol = 6, cols.use = rainbow(8)[6:1], max.cutoff = 0.1);
dev.off();

####
flist <- list.files('objects/reanalysis_mast_CAF/Fibroblast/7.QuSAGE/41_human_GeneSets/', pattern = 'xlsx', full.names = TRUE);
i <- 1;	
ifile <- data.frame(read_excel(flist[i]));
mypath <- ifile[, 1:2]
colnames(mypath)[2] <- paste0('C', gsub('\\..*|GeneSetsInfo_Cluster_|_41_human_GeneSets', '', gsub('.*/', '', flist[i])));
for(i in seq(2, length(flist))){
	ifile <- data.frame(read_excel(flist[i]));
	assign(paste0('mygo.', gsub('\\..*', '', gsub('.*/', '', flist[i]))), ifile);
	mypath[, i+1] <- ifile[match(mypath[, 1], ifile[, 1]), 2];
	colnames(mypath)[i+1] <- paste0('C', gsub('\\..*|GeneSetsInfo_Cluster_|_41_human_GeneSets', '', gsub('.*/', '', flist[i])));
};
rownames(mypath) <- mypath[, 1];
mypath <- mypath[, -1];
mypath <- mypath[, c(1,2,6,3,4,7,5)];
####
mscore <- data.frame(read_excel('objects/reanalysis_mast_CAF/Fibroblast/9.ScGeneModule/GraphClust.ModuleExp.xlsx'));
rownames(mscore) <- gsub(' ', '_', mscore$Module);
mscore <- mscore[, -1]
mscore <- mscore[, c(1, 2, 6, 3, 4, 7, 5)];
mymodule <- c(
	rownames(mscore[rowSums(mscore[, 1:3]>0)==3&rowSums(mscore[, 4:7]<0)==4, ]),
	rownames(mscore[rowSums(mscore[, 4:6]>0)==3&rowSums(mscore[, -c(4:6)]<0)==4, ]),
	rownames(mscore[mscore$Cluster4>0&rowSums(mscore[, -7]<0)==6, ])
	);
to.plot <- mscore[mymodule, ];
pheatmap(to.plot, cluster_row = FALSE, cluster_col = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(100),
        breaks = seq(-max(to.plot), max(to.plot), length.out = 100), border_color = NA);
dev.off();
#sscore <- data.frame(read_excel('objects/reanalysis_mast_CAF/Fibroblast/9.ScGeneModule/GraphClust.SuperModuleExp.xlsx'));
mgene <- data.frame(read_excel('objects/reanalysis_mast_CAF/Fibroblast/9.ScGeneModule/GraphClust.module.xlsx'));
test_module <- function(mymodule, mgene, annot.keg, annot.go){
	gene.mymodule <- as.vector(mgene[mgene$module%in%gsub('Module_', '', mymodule), ]$id);
	keg.mymodule <- test_enrich(annot.keg, gene.mymodule, mybg = NULL);
	go.mymodule <- test_enrich(annot.go, gene.mymodule, mybg = NULL);
	myresult <- list(keg = keg.mymodule, go = go.mymodule);
	return(myresult);
};
keg.all <- go.all <- list();
for(i in rownames(to.plot)){
	print(i);
	iresult <- test_module(i, mgene, annot.keg, annot.go);
	keg.all[[i]] <- iresult[['keg']];
	go.all[[i]] <- iresult[['go']];
};
save(keg.all, go.all, file = generate.filename('enrich_modules', name, 'rda'));
filter_term <- function(mytest){
	itest <- mytest[mytest$FDR<0.25&mytest$nPath>10&mytest$P.val<0.01, ];
	itest <- itest[order(-itest$Enrichment), ];
	return(itest);
};

for(i in names(go.all)){
	go.all[[i]] <- filter_term(go.all[[i]]);
	print(i)
	print(go.all[[i]][1:10, ])
};

for(i in names(go.all)){
	igo <- go.all[[i]]
	ifelse(min(igo$FDR<0.05), print(i), print(''))
};
#####
rag <- read.table('/.mounts/labs/cpcgene/private/projects/Hypoxia_Genomics/hypoxia_signature_scores/2016-09-06_ragnum_gene_list.txt', as.is = TRUE);
iseurat@meta.data$hypoxia <- colSums(na.omit(t(scale(t(iseurat@data[rownames(iseurat@data)%in%rag$ragnum_gene_list, ])))));

pdf(generate.filename('plotgene_hypoxia', name, 'pdf'), width = 5, height = 5);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = 'hypoxia', pt.size = 2, vector.friendly = TRUE, no.axes = TRUE, no.legend = TRUE, nCol = 1);
dev.off();

