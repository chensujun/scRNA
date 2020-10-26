library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(readxl);
library(viridis);
library(dendextend);
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/fibroblast');
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/grep_term.R');
source('~/svn/singleCell/myfunctions/plot_enrich_list.R');
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
exp.avg <- sapply(seq(0, 4), function(x) rowMeans(iseurat@data[, iseurat@ident==x]));
colnames(exp.avg) <- paste0('C', seq(0, 4));
myclust <- hclust(dist(t(exp.avg[unique(jm.out$gene), ])));
mytype <- cutree(myclust, 3);
myclust <- as.dendrogram(myclust);
myclust <- color_labels(myclust, k = 3, col = c('#FFC000', '#00B050', '#00B0F0'));
myclust <- color_branches(myclust, k = 3, col = c('#FFC000', '#00B050', '#00B0F0'));
pdf(generate.filename('average_profile_subtype', name, 'pdf'), width = 4);
par(mar = c(3,1,1,10));
plot(myclust, horiz = TRUE);
dev.off();
####
for(i in seq(0, 5)){
	idiff <- jm.out[jm.out$cluster==i, ]$gene[1:150];
	igo <- test_enrich(annot.go, idiff);
	ikeg <- test_enrich(annot.keg, idiff);
	assign(paste0('mygo', i), igo);
	assign(paste0('mykeg', i), ikeg);
};
save(mygo0, mygo1, mygo2, mygo3, mygo4, mykeg0, mykeg1, mykeg2, mykeg3, mykeg4, file = generate.filename('diff_cluster_enrich', name, 'rda'));

myego <- list(C0 = mygo0, C1 = mygo1, C2 = mygo2, C3 = mygo3, C4 = mygo4);
to.plot <- grep_data(myego, cut.fdr = 0.05, Nmin = 10, Nmax = 500, nIntersect = 10);
plot_enrich(to.plot, name, width = 12, height = 12, spot.size.function = function(x) {ifelse(x>10, 2.5, abs(x)/4)}, 
	dot.key = TRUE, dotkeymax = 10, colourkey = TRUE);

###
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
save(mygo1, mygo2, mygo3, mykeg1, mykeg2, mykeg3, file = generate.filename('diff_subtype_enrich', name, 'rda'));

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

myego <- list(S1 = mykeg1, S2 = mykeg2, S3 = mykeg3);
to.plot <- grep_data(myego, cut.fdr = 0.05, Nmin = 10, Nmax = 500, nIntersect = 5);
plot_enrich(to.plot, name, width = 6, height = 6, spot.size.function = function(x) {ifelse(x>10, 2.5, abs(x)/4)}, 
	dot.key = TRUE, dotkeymax = 10, colourkey = TRUE);

myego <- list(S1 = mygo1, S2 = mygo2, S3 = mygo3);
to.plot <- grep_data(myego, cut.fdr = 0.05, Nmin = 10, Nmax = 500, nIntersect = 10, rankterm = 'FDR', decreasing = FALSE);
plot_enrich(to.plot, name, width = 5.5, height = 6, spot.size.function = function(x) {ifelse(x>10, 2.5, abs(x)/4)}, 
        dot.key = TRUE, dotkeymax = 10, colourkey = TRUE);
######
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
library(qusage);
annot.41 <- read.table(conf$marker_geneset);
annot.41 <- data.frame(ont = annot.41$V2, gene = annot.41$V1, name = annot.41$V2);
annot.41 <- annot.41[annot.41$gene%in%rownames(iseurat@data), ];
gs.41 <- list();
for(i in unique(annot.41$name)){
        gs.41[[i]] <- annot.41[annot.41$name==i, ]$gene;
}
iseurat@meta.data$cluster <- iseurat@meta.data$subtype;
myresult <- run_qusage(iseurat, paste0(name, '_41'), gs.41);
mypath <- myresult[[1]]
fc <- reshape::cast(mypath, pathway.name~Cluster, mean, value = 'log.fold.change');
p <- reshape::cast(mypath, pathway.name~Cluster, mean, value = 'FDR');

for(i in seq(3)){
	idiff <- jm.out[jm.out$cluster==i, ]$gene[1:150];
	igo <- test_enrich(annot.41, idiff);
	assign(paste0('mygs', i), igo);
};
myego <- list(S1 = mygs1, S2 = mygs2, S3 = mygs3);
to.plot <- grep_data(myego, cut.fdr = 0.25, Nmin = 10, Nmax = 1000, nIntersect = 1, rankterm = 'FDR', decreasing = FALSE);

######
assigned <- readRDS(conf$cellcycle);
iseurat@meta.data$phase <- assigned$phases[match(rownames(iseurat@meta.data), assigned$name)];
####
exp.avg.sub <- sapply(seq(3), function(x) rowMeans(iseurat@data[, iseurat@meta.data$subtype==x]));
colnames(exp.avg.sub) <- paste0('S', seq(3))
myref2 <- exp.avg.sub[unique(jm.out[jm.out$p_val_adj<0.05, ]$gene), ]
mydata <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-01-05_tcga_fpkm.rds');
estm <- read.table('objects/immuneEstimation.txt', header = TRUE, row.names = 1);
myref2 <- myref2[rownames(myref2)%in%rownames(mydata), ];
mycor <- apply(mydata[rownames(myref2), ], 2, function(x) cor(x, myref2[, 2]));
names(mycor) <- gsub('\\.', '-', names(mycor));
estm <- estm[names(mycor), ];
estm$cor <- mycor;
sapply(seq(6), function(x) cor(estm[, x], estm$cor));
gene.matrix <- mydata[jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>log(2)&rowSums(jm.out[, 3:4])>0.1&jm.out$cluster==1, ]$gene, ];
gene.matrix <- na.omit(data.frame(t(scale(t(gene.matrix)))));
rownames(estm) <- gsub('-', '.', rownames(estm))
estm$mean1 <- colMeans(gene.matrix[, rownames(estm)]);
gene.matrix <- mydata[jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>log(2)&rowSums(jm.out[, 3:4])>0.1&jm.out$cluster==2, ]$gene, ];
gene.matrix <- na.omit(data.frame(t(scale(t(gene.matrix)))));
rownames(estm) <- gsub('-', '.', rownames(estm))
estm$mean2 <- colMeans(gene.matrix[, rownames(estm)]);
gene.matrix <- mydata[jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>log(2)&rowSums(jm.out[, 3:4])>0.1&jm.out$cluster==3, ]$gene, ];
gene.matrix <- na.omit(data.frame(t(scale(t(gene.matrix)))));
rownames(estm) <- gsub('-', '.', rownames(estm))
estm$mean3 <- colMeans(gene.matrix[, rownames(estm)]);
####
sapply(seq(6), function(x) cor(estm[, x], estm$mean, method = 'spearman'));
####
gene.matrix <- mydata[annot.keg[annot.keg$name=='Complement and coagulation cascades', ]$gene, ];
gene.matrix <- na.omit(data.frame(t(scale(t(gene.matrix)))));
estm$ccc <- colMeans(gene.matrix[, rownames(estm)]);
gene.matrix <- mydata[annot.keg[annot.keg$name=='ECM-receptor interaction', ]$gene, ];
gene.matrix <- na.omit(data.frame(t(scale(t(gene.matrix)))));
estm$ecm <- colMeans(gene.matrix[, rownames(estm)]);
#####
pval.cor <- get.corr.key(x = estm$mean2,
        y = estm$ccc,
        label.items = c('spearman', 'spearman.p'),
        alpha.background = 0,
        key.cex = 1
        );

create.scatterplot(
        ccc~mean2,
        estm,
        cex = 0.5,
        #add.xyline = TRUE,
        xyline.lty = 2,
        ylab.label = 'Compl. and coag. signature',
        xlab.label = 'S2 signature',
        file = generate.filename('sig_complement_S2', name, 'pdf'),
        style = 'Nature',
        legend = list(
                inside = list(
                        fun = draw.key,
                        args = list(
                                key = pval.cor
                                ),
                        x = 0.35,
                        y = 0.99
                        )
                )
        );

pval.cor <- get.corr.key(x = estm$mean2,
        y = estm$Macrophage,
        label.items = c('spearman', 'spearman.p'),
        alpha.background = 0,
        key.cex = 1
        );

create.scatterplot(
        Macrophage~mean2,
        estm,
        cex = 0.5,
        #add.xyline = TRUE,
        xyline.lty = 2,
        ylab.label = 'Macrophage estimate',
        xlab.label = 'S2 signature',
        file = generate.filename('sig_macro_S2', name, 'pdf'),
        style = 'Nature',
        legend = list(
                inside = list(
                        fun = draw.key,
                        args = list(
                                key = pval.cor
                                ),
                        x = 0.35,
                        y = 0.99
                        )
                )
        );

pval.cor <- get.corr.key(x = estm$mean2,
        y = estm$Dendritic,
        label.items = c('spearman', 'spearman.p'),
        alpha.background = 0,
        key.cex = 1
        );

create.scatterplot(
        Dendritic~mean2,
        estm,
        cex = 0.5,
        #add.xyline = TRUE,
        xyline.lty = 2,
        ylab.label = 'Dendritic estimate',
        xlab.label = 'S2 signature',
        file = generate.filename('sig_dc_S2', name, 'pdf'),
        style = 'Nature',
        legend = list(
                inside = list(
                        fun = draw.key,
                        args = list(
                                key = pval.cor
                                ),
                        x = 0.35,
                        y = 0.99
                        )
                )
        );


#####
genes.marker <- intersect(annot.go[annot.go$name%in%rownames(to.plot[[1]])[1:3], ]$gene, jm.out[jm.out$pct.1>0.5, ]$gene);
pdf(generate.filename('plotviolin', paste0(name, '_marker_tam'), 'pdf'), width = 1.8, height = 12);
VlnPlot(iseurat, features.plot = genes.marker, point.size.use = 0, nCol = 1, cols.use = c('#00B0F0', '#00B050', '#FFC000'));
dev.off();


myego <- list(S1 = mykeg1, S2 = mykeg2, S3 = mykeg3);
to.plot <- grep_data(myego, cut.fdr = 0.05, Nmin = 0, Nmax = 500, nIntersect = 0);
###### 
###### TF
tf.m <- read.table('objects/13_sample_fibroblast_20190718/6.SCENIC/GraphClust.aucell_output.tsv', header = TRUE, row.names = 1);
rownames(tf.m) <- gsub('\\(.*', '', rownames(tf.m));
rownames(tf.m) <- gsub('NKX3.1', 'NKX3-1', rownames(tf.m));
test.df <- data.frame(t(tf.m));
rownames(test.df) <- gsub('\\.', '-', rownames(test.df));
test.df$group <- factor(iseurat@meta.data[rownames(test.df), ]$subtype);
test.pval <- data.frame(pval = apply(test.df[, 1:nrow(tf.m)], 2, function(x) summary(aov(x~test.df$group))[[1]][[5]][1]));
test.pval$fdr <- p.adjust(test.pval$pval);
rownames(test.pval) <- gsub('NKX3.1', 'NKX3-1', rownames(test.pval));

test.pval$S1 <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='1'], 1, median);
test.pval$S2 <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='2'], 1, median);
test.pval$S3 <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='3'], 1, median);
test.pval$S1.mean <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='1'], 1, mean);
test.pval$S2.mean <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='2'], 1, mean);
test.pval$S3.mean <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='3'], 1, mean);

test.pval$S1.pct <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='1'], 1, 
	function(x) length(x[x>0])/length(x));
test.pval$S2.pct <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='2'], 1, 
	function(x) length(x[x>0])/length(x));
test.pval$S3.pct <- apply(iseurat@data[rownames(test.pval), iseurat@meta.data$subtype=='3'], 1, 
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

test.diff <- na.omit(test.pval[test.pval$fdr< 0.05, ]);
dim(test.diff[test.diff$auc.1>test.diff$auc.2&test.diff$auc.1>test.diff$auc.3, ]);
dim(test.diff[test.diff$auc.2>test.diff$auc.1&test.diff$auc.2>test.diff$auc.3, ]);
dim(test.diff[test.diff$auc.3>test.diff$auc.2&test.diff$auc.3>test.diff$auc.1, ]);

test.comm <- test.pval[test.pval$fdr>0.5, ];
#test.diff <- test.diff[rowSums(test.diff[, grep('auc', colnames(test.diff))]>0.1)>0&rownames(test.diff)%in%jm.out$gene, ];
test.diff <- test.diff[rownames(test.diff)%in%jm.out$gene, ];
for(i in rownames(test.diff)){
	a1 <- aov(test.df[, i]~test.df$group);
	ph <- TukeyHSD(a1, 'test.df$group', conf.level = 0.95);
	assign(paste0('ph.', i), ph)
}
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
####
comm <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pvalues_0.05_communication.txt', header = TRUE, sep = '\t', as.is = TRUE);
comm$fdr <- p.adjust(comm$P.Value);
#comm.e <- comm[grepl('Fibroblast', comm$a)&!comm$b%in%c('Endothelial6', 'Endothelial7')&grepl('TAM|^CD|Monoyte', comm$b), ];
pair.1 <- unique(comm.e[grepl('Fibroblast0|Fibroblast1|Fibroblast3', comm.e$a), ]$interacting_pair);
pair.2 <- unique(comm.e[grepl('Fibroblast2', comm.e$a), ]$interacting_pair);
pair.3 <- unique(comm.e[grepl('Fibroblast4', comm.e$a), ]$interacting_pair);

#comm.e <- comm[grepl('Fibroblast', comm$a)&!comm$b%in%c('Endothelial6', 'Endothelial7')&grepl('Fibroblast', comm$b), ];
comm.e <- comm[grepl('Fibroblast', comm$a)&grepl('Fibroblast', comm$b), ];
#comm.e$fdr <- p.adjust(comm.e$P.Value);
#comm.e <- comm.e[comm.e$fdr<=0.05, ];
#comm.e <- comm;
#comm.e$a <- paste0('Fibroblast', comm.e$a);
#comm.e$b <- paste0('Fibroblast', comm.e$b);

pair.20 <- unique(comm.e[comm.e$a == 'Fibroblast2'&comm.e$b == 'Fibroblast0', ]$interacting_pair);
pair.21 <- unique(comm.e[comm.e$a == 'Fibroblast2'&comm.e$b == 'Fibroblast1', ]$interacting_pair);
pair.23 <- unique(comm.e[comm.e$a == 'Fibroblast2'&comm.e$b == 'Fibroblast3', ]$interacting_pair);
#####
pair.22 <- unique(comm.e[comm.e$a == 'Fibroblast2'&comm.e$b == 'Fibroblast2', ]$interacting_pair);
pair.24 <- unique(comm.e[comm.e$a == 'Fibroblast2'&comm.e$b == 'Fibroblast4', ]$interacting_pair);

pair.02 <- unique(comm.e[comm.e$b == 'Fibroblast2'&comm.e$a == 'Fibroblast0', ]$interacting_pair);
pair.12 <- unique(comm.e[comm.e$b == 'Fibroblast2'&comm.e$a == 'Fibroblast1', ]$interacting_pair);
pair.32 <- unique(comm.e[comm.e$b == 'Fibroblast2'&comm.e$a == 'Fibroblast3', ]$interacting_pair);

pair.22 <- unique(comm.e[comm.e$b == 'Fibroblast2'&comm.e$a == 'Fibroblast2', ]$interacting_pair);
pair.42 <- unique(comm.e[comm.e$b == 'Fibroblast2'&comm.e$a == 'Fibroblast4', ]$interacting_pair);
#####
comm.013 <- intersect(pair.02, intersect(pair.12, pair.32));
pair.s1 <- comm.013[!comm.013%in%pair.42];

comm.013r <- intersect(pair.20, intersect(pair.21, pair.23));
pair.s2 <- comm.013r[!comm.013r%in%pair.24];
####
comm.013r <- table(c(pair.20, pair.21, pair.23));
comm.013r <- names(comm.013r[comm.013r>0]);
pair.s2 <- comm.013r[!comm.013r%in%pair.24];

comm.013 <- table(c(pair.02, pair.12, pair.32));
comm.013 <- names(comm.013[comm.013>0]);
pair.s1 <- comm.013[!comm.013%in%pair.42];
####
pair.0 <- unique(comm.e[comm.e$a=='Fibroblast0'|comm.e$b=='Fibroblast0'&comm.e$a!=comm.e$b, ]$interacting_pair)
pair.1 <- unique(comm.e[comm.e$a=='Fibroblast1'|comm.e$b=='Fibroblast1'&comm.e$a!=comm.e$b, ]$interacting_pair)
pair.2 <- unique(comm.e[comm.e$a=='Fibroblast2'|comm.e$b=='Fibroblast2'&comm.e$a!=comm.e$b, ]$interacting_pair)
pair.3 <- unique(comm.e[comm.e$a=='Fibroblast3'|comm.e$b=='Fibroblast3'&comm.e$a!=comm.e$b, ]$interacting_pair)
pair.4 <- unique(comm.e[comm.e$a=='Fibroblast4'|comm.e$b=='Fibroblast4'&comm.e$a!=comm.e$b, ]$interacting_pair)

write.table(c(pair.s1, '0-2', '1-2', '3-2', '4-2'), 'Fibro_unique_pairs_S1-S2.txt', col.names = FALSE, row.names = FALSE);
write.table(c(pair.s2, '2-0', '2-1', '2-3', '2-4'), 'Fibro_unique_pairs_S2-S1.txt', col.names = FALSE, row.names = FALSE);
write.csv(comm.e, generate.filename('comm', name, 'csv'))；
####
pair.40 <- unique(comm.e[comm.e$a=='Fibroblast4'&comm.e$b=='Fibroblast0', ]$interacting_pair)
pair.41 <- unique(comm.e[comm.e$a=='Fibroblast4'&comm.e$b=='Fibroblast1', ]$interacting_pair)
pair.42 <- unique(comm.e[comm.e$a=='Fibroblast4'&comm.e$b=='Fibroblast2', ]$interacting_pair)
pair.43 <- unique(comm.e[comm.e$a=='Fibroblast4'&comm.e$b=='Fibroblast3', ]$interacting_pair)
####
comm.t <- comm
comm.t$type <- 'luminal';
comm.t[grepl('^CD', comm.t$b), ]$type <- 'tcell';
comm.t[grepl('Fibro', comm.t$b), ]$type <- 'fibro';
comm.t[grepl('TAM|Mono|DC', comm.t$b), ]$type <- 'myeloid';
comm.t[grepl('vCAF|Endo', comm.t$b), ]$type <- 'endo';
comm.t[grepl('Basal', comm.t$b), ]$type <- 'basal';
comm.t[grepl('Mast', comm.t$b), ]$type <- 'mast';

comm.t$type.a <- 'luminal';
comm.t[grepl('^CD', comm.t$a), ]$type.a <- 'tcell';
comm.t[grepl('Fibro', comm.t$a), ]$type.a <- 'fibro';
comm.t[grepl('TAM|Mono|DC', comm.t$a), ]$type.a <- 'myeloid';
comm.t[grepl('vCAF|Endo', comm.t$a), ]$type.a <- 'endo';
comm.t[grepl('Basal', comm.t$a), ]$type.a <- 'basal';
comm.t[grepl('Mast', comm.t$a), ]$type.a <- 'mast';
comm.t[grep('TAM', comm.t$a), ]$type.a <- 'TAM';
table(comm.t[comm.t$type==comm.t$type.a&!comm.t$a==comm.t$b, ]$ab);
###
to.plot <- exp.avg[intersect(annot.go[grep('extracellular matrix', annot.go$name), ]$gene, jm.out$gene), ];
to.plot <- data.frame(t(scale(t(to.plot))));
to.plot$max <- apply(to.plot, 1, function(x) grep(max(x), x));
pdf(generate.filename(paste0('avg_', name), 'ecm', 'pdf'), width = 4);
pheatmap(to.plot[, 1:5], col = colorRampPalette(c('blue', 'white', 'red'))(100), show_rownames = FALSE, angle_col = 0);
dev.off(); 

for(i in c(5, 3)){
idiff <- rownames(to.plot[to.plot$max==i, ]);
igo <- test_enrich(annot.go, idiff);
ikeg <- test_enrich(annot.keg, idiff);
assign(paste0('myecm', i-1), igo);
assign(paste0('myecm_keg', i-1), ikeg);
};

for(i in c(5, 3)){
idiff <- rownames(to.plot[to.plot$max%in%i, ]);
igo <- test_enrich(annot.go, idiff, unique(annot.go[grep('extracellular matrix', annot.go$name), ]$gene));
ikeg <- test_enrich(annot.keg, idiff, unique(annot.go[grep('extracellular matrix', annot.go$name), ]$gene));
assign(paste0('myecm', i-1), igo);
assign(paste0('myecm_keg', i-1), ikeg);
};
i <- c(1, 2, 4);
igo <- test_enrich(annot.go, idiff, unique(annot.go[grep('extracellular matrix', annot.go$name), ]$gene));
ikeg <- test_enrich(annot.keg, idiff, unique(annot.go[grep('extracellular matrix', annot.go$name), ]$gene));
assign(paste0('myecm', 1), igo);
assign(paste0('myecm_keg', 1), ikeg);
save(myecm1, myecm2, myecm4, myecm_keg1, myecm_keg2, myecm_keg4, file = generate.filename('enrich_ecm', name, 'rda'));

myego <- list(S1 = myecm1, S2 = myecm2, S3 = myecm4);
to.plot <- grep_data(myego, cut.fdr = 0.25, Nmin = 10, Nmax = 500, nIntersect = 0, rankterm = 'P.val', decreasing = FALSE);
plot_enrich(to.plot, name, width = 10, height = 6, spot.size.function = function(x) {abs(x)}, 
        dot.key = TRUE, dotkeymax = 5, colourkey = TRUE);

to.plot <- grep_data(myego, cut.fdr = 0.25, Nmin = 10, Nmax = 500, nIntersect = 0, rankterm = 'Enrichment', decreasing = TRUE);
plot_enrich(to.plot, paste0('ecm_or_', name), width = 8, height = 6, cut1 = 0.25, spot.size.function = function(x) {abs(x)/2}, 
        dot.key = TRUE, dotkeymax = 5, colourkey = TRUE);

gs.ecm <- list();
for(i in droplevels(annot.go[grepl('extracellular matrix|angiogenesis', annot.go$name)&annot.go$gene%in%rownames(iseurat@data), ]$name)){
        igene <- intersect(annot.go[annot.go$name==i, ]$gene, rownames(iseurat@data));
        #if(length(igene)>5)
        gs.ecm[[i]] <- igene;
}
iseurat@meta.data$cluster <- iseurat@meta.data$subtype;
myresult <- run_qusage(iseurat, paste0(name, '_ecm'), gs.ecm);
mypath <- myresult[[1]];
#mypath <- mypath[mypath$FDR<0.05, ];
fc <- reshape::cast(mypath, pathway.name~Cluster, mean, value = 'log.fold.change');
p <- reshape::cast(mypath, pathway.name~Cluster, mean, value = 'FDR');

rownames(fc) <- fc$pathway.name;
fc <- fc[, -1];
rownames(p) <- p$pathway.name;
p <- p[, -1];
