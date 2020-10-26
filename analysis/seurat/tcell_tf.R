library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
iseurat <- readRDS(conf$seurat_tcell);
name <- 'tcell';
tf.auc <- read.table('raw_data/TCell/7.SCENIC/GraphClust.aucell_output.tsv', header = TRUE, row.names = 1);
new.ident <- gsub('6', 'Treg', iseurat@ident);
new.ident <- gsub('0', 'Tconv', new.ident);
new.ident <- gsub('1|4', 'Naive', new.ident);
new.ident <- gsub('2|3|5', 'Effector', new.ident);
names(new.ident) <- names(iseurat@ident);
new.ident <- factor(new.ident);
iseurat@ident <- new.ident;
rownames(tf.auc) <- gsub('\\(.*', '', rownames(tf.auc));
test.df <- data.frame(t(tf.auc));
rownames(test.df) <- gsub('\\.', '-', rownames(test.df));
test.df$group <- new.ident[rownames(test.df)];
test.pval <- data.frame(pval = apply(test.df[, 1:150], 2, function(x) summary(aov(x~test.df$group))[[1]][[5]][1]));
test.pval$fdr <- p.adjust(test.pval$pval);
jm.out <- FindAllMarkers(iseurat, logfc.threshold = 0, only.pos = TRUE);
saveRDS(test.pval, file = generate.filename('subtype', 'diff_reg', 'rds'));
saveRDS(jm.out, file = generate.filename('subtype', 'diff_gene', 'rds'));
#####
diff.tf <- intersect(jm.out$gene, rownames(test.pval[test.pval$fdr<0.05, ]));
pdf(generate.filename('diff', 'tf', 'pdf'), width = 9, height = 12);
VlnPlot(iseurat, features.plot = diff.tf, nCol = 3, point.size = 0);
dev.off();
######
test.pval[, 3:6] <- data.frame(effector = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Effector'], 1, median),
	naive = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Naive'], 1, median),
	tconv = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Tconv'], 1, median),
	treg = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Treg'], 1, median)
	);
test.pval$max <- apply(test.pval[, 3:6], 1, max);
test.pval[, 8:11] <- data.frame(effector.p = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Effector'], 1, function(x) length(x[x>0])/length(x)),
	naive.p = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Naive'], 1, function(x) length(x[x>0])/length(x)),
	tconv.p = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Tconv'], 1, function(x) length(x[x>0])/length(x)),
	treg.p = apply(iseurat@data[rownames(test.pval), iseurat@ident=='Treg'], 1, function(x) length(x[x>0])/length(x))
	);
test.pval$max.p <- apply(test.pval[, 8:11], 1, max);
diff.tf <- rownames(test.pval[test.pval$max.p>0.1&test.pval$fdr<0.05, ]);
pdf(generate.filename('diff', 'tf_p10', 'pdf'), width = 9, height = 12);
VlnPlot(iseurat, features.plot = diff.tf, nCol = 3, point.size = 0.5);
dev.off();
####
library(mclust);
iseurat@meta.data$orig.ident <- seurat.all@meta.data[rownames(iseurat@meta.data), ]$orig.ident;
adjustedRandIndex(seurat.all@meta.data[seurat.all@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]$cluster, seurat.all@meta.data[seurat.all@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]$orig.ident)                      
adjustedRandIndex(seurat.all@meta.data[!seurat.all@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]$cluster, seurat.all@meta.data[!seurat.all@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]$orig.ident)                      
