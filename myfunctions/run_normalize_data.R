library(BoutrosLab.plotting.general);
library(scran);
library(Seurat);
library(scater);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
####
args <- commandArgs(trailingOnly = TRUE);
mydat <- args[1];
name <- args[2];
####
epi.ne <- readRDS(mydat);
#name <- 'NEDN'
count.orig <- as.matrix(epi.ne@raw.data[, colnames(epi.ne@data)]);
print(dim(count.orig));
mycount <- SingleCellExperiment(list(counts = count.orig));
#clusters <- quickCluster(mycount, min.size = 100);
mycount <- computeSumFactors(mycount);
mycount <- normalize(mycount);
###
saveRDS(mycount, file = paste0(Sys.Date(), '_', name, '_normalized.rds'));

### block batch effect
batch <- as.numeric(as.factor(epi.ne@meta.data$orig.ident));
#fit <- trendVar(mycount, parametric = TRUE, use.spikes = FALSE, block = batch);
#decomp <- decomposeVar(mycount, fit);
#top.hvgs <- order(decomp$bio, decreasing = TRUE);
######
#pdf(generate.filename('decompose', paste0('variance_', name), 'pdf'));
#plot(decomp$mean, decomp$total, pch = 16, cex = 0.6, xlab="Mean log-expression", 
#    ylab="Variance of log-expression");
#curve(fit$trend(x), col="dodgerblue", lwd=2, add=TRUE);
#dev.off();	

#mycol <- rainbow(2*length(batch))[seq(1, 2*length(batch), 2)]#

#pdf(generate.filename('decompose', paste0('variance_patient', name), 'pdf'));
#matplot(fit$mean, fit$vars, col = mycol, 
#	xlab="Mean log-expression", ylab="Variance of log-expression");
#curve(fit$trend(x), add=TRUE, col="red");
#dev.off();
#####remove batch effect using MNN method
count1 <- logcounts(mycount)[, batch == 1];
count2 <- logcounts(mycount)[, batch == 2];
count3 <- logcounts(mycount)[, batch == 3];
count4 <- logcounts(mycount)[, batch == 4];
count5 <- logcounts(mycount)[, batch == 5];
count6 <- logcounts(mycount)[, batch == 6];
count7 <- logcounts(mycount)[, batch == 7];
count8 <- logcounts(mycount)[, batch == 8];
count9 <- logcounts(mycount)[, batch == 9];
count10 <- logcounts(mycount)[, batch == 10];
count11 <- logcounts(mycount)[, batch == 11];
count12 <- logcounts(mycount)[, batch == 12];
count13 <- logcounts(mycount)[, batch == 13, drop = FALSE];
####
original <- list(count1, count2, count3, count4, count5, count6, count7, count8,
	count9, count10, count11, count12, count13);
names(original) <- levels(as.factor(epi.ne@meta.data$orig.ident));
####run fastMNN
out <- do.call(fastMNN, c(original, list(k = 20, d = 50, approximate = TRUE)));
saveRDS(out, file = generate.filename(name, 'mnn', 'rds'));
####
combined <- do.call(cbind, original);
sce <- SingleCellExperiment(list(logcounts = combined));
reducedDim(sce, 'MNN') <- out$corrected;
#sce$Batch <- paste0('p', batch[order(batch)]);
sce$Batch <- epi.ne@meta.data[rownames(sce@colData), ]$orig.ident;

####run pca
set.seed(100);
osce <- runPCA(sce, ntop=Inf, method="irlba");
osce <- runTSNE(osce, use_dimred="PCA")
set.seed(100)
csce <- runTSNE(sce, use_dimred="MNN");
csce$group <- epi.ne@meta.data[rownames(csce@colData), ]$type;
saveRDS(sce, file = generate.filename(name, 'sce', 'rds'));
saveRDS(osce, file = generate.filename(name, 'tsne_orig', 'rds'));
saveRDS(csce, file = generate.filename(name, 'tsne_mnn', 'rds'));
####
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original") 
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
pdf(generate.filename('plottsne', paste0('batch_compare_', name), 'pdf'), width = 15)
multiplot(ot, ct, cols=2);
dev.off();
###
csce$group <- epi.ne@meta.data[rownames(csce@colData), ]$celltype;
pdf(generate.filename('plottsne', paste0('corrected_group_', name), 'pdf'), width = 8)
plotTSNE(csce, colour_by="group") + ggtitle("Corrected")
dev.off();

pdf(generate.filename('plottsne', paste0('seurat_group_', name), 'pdf'), width = 15)
TSNEPlot(epi.ne, group.by = 'fig.celltype', do.label = TRUE);
dev.off();
pdf(generate.filename('plottsne', paste0('seurat_patient_', name), 'pdf'), width = 15)
TSNEPlot(epi.ne, group.by = 'orig.ident', do.label = TRUE);
dev.off();