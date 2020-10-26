library(BoutrosLab.plotting.general);
library(scran);
library(scater);
library(gdata);
#### in running this script, you might need to restart R for a couple of times
#### I feel lazy so just write combine all actions related to preprocessing, cluster and annotation together in this one script 
#### unfortunately this makes the script not qsub-able -- unless you've got big memories! 
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
#### create scran object
count.orig <- read.table('raw_data/AllSample.counts.tsv', header = TRUE, row.names = 1);
sce <- SingleCellExperiment(list(counts = as.matrix(count.orig)));
sce$batch <- gsub('_.*', '', colnames(count.orig));
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=grep('^MT-', rownames(sce), value = TRUE)));
name <- 'all11';


count.orig <- as.matrix(epi.ne@raw.data[, my.cells]);
sce <- SingleCellExperiment(list(counts = count.orig));
sce$batch <- epi.ne@meta.data[colnames(sce), ]$orig.ident;
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(rowData(sce)$Chr=="chrM")));
###
pdf(generate.filename('preprocess', 'qc', 'pdf'), height = 4, width = 6)
par(mfrow=c(1,2));
hist(sce$log10_total_counts, xlab=expression("Log"[10]~"library size"), col="grey80", main = 'Total counts')
hist(sce$log10_total_features_by_counts, xlab=expression("Log"[10]~"number of genes expressed"), 
	col="grey80", main = 'Number of genes');
dev.off();
#### Discard low-quality cells from the object, we skip this as this was already done with Seurat
#low.libsize <- isOutlier(sce$log10_total_counts, nmad=3, batch=sce$batch, type="lower");
#low.ngenes <- isOutlier(sce$log10_total_features_by_counts, nmad=3, batch=sce$batch, type="lower");
#discard <- low.libsize | low.ngenes;
###
ave.count <- calcAverage(sce)
rowData(sce)$AverageCount <- ave.count
n.cells <- nexprs(sce, byrow=TRUE)
rowData(sce)$NumCells <- n.cells
pdf(generate.filename('preprocess', 'qc_gene', 'pdf'), height = 4, width = 6)
par(mfrow=c(1,2));
hist(log10(ave.count), xlab=expression("Log"[10]~"average count"), col="grey80")
hist(log10(n.cells), xlab=expression("Log"[10]~"number of cells"), col="grey80")
dev.off();
### filter using house keeping gene
house.keeping <- read.xls('housekeeping_tirosh_aad0501_Table_S16.xlsx');
colnames(house.keeping) <- 'genes';
my.stat <- data.frame(ngene.all = colSums(count.orig>0), ngene.hk = colSums(count.orig[intersect(rownames(count.orig), house.keeping$genes), ]>0));
my.stat$pass.qc <- ifelse(rownames(my.stat)%in%rownames(sce@colData[sce@colData$pct_counts_Mito<20, ]), 'Y', 'N');
my.stat$exp.hk <- colMeans(count.orig[intersect(rownames(count.orig), house.keeping$genes), ]);
####
my.sd <- sd(my.stat[my.stat$ngene.all>2000, ]$exp.hk);
my.mean <- mean(my.stat[my.stat$ngene.all>2000, ]$exp.hk);

my.cof <- length(intersect(house.keeping$genes, rownames(sce)))*0.2
pdf(generate.filename('distribution', 'ngene_hk', 'pdf'));
plot(density(my.stat$ngene.all), main = 'nGene');
plot(density(my.stat$ngene.hk), main = '# Housekeeping genes');
abline(v = my.cof, lty = 2)
plot(density(log2(my.stat$exp.hk)), main = 'Mean housekeeping genes');
abline(v = log2(my.mean-my.sd), lty = 2)
dev.off();
###
my.sd <- sd(my.stat[my.stat$ngene.all>2000, ]$exp.hk);
my.mean <- mean(my.stat[my.stat$ngene.all>2000, ]$exp.hk);
####
####
my.cells <- rownames(my.stat[my.stat$ngene.hk>my.cof&my.stat$pass.qc=='Y'&my.stat$ngene.all>200, ]);
sce <- sce[, colnames(sce)%in%my.cells];
by.lib <- split(seq_len(ncol(sce)), sce$batch);
cluster.id <- character(ncol(sce));
for (lib in names(by.lib)) { 
	print(lib)
    current <- by.lib[[lib]]
    cur.exprs <- realize(counts(sce)[,current]) # for speed; avoid multiple file reads here.
    ids <- quickCluster(cur.exprs, min.mean=0.1, method="igraph", 
        max.size=3000, irlba.args=list(maxit=1000))
    cluster.id[current] <- paste0(lib, ".", ids)
}
saveRDS(cluster.id, file = paste0('objects/', generate.filename(name, 'quick_cluster_patient', 'rds')));
#####
sce <- computeSumFactors(sce, cluster=cluster.id, min.mean=0.1);
sce <- normalize(sce);
saveRDS(sce, file = paste0('objects/', generate.filename(name, 'clustered_normalized', 'rds')));
####
pdf(generate.filename('library', 'size_factor', 'pdf'))
plot(sce$total_counts, sizeFactors(sce), log="xy", xlab="Library size", 
    ylab="Size factors", cex=0.2, pch=16)
dev.off();
##### batch removal and dim reduction with MNN
batch <- as.numeric(as.factor(sce$batch));
#####
count.log <- logcounts(sce);
count1 <- count.log[, batch == 1];
count2 <- count.log[, batch == 2];
count3 <- count.log[, batch == 3];
count4 <- count.log[, batch == 4];
count5 <- count.log[, batch == 5];
count6 <- count.log[, batch == 6];
count7 <- count.log[, batch == 7];
count8 <- count.log[, batch == 8];
count9 <- count.log[, batch == 9];
count10 <- count.log[, batch == 10];
count11 <- count.log[, batch == 11];

original <- list(count1, count2, count3, count4, count5, count6, count7, count8,
	count9, count10, count11);
names(original) <- levels(as.factor(sce$batch));
save(original, batch, file = generate.filename('normalize', 'tmp', 'rda'));
#### memory may not allow for running through this, try restart R and load the saved rda from above
load('objects/2019-06-25_normalize_tmp.rda');
name <- 'all11_auto_k5';
set.seed(20190317);
out <- do.call(fastMNN, c(original, list(k = 5, d = 50, approximate = TRUE, auto.order = TRUE)));
saveRDS(out, file = paste0('objects/', Sys.Date(), '_', name, '_mnn.rds'));
combined <- do.call(cbind, original);
reducedDim(sce, 'MNN') <- out$corrected[match(colnames(sce), colnames(combined)), ];
saveRDS(sce, file = paste0('objects/', Sys.Date(), '_', name, '_mnn_sce.rds'));
rownames(out$rotation) <- rownames(combined);
rownames(out$corrected) <- colnames(combined);
saveRDS(out, file = paste0('objects/', Sys.Date(), '_', name, '_mnn_names.rds'));
##### run TSNE and find cluster
csce <- run_snn(sce, 'all');
csce <- runTSNE(sce, use_dimred="MNN", method = 'irlba', max_iter = 2000)
save.session.profile(generate.filename('session', 'find_cluster_snn', 'txt'));
##### annotate cell type
