library(BoutrosLab.plotting.general);
library(scran);
library(Seurat);
library(scater);

setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
sce <- readRDS(conf$scran_norm);
#####
batch <- as.numeric(as.factor(sce$batch));
#####
count1 <- logcounts(sce)[, batch == 1];
count2 <- logcounts(sce)[, batch == 2];
count3 <- logcounts(sce)[, batch == 3];
count4 <- logcounts(sce)[, batch == 4];
count5 <- logcounts(sce)[, batch == 5];
count6 <- logcounts(sce)[, batch == 6];
count7 <- logcounts(sce)[, batch == 7];
count8 <- logcounts(sce)[, batch == 8];
count9 <- logcounts(sce)[, batch == 9];
count10 <- logcounts(sce)[, batch == 10];
count11 <- logcounts(sce)[, batch == 11];
count12 <- logcounts(sce)[, batch == 12];
count13 <- logcounts(sce)[, batch == 13];
#####
original <- list(count1, count2, count3, count4, count5, count6, count7, count8,
	count9, count10, count11, count12, count13);
names(original) <- levels(as.factor(sce$batch));
save(original, batch, file = generate.filename('normalize', 'tmp', 'rda'));
####
name <- 'all_auto_k5';
set.seed(20190317);
out <- do.call(fastMNN, c(original, list(k = 5, d = 50, approximate = TRUE, auto.order = TRUE)));
#saveRDS(out, file = generate.filename(name, 'mnn', 'rds'));
saveRDS(out, file = paste0('objects/', Sys.Date(), '_', name, '_mnn.rds'));

####
combined <- do.call(cbind, original);
reducedDim(sce, 'MNN') <- out$corrected[match(colnames(sce), colnames(combined)), ];
#saveRDS(sce, file = generate.filename(name, 'sce_mnn', 'rds'));
saveRDS(sce, file = paste0('objects/', Sys.Date(), '_', name, '_mnn_sce.rds'));
#### add MNN corrected gene-level value to sce data
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
load('2019-03-18_normalize_tmp.rda');
combined <- do.call(cbind, original);
sce <- readRDS(conf$scran_mnn);
out <- readRDS('objects/2019-03-20_all_auto_k5_mnn.rds');
cor.exp <- tcrossprod(out$rotation, out$corrected);
cor.exp <- cor.exp[, match(colnames(sce), colnames(combined))];
rownames(cor.exp) <- rownames(sce);
colnames(cor.exp) <- colnames(sce);
normcounts(sce) <- cor.exp;
saveRDS(sce, file = paste0('objects/', Sys.Date(), '_', name, '_mnn_sce_norm.rds'));
saveRDS(cor.exp, file = paste0('objects/', Sys.Date(), '_', name, '_mnn_corrected_exp.rds'));
####
rownames(out$rotation) <- rownames(combined);
rownames(out$corrected) <- colnames(combined);
saveRDS(out, file = paste0('objects/', Sys.Date(), '_', name, '_mnn_names.rds'));
####
#qsub -cwd -b y -N dimred_tmp -l h_vmem=125G "module load R; Rscript run_dimred_tmp.R"