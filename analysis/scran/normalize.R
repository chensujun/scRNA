library(BoutrosLab.plotting.general);
library(scran);
library(Seurat);
library(scater);

setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
####
my.cells <- readRDS(conf$filtered_cells);
epi.ne <- readRDS(conf$seurat_all);
name <- 'all'
count.orig <- as.matrix(epi.ne@raw.data[, my.cells]);
sce <- SingleCellExperiment(list(counts = count.orig));
sce$batch <- epi.ne@meta.data[colnames(sce), ]$orig.ident;
###

#pre-clustering to split the cells into sensible clusters, to avoid violating the non-DE assumption and distorting the size factors.
####
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
#### Calculating size factors
sce <- computeSumFactors(sce, cluster=cluster.id, min.mean=0.1);
sce <- normalize(sce);
saveRDS(sce, file = paste0('objects/', generate.filename(name, 'clustered_normalized', 'rds')));
####
pdf(generate.filename('library', 'size_factor', 'pdf'))
plot(sce$total_counts, sizeFactors(sce), log="xy", xlab="Library size", 
    ylab="Size factors", cex=0.2, pch=16)
dev.off();
