library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(qusage);
library(viridis);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/pub');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$seurat_tcell_bcc);
seurat.all@meta.data$type <-  gsub('.*pre|\\.|_.*|.*post', '', rownames(seurat.all@meta.data));
seurat.all <- SetIdent(seurat.all, ident.use = seurat.all@meta.data$time)
iseurat <- SubsetData(seurat.all, cells.use = rownames(seurat.all@meta.data[grepl('tcell', seurat.all@meta.data$type), ]))
name <- 'tcell_bcc';
#jm.out <- readRDS('2019-08-29_diff_pub_tcell.rds');
jm.out <- FindMarkers(iseurat, ident.1 = 'pre', ident.2 = 'post', logfc.threshold = 0);
saveRDS(jm.out, file = generate.filename('diff_t', name, 'rds'));
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
####
keg.up <- test_enrich(annot.keg, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]), mybg = NULL);
keg.dn <- test_enrich(annot.keg, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]), mybg = NULL);

ego.up <- test_enrich(annot.go, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]), mybg = NULL);
ego.dn <- test_enrich(annot.go, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]), mybg = NULL);
save(keg.up, keg.dn, ego.up, ego.dn, file = generate.filename('diff_enrich', name, 'rda'));
####
plot.gene <- grep('MLANA|SOX10|S100|PMEL', rownames(seurat.all@data), value = T)
pdf(generate.filename('marker_gene', name, 'pdf'), height = 12);
VlnPlot(seurat.all, features.plot = plot.gene, nCol = 3, point.size.use = 0);
dev.off();

plot.gene <- grep('^CEA', rownames(seurat.all@data), value = T)
pdf(generate.filename('marker_gene', name, 'pdf'), height = 3, width = 8);
VlnPlot(seurat.all, features.plot = plot.gene, nCol = 3, point.size.use = 0);
dev.off();
aa <- data.frame(cluster = seurat.all@meta.data$cluster, gzma = seurat.all@data['GZMA', ], time = seurat.all@ident)
aa <- ddply(aa, 'cluster', numcolwise(mean));
aa <- cbind(aa, as.data.frame.matrix(table(seurat.all@meta.data$cluster, seurat.all@ident)));

seurat.all <- SetIdent(seurat.all, ident.use = seurat.all@meta.data$cluster);
jm.out <- FindMarkers(seurat.all, ident.1 = '16', ident.2 = '10', logfc.threshold = 0);

keg.up <- test_enrich(annot.keg, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]), mybg = NULL);
keg.dn <- test_enrich(annot.keg, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]), mybg = NULL);

ego.up <- test_enrich(annot.go, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]), mybg = NULL);
ego.dn <- test_enrich(annot.go, rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]), mybg = NULL);
save(keg.up, keg.dn, ego.up, ego.dn, file = generate.filename('diff_enrich_cluster16v10', name, 'rda'));
####
gs <- list();
gs[['exosome']] <- annot.go[annot.go$name=='extracellular exosome', ]$gene;
flist <- list.files('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pathway_from_brca', full.name = TRUE);
for(i in flist){
        ifile <- read.table(i, header = TRUE);
        gs[[gsub('.*/|\\..*', '', i)]] <- ifile[, 1]
};
myresult <- run_qusage(seurat.all, paste0(name, '_exosome_curated_brca'), gs);
