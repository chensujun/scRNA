library(BoutrosLab.plotting.general);
library(Seurat);
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/seurat/");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
source('~/svn/singleCell/myfunctions/enrich_gprofiler.R');
source('~/svn/singleCell/myfunctions/plot_enrich.R');
source('~/svn/singleCell/myfunctions/go_simplify.R');

seurat.all <- readRDS(conf$sseurat_all);
colData.epi <- readRDS('../convert/epiTumor.basal.meta.data.MNN50.seuratClustered.rds');
name <- 'epi'
mycount <- seurat.all@data;
genes.marker <- readRDS(conf$Epi_subtype);
genes.ne <- colnames(genes.marker)[1:7]
#iseurat <- SubsetData(seurat.all, cells.use = rownames(seurat.all@meta.data[!seurat.all@meta.data$type%in%c('Unknown'), ]))
iseurat <- SubsetData(seurat.all, cells.use = rownames(colData.epi));

genes.marker <- data.frame(t(mycount[c(genes.ne, "KLK3", "AR", "ACPP", "KRT8"), rownames(iseurat@meta.data)]));
genes.marker$patient <- iseurat@meta.data[rownames(genes.marker), ]$orig.ident;
genes.marker$group <- 'Others';
genes.marker[rowSums(genes.marker[, genes.ne])>0&genes.marker$AR==0, ]$group <- 'NE';
genes.marker[rowSums(genes.marker[, genes.ne])==0&genes.marker$AR==0, ]$group <- 'DN';
genes.marker[rowSums(genes.marker[, genes.ne])==0&genes.marker$KLK3>3.7&genes.marker$AR>0, ]$group <- 'Lum';
####
genes.marker$group_ar <- 'Others';
genes.marker[rowSums(genes.marker[, c('KLK3', 'AR')])==0, ]$group_ar <- 'neg'
genes.marker[genes.marker$KLK3>quantile(genes.marker$KLK3, probs = 0.5)&genes.marker$AR>quantile(genes.marker$AR, probs = 0.5), ]$group_ar <- 'pos';
####
jseurat <- SubsetData(iseurat, cells.use = rownames(genes.marker[genes.marker$group_ar!='Others', ]));
jseurat <- SetIdent(jseurat, ident.use = genes.marker[genes.marker$group_ar!='Others', ]$group_ar);
jm.out <- FindAllMarkers(jseurat);
jm.out <- jm.out[jm.out$cluster == 'neg', ];
jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
ego.up <- enrich_gprofiler(jup, rownames(seurat.all@data));
ego.dn <- enrich_gprofiler(jdn, rownames(seurat.all@data));
red.ego.up <- run_simplify(ego.up);
red.ego.dn <- run_simplify(ego.dn);
save(jm.out, ego.up, ego.dn, red.ego.up, red.ego.dn, file = generate.filename('diff_AR_pull', name, 'rda'));
plot_enrich(list(up = red.ego.up, dn = red.ego.dn), gene.go, 'BP', 'combine_AR_pull', 9, 10, mybg = rownames(iseurat@data), xrot = 0);
plot_enrich(list(up = ego.up[ego.up$domain=='keg', ], dn = ego.dn[ego.dn$domain=='keg', ]), gene.go, 'keg', 'combine_AR_pull', 6, 8, mybg = rownames(iseurat@data), xrot = 0, noOR = TRUE, spot.size.function = function(x) {abs(x)*8});
####
jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>1, ]);
jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<(-1), ]);
ego.up <- enrich_gprofiler(jup, rownames(seurat.all@data));
ego.dn <- enrich_gprofiler(jdn, rownames(seurat.all@data));
red.ego.up <- run_simplify(ego.up);
red.ego.dn <- run_simplify(ego.dn);
save(jm.out, ego.up, ego.dn, red.ego.up, red.ego.dn, file = generate.filename('diff_AR_pull_fc2', name, 'rda'));
plot_enrich(list(up = red.ego.up, dn = red.ego.dn), gene.go, 'BP', 'combine_AR_pull_fc2', 9, 10, mybg = rownames(iseurat@data), xrot = 0);
plot_enrich(list(up = ego.up[ego.up$domain=='keg', ], dn = ego.dn[ego.dn$domain=='keg', ]), gene.go, 'keg', 'combine_AR_pull_fc2', 6, 8, mybg = rownames(iseurat@data), xrot = 0, noOR = TRUE, spot.size.function = function(x) {abs(x)*8});
###