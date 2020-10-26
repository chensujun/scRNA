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
saveRDS(genes.marker, file = generate.filename('group', 'NE', 'rds'));
genes.marker$type <- seurat.all@meta.data[rownames(genes.marker), ]$type;
###
iseurat@meta.data$group <- ifelse(genes.marker[rownames(iseurat@meta.data), ]$group=='NE', 'NE', 'others');
pdf(generate.filename(name, 'group_NE', 'pdf'), width = 8);
DimPlot(iseurat, reduction.use = 'TSNE', group.by = 'group', plot.order = 'NE', cols.use = c('grey', 'blue'));
dev.off();
### difference between NE and Lum, pulled data
gene.go <- readRDS('~/circRNA/star-circ/circRNA_landscape/rnaseq_landscape/scRNA/2019-01-05_gp_hsapiens.GO.Name.rds');
to.plot <- as.data.frame.matrix(table(genes.marker$patient,genes.marker$group));
to.plot$ratio <- 100*to.plot$NE/to.plot$Lum;
iseurat.1 <- SubsetData(iseurat, cells.use = rownames(genes.marker[genes.marker$group%in%c('Lum', 'NE'), ]));
###
iseurat.1 <- SetIdent(iseurat.1, ident.use = genes.marker[iseurat.1@cell.names, ]$group);
jm.out <- FindAllMarkers(iseurat.1);
jm.out <- jm.out[jm.out$cluster == 'Lum', ];
jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
ego.up <- enrich_gprofiler(jup, rownames(seurat.all@data));
ego.dn <- enrich_gprofiler(jdn, rownames(seurat.all@data));
red.ego.up <- run_simplify(ego.up);
red.ego.dn <- run_simplify(ego.dn);
save(jm.out, ego.up, ego.dn, red.ego.up, red.ego.dn, file = generate.filename('diff_NE_Lum_pull', j, 'rda'));

###
jseurat <- SubsetData(iseurat, cells.use = rownames(genes.marker[genes.marker$group%in%c('NE'), ]));
### difference between NE and Lum, per patient
for(j in rownames(to.plot[to.plot$NE>10&to.plot$Lum>10, ])){
	jseurat <- SubsetData(iseurat.1, cells.use = rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==j, ]));
	jseurat <- SetIdent(jseurat, ident.use = genes.marker[jseurat@cell.names, ]$group);
	jm.out <- FindAllMarkers(jseurat);
	jm.out <- jm.out[jm.out$cluster == 'Lum', ];
	jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
	ego.up <- enrich_gprofiler(jup, rownames(seurat.all@data));
	ego.dn <- enrich_gprofiler(jdn, rownames(seurat.all@data));
	save(jm.out, ego.up, ego.dn, file = generate.filename('diff_NE_Lum', j, 'rda'))
	#####
};

mydiff.up <- mydiff.dn <- mydiff <- myego.up <- myego.dn <- {};
for(j in rownames(to.plot[to.plot$NE>10&to.plot$Lum>10, ])){
	load(generate.filename('diff_NE_Lum', j, 'rda'));
	red.ego.up <- run_simplify(ego.up);
	red.ego.dn <- run_simplify(ego.dn);
	myego.up[[j]] <- red.ego.up;	
	myego.dn[[j]] <- red.ego.dn;
	mydiff[[j]] <- jm.out
	mydiff.up[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	mydiff.dn[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
};
#####
plot_enrich(myego.up, gene.go, 'BP', 'combine_up', 10, 15, mybg = rownames(iseurat@data), xrot = 90);
plot_enrich(myego.dn, gene.go, 'BP', 'combine_dn', 10, 15, mybg = rownames(iseurat@data), xrot = 90);
#### plot keg
mydomain <- 'keg'
mydiff.up <- mydiff.dn <- mydiff <- myego.up <- myego.dn <- {};
for(j in rownames(to.plot[to.plot$NE>10&to.plot$Lum>10, ])){
	load(paste0('2019-06-07_diff_NE_Lum_', j, '.rda'));
	red.ego.up <- ego.up[ego.up$domain==mydomain, ];
	red.ego.dn <- ego.dn[ego.dn$domain==mydomain, ];
	myego.up[[j]] <- red.ego.up;	
	myego.dn[[j]] <- red.ego.dn;
	mydiff[[j]] <- jm.out
	mydiff.up[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	mydiff.dn[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
};
####
plot_enrich(myego.up, gene.go, 'keg', 'combine_up_NE_Lum', 10, 8, mybg = rownames(iseurat@data), xrot = 90, noOR = TRUE, spot.size.function = function(x) {abs(x)*2});
plot_enrich(myego.dn, gene.go, 'keg', 'combine_dn_NE_Lum', 10, 8, mybg = rownames(iseurat@data), xrot = 90, noOR = TRUE, spot.size.function = function(x) {abs(x)*2});

#### compare AR+/-
for(j in rownames(to.plot)){
	jseurat <- SubsetData(iseurat, cells.use = rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==j, ]));
	jseurat <- SetIdent(jseurat, ident.use = ifelse(jseurat@data['AR', ]>0, 1, 0));
	jm.out <- FindAllMarkers(jseurat);
	jm.out <- jm.out[jm.out$cluster == '1', ];
	jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
	ego.up <- enrich_gprofiler(jup, rownames(seurat.all@data));
	ego.dn <- enrich_gprofiler(jdn, rownames(seurat.all@data));
	save(jm.out, ego.up, ego.dn, file = generate.filename('diff_AR', j, 'rda'))
	#####
};


mydiff.up <- mydiff.dn <- mydiff <- myego.up <- myego.dn <- {};
mydomain <- 'BP'
for(j in rownames(to.plot)){
	load(paste0('2019-06-07_diff_AR_', j, '.rda'));
	red.ego.up <- run_simplify(ego.up);
	red.ego.dn <- run_simplify(ego.dn);
	#rownames(jm.out) <- gsub('\\..*', '', rownames(jm.out))
	#jup <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	#jdn <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
	#ego.up <- enrich_gprofiler(jup, rownames(seurat.all@data));
	#ego.dn <- enrich_gprofiler(jdn, rownames(seurat.all@data));
	#save(jm.out, ego.up, ego.dn, file = generate.filename('diff_AR', j, 'rda'))

	myego.up[[j]] <- red.ego.up;	
	myego.dn[[j]] <- red.ego.dn;
	mydiff[[j]] <- jm.out
	mydiff.up[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	mydiff.dn[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
};
#####
plot_enrich(myego.up, gene.go, 'BP', 'combine_up_AR', 10, 15, mybg = rownames(iseurat@data), xrot = 90);
plot_enrich(myego.dn, gene.go, 'BP', 'combine_dn_AR', 10, 15, mybg = rownames(iseurat@data), xrot = 90);
####
mydomain <- 'keg'
mydiff.up <- mydiff.dn <- mydiff <- myego.up <- myego.dn <- {};
for(j in rownames(to.plot)){
	load(paste0('2019-06-07_diff_AR_', j, '.rda'));
	red.ego.up <- ego.up[ego.up$domain==mydomain, ];
	red.ego.dn <- ego.dn[ego.dn$domain==mydomain, ];
	myego.up[[j]] <- red.ego.up;	
	myego.dn[[j]] <- red.ego.dn;
	mydiff[[j]] <- jm.out
	mydiff.up[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0, ]);
	mydiff.dn[[j]] <- rownames(jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC<0, ]);
};
####
plot_enrich(myego.up, gene.go, 'keg', 'combine_up_AR', 10, 8, mybg = rownames(iseurat@data), xrot = 90, noOR = TRUE, spot.size.function = function(x) {abs(x)*2});
plot_enrich(myego.dn, gene.go, 'keg', 'combine_dn_AR', 10, 8, mybg = rownames(iseurat@data), xrot = 90, noOR = TRUE, spot.size.function = function(x) {abs(x)*2});

####
mykeg.up <- cbind(myego.up[[1]], group=names(myego.up)[1]);
for(i in names(myego.up)){
	if(nrow(myego.up[[i]])>0){
	mykeg.up <- rbind(mykeg.up, cbind(myego.up[[i]], group = i));
}
};
to.plot <- mykeg.up[, -14]
create.barplot()