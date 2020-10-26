library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(UpSetR);
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/vCAF');
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
iseurat <- readRDS(conf$seurat_endo);
name <- 'vCAF';
seurat.all@ident <- factor(ifelse(rownames(seurat.all@meta.data)%in%rownames(iseurat@meta.data[iseurat@ident%in%c(2,3,4,5), ]), 1, 0));
names(seurat.all@ident) <- rownames(seurat.all@meta.data);
jm.out <- FindAllMarkers(seurat.all, logfc.threshold = 0);
saveRDS(jm.out, file = generate.filename('diff_all', name, 'rds'));
#####
sig.vcaf <- jm.out[jm.out$cluster==1&jm.out$p_val_adj<0.05&jm.out$avg_logFC>log(2), ]$gene;
####
#run part of ~/svn/singleCell/archive/cpcg_gene.R and save the RPKM to file
#saveRDS(rpkm.all.ln, file = generate.filename('rpkm', 'all', 'rds'));
load('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/data/2019-08-24_rpkm_all.rds');
load('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/data/2019-08-24_rpkm_all_rank.rds');
load("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/ML/2018-01-23_features_parent_gene_circRNA.rda");
to.plot <- data.frame(exp = colMeans(rpkm.all.ln[rownames(rpkm.all.ln)%in%rownames(allEns[allEns$gene%in%sig.vcaf, ]), ]),
        cohorts = rep(c('NGS','CH', 'CPC', 'NGS', 'CH', 'MET', 'NEPC', 'CRPC'),
                times = c(ncol(g.ln.n), ncol(ch.ln.n), ncol(rpkm.ln.cpc), ncol(g.ln.t), ncol(ch.ln.t), ncol(rpkm.ml), ncol(rpkm.nepc), ncol(rpkm.crpc))),
        status =  rep(c('N','N', 'P', 'P', 'P', 'M', 'NE', 'CR'),
                times = c(ncol(g.ln.n), ncol(ch.ln.n), ncol(rpkm.ln.cpc), ncol(g.ln.t), ncol(ch.ln.t), ncol(rpkm.ml), ncol(rpkm.nepc), ncol(rpkm.crpc))));
to.plot$groups <- factor(paste0(to.plot$cohorts, '_', to.plot$status), levels = c('CH_N', 'NGS_N', 'CPC_P', 'CH_P', 'NGS_P', 'MET_M', 'CRPC_CR','NEPC_NE'));
to.plot$actb <- unlist(rpkm.all.ln['ENSG00000075624.13', ]);
to.plot$ratio <- to.plot$exp/to.plot$actb*100;
create.boxplot(
        data = to.plot,
        formula = ratio ~groups,
        xlab.label = '',
        ylab.label = expression('%ACTB'),
        filename = generate.filename('sig_vcaf_actb_ratio', 'types', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(to.plot$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 1, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
####
#rpkm.rank <- apply(rpkm.all.ln, 2, function(x) rank(-x));
#saveRDS(rpkm.rank, file = generate.filename('rpkm_all', 'rank', 'rds'));
to.plot$mrank <- colMeans(rpkm.rank[rownames(rpkm.rank)%in%rownames(allEns[allEns$gene%in%sig.vcaf, ]), ]);
create.boxplot(
        data = to.plot,
        formula = mrank ~groups,
        xlab.label = '',
        ylab.label = expression('Rank mean'),
        filename = generate.filename('sig_vcaf_rank', 'types', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(to.plot$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 1, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
####
sig.vcaf <- jm.out[jm.out$cluster==1&jm.out$p_val_adj<0.05&jm.out$avg_logFC<(-log(2)), ]$gene;
to.plot$mrank_dn <- colMeans(rpkm.rank[rownames(rpkm.rank)%in%rownames(allEns[allEns$gene%in%sig.vcaf, ]), ]);
create.boxplot(
        data = to.plot,
        formula = mrank_dn ~groups,
        xlab.label = '',
        ylab.label = expression('Rank mean'),
        filename = generate.filename('sig_vcaf_rank_dn', 'types', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(to.plot$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 1, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
###
create.boxplot(
        data = to.plot,
        formula = n0 ~groups,
        xlab.label = '',
        ylab.label = expression('Rank mean'),
        filename = generate.filename('group', 'n0', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(to.plot$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 1, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
####
sig.up <- jm.out[jm.out$cluster==1&jm.out$p_val_adj<0.05&jm.out$avg_logFC>log(2), ]$gene;
sig.dn <- jm.out[jm.out$cluster==1&jm.out$p_val_adj<0.05&jm.out$avg_logFC<(-log(2)), ]$gene;
rpkm.rank.s <- rpkm.su;
rpkm.rank.s[, -1] <- apply(rpkm.rank.s[, -1], 2, function(x) rank(-x));
rpkm.su <- read.table('~/circRNA/star-circ/lrf/expression/data_RNA_Seq_expression_median.txt', header = TRUE);
dat.su <- data.frame(exp = colMeans(rpkm.su[rpkm.su$Hugo_Symbol%in%sig.up, -1]), 
	cohorts = 'SU2C', status = 'M', groups = 'SU2C_M', actb = unlist(rpkm.su[rpkm.su$Hugo_Symbol=='ACTB', -1]),
	ratio = NA, mrank = colMeans(rpkm.rank.s[rpkm.rank.s$Hugo_Symbol%in%sig.up, -1]), 
	mrank_dn = colMeans(rpkm.rank.s[rpkm.rank.s$Hugo_Symbol%in%sig.dn, -1]));
dat.su$ratio <- 100*(dat.su$exp/dat.su$actb);
plot.all <- rbind(to.plot[, -9], dat.su);
plot.all$groups <- factor(plot.all$groups, levels = c('CH_N', 'NGS_N', 'CPC_P', 'CH_P', 'NGS_P', 'MET_M', 'SU2C_M', 'CRPC_CR','NEPC_NE'));

plot.all$exp_up <- c(colMeans(rpkm.all.ln[rownames(rpkm.all.ln)%in%rownames(allEns[allEns$gene%in%sig.up, ]), ]),
	colMeans(rpkm.su[rpkm.su$Hugo_Symbol%in%sig.up, -1]));
plot.all$exp_dn <- c(colMeans(rpkm.all.ln[rownames(rpkm.all.ln)%in%rownames(allEns[allEns$gene%in%sig.dn, ]), ]),
	colMeans(rpkm.su[rpkm.su$Hugo_Symbol%in%sig.dn, -1]));
plot.all$ratio_up <- 100*(plot.all$exp_up/plot.all$actb);
plot.all$ratio_dn <- 100*(plot.all$exp_dn/plot.all$actb);

create.boxplot(
        data = plot.all,
        formula = mrank ~groups,
        xlab.label = '',
        ylab.label = expression('Rank mean'),
        filename = generate.filename('sig_vcaf_rank', 'all_su2c', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(plot.all$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 2, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );

create.boxplot(
        data = plot.all,
        formula = mrank_dn ~groups,
        xlab.label = '',
        ylab.label = expression('Rank mean'),
        filename = generate.filename('sig_vcaf_rank_dn', 'all_su2c', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(plot.all$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 2, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
###
create.boxplot(
        data = plot.all,
        formula = ratio_up ~groups,
        xlab.label = '',
        ylab.label = expression('%ACTB'),
        filename = generate.filename('sig_vcaf_up', 'all_su2c', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(plot.all$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 2, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );

create.boxplot(
        data = plot.all,
        formula = ratio_dn ~groups,
        xlab.label = '',
        ylab.label = expression('%ACTB'),
        filename = generate.filename('sig_vcaf_dn', 'all_su2c', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(plot.all$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 2, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
### reference signature
dat.ref <- rowMeans(iseurat@data[jm.out[jm.out$cluster==0&jm.out$p_val_adj<0.05&abs(jm.out$avg_logFC)>log(2), ]$gene, 
	iseurat@ident%in%c(2,3,4,5)]);
dat.ref <- dat.ref[names(dat.ref)%in%rpkm.su$Hugo_Symbol];
dat.ref.ens <- dat.ref;
names(dat.ref.ens) <- rownames(allEns[match(names(dat.ref.ens), allEns$gene), ]);
dat.ref.ens <- dat.ref.ens[names(dat.ref.ens)%in%rownames(rpkm.all.ln)];
plot.all$cor <- c(
	apply(rpkm.all.ln[names(dat.ref.ens), ], 2, function(x) cor(x, dat.ref.ens)),
	apply(rpkm.su[match(names(dat.ref), rpkm.su$Hugo_Symbol), -1], 2, function(x) cor(x, dat.ref))	
	);
create.boxplot(
        data = plot.all,
        formula = cor ~groups,
        xlab.label = '',
        ylab.label = expression('Correlation'),
        filename = generate.filename('sig_vcaf_cor', 'all_su2c', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(plot.all$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 2, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );
#### use only signature from cluster 3,4,5
seurat.all@ident <- factor(ifelse(rownames(seurat.all@meta.data)%in%rownames(iseurat@meta.data[iseurat@ident%in%c(3,4,5), ]), 1, 0));
names(seurat.all@ident) <- rownames(seurat.all@meta.data);
jm.out <- FindMarkers(seurat.all, logfc.threshold = 0, ident.1 = 1);
saveRDS(jm.out, file = generate.filename('diff_allv345', name, 'rds'));
#####
dat.ref <- rowMeans(iseurat@data[rownames(jm.out[jm.out$p_val_adj<0.05&abs(jm.out$avg_logFC)>log(2), ]), 
	iseurat@ident%in%c(3,4,5)]);
dat.ref <- dat.ref[names(dat.ref)%in%rpkm.su$Hugo_Symbol];
dat.ref.ens <- dat.ref;
names(dat.ref.ens) <- rownames(allEns[match(names(dat.ref.ens), allEns$gene), ]);
dat.ref.ens <- dat.ref.ens[names(dat.ref.ens)%in%rownames(rpkm.all.ln)];
plot.all$cor345 <- c(
	apply(rpkm.all.ln[names(dat.ref.ens), ], 2, function(x) cor(x, dat.ref.ens)),
	apply(rpkm.su[match(names(dat.ref), rpkm.su$Hugo_Symbol), -1], 2, function(x) cor(x, dat.ref))	
	);
create.boxplot(
        data = plot.all,
        formula = cor345 ~groups,
        xlab.label = '',
        ylab.label = expression('Correlation'),
        filename = generate.filename('sig_vcaf_cor345', 'all_su2c', 'pdf'),
        xaxis.lab = gsub('_N|_P|_M|_NE|_CR', '', levels(plot.all$group)),
        add.stripplot = TRUE,
        lwd = 1,
        style = 'Nature',
        border.col = rep(default.colours(5), times = c(2, 3, 2, 1, 1)),
        height = 4,
        width = 8,
        key = list(
                lines = list(
                        lty = 1,
                        col = default.colours(5)
                        ),
                text = list(
                        col = default.colours(5),
                        lab = c('Normal', 'Primary', 'Metastatic', 'NEPC', 'CRPC')
                        ),
                x = 0.01,
                y = 0.99
                ),
        abline.h = 5,
        abline.lty = 2
        );