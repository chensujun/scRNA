library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(Seurat);
library(pamr);
library(plyr);
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/sig/pam');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/epitumour/objects/2019-08-27_sample13_epi.rds');
###
typeb <- readRDS('basal.cells.list.rds');
typec <- readRDS('cycle.cells.list.rds');
typel <- readRDS('luminal.cells.list.rds');
seurat.all@meta.data$type <- 'lum';
seurat.all@meta.data[typeb, ]$type <- 'basal';
seurat.all@meta.data[typec, ]$type <- 'cycle';
seurat.all <- SetIdent(seurat.all, ident.use = seurat.all@meta.data$type);
jm.out <- FindAllMarkers(seurat.all, logfc.threshold = 0);
saveRDS(jm.out, file = generate.filename('diff', 'epi_types', 'rds'));
jm.out <- readRDS('2020-08-18_diff_epi_types.rds');
mydiff <- jm.out[jm.out$p_val_adj<0.05&jm.out$avg_logFC>0.58, ];
mydata <- list();
mydata$x <- as.matrix(seurat.all@data[rownames(seurat.all@data)%in%mydiff$gene, ]);
mydata$y <- factor(seurat.all@meta.data$type);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=20)
pamr.plotcen(mydata.train, mydata, threshold=20)
#pamr.geneplot(mydata.train, mydata, threshold=20);
dev.off();

fdr.obj<- pamr.fdr(mydata.train, mydata)
pamr.plotfdr(fdr.obj)
mydata.train2 <- pamr.train(mydata,hetero="BL")
mydata.results2 <-  pamr.cv(mydata.train2, mydata)
mydata.scales <- pamr.adaptthresh(mydata.train)
mydata.train3 <- pamr.train(mydata, threshold.scale=mydata.scales)
mydata.results3 <-  pamr.cv(mydata.train3, mydata);
####
#### randome sample from lum
set.seed(100)
lum <- sample(colnames(mydata$x[, mydata$y=='lum']), 364, replace = FALSE);
mydata <- list();
mydata$x <- as.matrix(seurat.all@data[rownames(seurat.all@data)%in%mydiff$gene, c(lum, rownames(seurat.all@meta.data[seurat.all@meta.data$type%in%c('basal', 'cycle'), ]))]);
mydata$y <- factor(seurat.all@meta.data$type);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_sample', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_sample', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=1)
pamr.plotcen(mydata.train, mydata, threshold=1)
pamr.geneplot(mydata.train, mydata, threshold=1);
dev.off();
##### use basal high low group in TCGA
mydata.raw <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
basal <- read.table('../basal.in.all.cell.markers.posi.genes.txt');
mydata.scale <- data.frame(t(scale(t(mydata.raw))));
groups <- data.frame(score = colMeans(mydata.scale[rownames(mydata.scale)%in%basal$V1, ]));
groups$group <- ifelse(groups$score>median(groups$score), 'high', 'low');
mydata <- list();
mydata$x <- as.matrix(na.omit(mydata.scale[rownames(mydata.scale)%in%mydiff$gene, ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, 	);

pdf(generate.filename('type_scRNA', 'error_tcga', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=2)
pamr.plotcen(mydata.train, mydata, threshold=2)
pamr.geneplot(mydata.train, mydata, threshold=5);
dev.off();

### using raw data
mydata$x <- as.matrix(na.omit(mydata.raw[rownames(mydata.raw)%in%mydiff$gene, ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_tcga_raw', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga_raw', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=3)
pamr.plotcen(mydata.train, mydata, threshold=3)
pamr.geneplot(mydata.train, mydata, threshold=7.5);
dev.off();
###
genes <- intersect(rownames(mydata.ch), rownames(mydata.train$centroids));
groups.ch <- data.frame(cor = apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.ch$group <- ifelse(groups.ch$cor>median(groups.ch$cor), 'high', 'low');
groups.ch$group <- factor(groups.ch$group, levels = c('low', 'high'));
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
groups.ch <- na.omit(groups.ch)
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai', 'pam_tcga_raw_diffsc', 'pdf')
	);
###
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_raw_diffsc', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_raw_diffsc', 'pdf')
	);

### diff gene 
#mydiff.t <- data.frame(pval = apply(mydata.raw, 1, function(x) wilcox.test(x~groups$group)$p.value));
#mydiff.t$lfc <- apply(mydata.raw+1, 1, function(x) log2(mean(x[groups$group=='high'])/mean(x[groups$group=='low'])))
#mydiff.t$padj <- p.adjust(mydiff.t$pval);
#saveRDS(mydiff.t, file = generate.filename('diff', 'tcga_basal2', 'rds'));
mydiff.t <- readRDS('2020-08-19_diff_tcga_basal2.rds');
mydiff.tgene <- mydiff.t[mydiff.t$padj<0.05&abs(mydiff.t$lfc)>0.58, ];

mydata <- list();
mydata$x <- as.matrix(na.omit(mydata.raw[rownames(mydata.raw)%in%rownames(mydiff.tgene), ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_tcga_raw_diffhl', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga_raw_diffhl', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=8)
pamr.plotcen(mydata.train, mydata, threshold=8)
pamr.geneplot(mydata.train, mydata, threshold=10);
dev.off();
####
#genes <- pamr.listgenes(mydata.train, mydata, threshold = 8);
#genes <- intersect(rownames(mydata.ch), genes[, 'id']);
genes <- intersect(rownames(mydata.ch), rownames(mydata.train$centroids));
groups.ch <- data.frame(cor = apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.ch$group <- ifelse(groups.ch$cor>median(groups.ch$cor), 'high', 'low');
groups.ch$group <- factor(groups.ch$group, levels = c('low', 'high'));
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
groups.ch <- na.omit(groups.ch)
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai', 'pam_tcga_raw', 'pdf')
	);
###
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_raw', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_raw', 'pdf')
	);
####
mydata.eb <- readRDS('../data/EBioMed/2019-02-10_EBioMedicine2015.rds');
mydata.eb <- 2^mydata.eb;
mydata.eb.raw <- mydata.eb
myclin.eb <- readRDS('../data/EBioMed/2019-02-10_EBioMedicine2015_survival.rds');
myclin.gs <- readRDS('../data/EBioMed/2019-05-02_EBioMedicine2015_gleason.rds');
myclin.t <- readRDS('../data/EBioMed/2019-06-07_EBioMedicine2015_tstage.rds');
myclin.eb$gs <- myclin.gs[rownames(myclin.eb)];
myclin.eb$t <- myclin.t[rownames(myclin.eb)]
genes <- intersect(rownames(mydata.eb), rownames(mydata.train$centroids));
groups.eb <- data.frame(cor = apply(mydata.eb[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.eb$group <- ifelse(groups.eb$cor>median(groups.eb$cor), 'high', 'low');
groups.eb$group <- factor(groups.eb$group, levels = c('low', 'high'));
groups.eb[, 3:4] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('bcr', 'time_to_bcr')];
groups.eb <- na.omit(groups.eb)
survobj <- Surv(groups.eb$time_to_bcr, groups.eb$bcr);
create.km.plot(survobj,
	groups.eb$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed', 'pam_tcga_raw', 'pdf')
	);

groups.eb[, 5:6] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('gs', 't')];
groups.eb.hi <- groups.eb[groups.eb$t%in%c("T3", "T4"), ];
groups.eb.hi$group <- ifelse(groups.eb.hi$cor>median(groups.eb.hi$cor), 'high', 'low');
groups.eb.hi$group <- factor(groups.eb.hi$group, levels = c('low', 'high'));
groups.eb.hi <- na.omit(groups.eb.hi)
survobj <- Surv(groups.eb.hi$time_to_bcr, groups.eb.hi$bcr);
create.km.plot(survobj,
	groups.eb.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed_T34', 'pam_tcga_raw', 'pdf')
	);

####
mydata.scale <- na.omit(mydata.scale);
#mydiff.s <- data.frame(pval = apply(mydata.scale, 1, function(x) wilcox.test(x~groups$group)$p.value));
#mydiff.s$lfc <- apply(mydata.scale+1, 1, function(x) log2(mean(x[groups$group=='high'])/mean(x[groups$group=='low'])))
#mydiff.s$padj <- p.adjust(mydiff.t$pval);
#mydiff.tgene <- mydiff.t[mydiff.t$padj<0.05&mydiff.t$lfc>0.58, ];

mydata <- list();
mydata$x <- as.matrix(na.omit(mydata.scale[rownames(mydata.scale)%in%rownames(mydiff.tgene), ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_tcga_scale_diffhl', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga_scale_diffhl', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=5)
pamr.plotcen(mydata.train, mydata, threshold=5)
pamr.geneplot(mydata.train, mydata, threshold=6);
dev.off();
###
#mytest <- pamr.predict(mydata.train, mydata$x, threshold = 5);
###


#mydata.ch <- read.table('/cluster/projects/hansengroup/Public_Dataset_Hub/changhai/rnaseq/filter_genes.FPKM.xls', header = TRUE, row.names = 1);
#mydata.ch <- mydata.ch[mydata.ch$gene_name%in%rownames(mydata.raw), ];
#mydata.ch <- ddply(mydata.ch, 'gene_name', numcolwise(mean));
#rownames(mydata.ch) <- mydata.ch$gene_name;
#mydata.ch <- mydata.ch[, -1];
#saveRDS(mydata.ch, file = generate.filename('filter_genes.FPKM', 'TCGA_match', 'rds'))
#mydata.ch <- mydata.ch[, grep('^T', colnames(mydata.ch))];
#colnames(mydata.ch) <- gsub('_WTS', '', colnames(mydata.ch))
#saveRDS(mydata.ch, file = generate.filename('filter_genes.FPKM', 'TCGA_match_T', 'rds'))
myclin.ch <- read.csv('/cluster/projects/hansengroup/Public_Dataset_Hub/changhai/clinical/clinical.csv', header = TRUE);
myclin.ch$id <- paste0('T', myclin.ch$Patient.ID);
myclin.ch <- myclin.ch[, -16];
myclin.ch$time_to_bcr <- as.numeric(gsub('m', '', myclin.ch$month.to.surgery.BCR));
myclin.ch$bcr <- ifelse(myclin.ch$BCR=='YES', TRUE, FALSE);
mydata.ch.raw <- readRDS('2020-08-19_filter_genes.FPKM_TCGA_match_T.rds');
mydata.ch.scale <- data.frame(t(scale(t(mydata.ch.raw))));
mytest <- pamr.predict(mydata.train, mydata.ch.scale, threshold = 5);
genes <- intersect(rownames(mydata.ch), rownames(mydata.train$centroids));
groups.ch <- data.frame(cor = apply(mydata.ch.scale[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.ch$group <- ifelse(groups.ch$cor>0, 'high', 'low');
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
groups.ch <- na.omit(groups.ch)
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
myplot <- create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);
pdf(generate.filename('KM_changhai', 'pam_tcga', 'pdf'));
myplot;
dev.off();
####
groups.ch <- data.frame(mean = colMeans(mydata.ch.scale[rownames(mydata.ch.scale)%in%basal$V1, ]));
groups.ch$group <- ifelse(groups.ch$mean>median(groups.ch$mean), 'high', 'low');
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai', 'mean', 'pdf')
	);
####
### using signature genes and raw data
### using raw data
mydata$x <- as.matrix(na.omit(mydata.raw[rownames(mydata.raw)%in%basal$V1, ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_tcga_raw_sig', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga_raw_sig', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=6)
pamr.plotcen(mydata.train, mydata, threshold=6)
pamr.geneplot(mydata.train, mydata, threshold=7.5);
dev.off();
###
genes <- intersect(rownames(mydata.ch), rownames(mydata.train$centroids));
groups.ch <- data.frame(cor = apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.ch$group <- ifelse(groups.ch$cor>median(groups.ch$cor), 'high', 'low');
groups.ch$group <- factor(groups.ch$group, levels = c('low', 'high'));
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
groups.ch <- na.omit(groups.ch)
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai', 'pam_tcga_raw_sig', 'pdf')
	);
###
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_raw_sig', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_raw_sig', 'pdf')
	);

genes <- intersect(rownames(mydata.eb), rownames(mydata.train$centroids));
groups.eb <- data.frame(cor = apply(mydata.eb[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.eb$group <- ifelse(groups.eb$cor>median(groups.eb$cor), 'high', 'low');
groups.eb$group <- factor(groups.eb$group, levels = c('low', 'high'));
groups.eb[, 3:4] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('bcr', 'time_to_bcr')];
groups.eb <- na.omit(groups.eb)
survobj <- Surv(groups.eb$time_to_bcr, groups.eb$bcr);
create.km.plot(survobj,
	groups.eb$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed', 'pam_tcga_raw_sig', 'pdf')
	);

groups.eb[, 5:6] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('gs', 't')];
groups.eb.hi <- groups.eb[groups.eb$t%in%c("T3", "T4"), ];
groups.eb.hi$group <- ifelse(groups.eb.hi$cor>median(groups.eb.hi$cor), 'high', 'low');
groups.eb.hi$group <- factor(groups.eb.hi$group, levels = c('low', 'high'));
groups.eb.hi <- na.omit(groups.eb.hi)
survobj <- Surv(groups.eb.hi$time_to_bcr, groups.eb.hi$bcr);
create.km.plot(survobj,
	groups.eb.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed_T34', 'pam_tcga_raw_sig', 'pdf')
	);
####
### using signature genes and scaled data
mydata <- list();
mydata$x <- as.matrix(na.omit(mydata.scale[rownames(mydata.scale)%in%basal$V1, ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_tcga_scale_sig', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga_scale_sig', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=0)
pamr.plotcen(mydata.train, mydata, threshold=0)
pamr.geneplot(mydata.train, mydata, threshold=5);
dev.off();
save(mydata.train, mydata.results, file = generate.filename('scale_sig', 'basal', 'rda'));
###
mydata.ch <- mydata.ch.scale;
genes <- intersect(rownames(mydata.ch), rownames(mydata.train$centroids));
groups.ch <- data.frame(cor = apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'], method = 'spearman')));
groups.ch$group <- ifelse(groups.ch$cor>median(groups.ch$cor), 'high', 'low');
groups.ch$group <- factor(groups.ch$group, levels = c('low', 'high'));
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
#groups.ch <- na.omit(groups.ch)
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai', 'pam_tcga_scale_sig', 'pdf')
	);
###
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
mydata.ch <- data.frame(t(scale(t(mydata.ch.raw[, rownames(groups.ch.hi)]))));
groups.ch.hi$cor <- apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'], method = 'spearman'));
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_scale_sig', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
mydata.ch <- data.frame(t(scale(t(mydata.ch.raw[, rownames(groups.ch.hi)]))));
groups.ch.hi$cor <- apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'], method = 'spearman'));
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_scale_sig', 'pdf')
	);
### pamr predict
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
mydata.ch <- data.frame(t(scale(t(mydata.ch.raw[, rownames(groups.ch.hi)]))));
#groups.ch.hi$group <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0);
groups.ch.hi$value <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0, 'posterior')[, 'high']
groups.ch.hi$group <- ifelse(groups.ch.hi$value>median(groups.ch.hi$value), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_scale_sig_prd', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
mydata.ch <- data.frame(t(scale(t(mydata.ch.raw[, rownames(groups.ch.hi)]))));
#groups.ch.hi$group <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0);
groups.ch.hi$value <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0, 'posterior')[, 'high']
groups.ch.hi$group <- ifelse(groups.ch.hi$value>median(groups.ch.hi$value), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_scale_sig_prd', 'pdf')
	);
###
mydata.eb <- data.frame(t(scale(t(mydata.eb))))
genes <- intersect(rownames(mydata.eb), rownames(mydata.train$centroids));
groups.eb <- data.frame(cor = apply(mydata.eb[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.eb$group <- ifelse(groups.eb$cor>median(groups.eb$cor), 'high', 'low');
groups.eb$group <- factor(groups.eb$group, levels = c('low', 'high'));
groups.eb[, 3:4] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('bcr', 'time_to_bcr')];
groups.eb <- na.omit(groups.eb)
survobj <- Surv(groups.eb$time_to_bcr, groups.eb$bcr);
create.km.plot(survobj,
	groups.eb$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed', 'pam_tcga_scale_sig', 'pdf')
	);

groups.eb[, 5:6] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('gs', 't')];
groups.eb.hi <- groups.eb[groups.eb$t%in%c("T3", "T4"), ];
groups.eb.hi$group <- ifelse(groups.eb.hi$cor>median(groups.eb.hi$cor), 'high', 'low');
groups.eb.hi$group <- factor(groups.eb.hi$group, levels = c('low', 'high'));
groups.eb.hi <- na.omit(groups.eb.hi)
survobj <- Surv(groups.eb.hi$time_to_bcr, groups.eb.hi$bcr);
create.km.plot(survobj,
	groups.eb.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed_T34', 'pam_tcga_scale_sig', 'pdf')
	);
### using signature genes and scaled data, column scaled 
mydata.scale.c <- data.frame(scale(mydata.raw));
mydata$x <- as.matrix(na.omit(mydata.scale.c[rownames(mydata.scale.c)%in%basal$V1, ]));
mydata$y <- factor(groups$group);
mydata$genenames <- mydata$geneids <- rownames(mydata$x);
mydata$samplelabels <- colnames(mydata$x);
###
mydata.train <- pamr.train(mydata);
mydata.results<- pamr.cv(mydata.train, mydata);

pdf(generate.filename('type_scRNA', 'error_tcga_scaleC_sig', 'pdf'))
pamr.plotcv(mydata.results);
dev.off();

pdf(generate.filename('type_scRNA', 'result_tcga_scaleC_sig', 'pdf'));
pamr.plotcvprob(mydata.results, mydata, threshold=0)
pamr.plotcen(mydata.train, mydata, threshold=0)
pamr.geneplot(mydata.train, mydata, threshold=5);
dev.off();
###
mydata.ch <- data.frame(scale(mydata.ch.raw));
genes <- intersect(rownames(mydata.ch), rownames(mydata.train$centroids));
groups.ch <- data.frame(cor = apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.ch$group <- ifelse(groups.ch$cor>median(groups.ch$cor), 'high', 'low');
groups.ch$group <- factor(groups.ch$group, levels = c('low', 'high'));
groups.ch[, 3:4] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('bcr', 'time_to_bcr')];
#groups.ch <- na.omit(groups.ch)
survobj <- Surv(groups.ch$time_to_bcr, groups.ch$bcr);
create.km.plot(survobj,
	groups.ch$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai', 'pam_tcga_scaleC_sig', 'pdf')
	);
###
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
groups.ch.hi$cor <- apply(mydata.ch[genes, rownames(groups.ch.hi)], 2, function(x) cor(x, mydata.train$centroids[genes, 'high']));
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_scaleC_sig', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
groups.ch.hi$cor <- apply(mydata.ch[genes, rownames(groups.ch.hi)], 2, function(x) cor(x, mydata.train$centroids[genes, 'high']));
groups.ch.hi$group <- ifelse(groups.ch.hi$cor>median(groups.ch.hi$cor), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_scaleC_sig', 'pdf')
	);
### pamr predict
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
#groups.ch.hi$group <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0);
groups.ch.hi$value <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0, 'posterior')[, 'high']
groups.ch.hi$group <- ifelse(groups.ch.hi$value>median(groups.ch.hi$value), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_GS910', 'pam_tcga_scaleC_sig_prd', 'pdf')
	);
groups.ch.hi <- groups.ch[grep('pT3|pT4', groups.ch$Pathology.T), ];
groups.ch.hi$value <- pamr.predict(mydata.train, mydata.ch[, rownames(groups.ch.hi)], 0, 'posterior')[, 'high']
groups.ch.hi$group <- ifelse(groups.ch.hi$value>median(groups.ch.hi$value), 'high', 'low');
groups.ch.hi$group <- factor(groups.ch.hi$group, levels = c('low', 'high'));
survobj <- Surv(groups.ch.hi$time_to_bcr, groups.ch.hi$bcr);
create.km.plot(survobj,
	groups.ch.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	filename = generate.filename('KM_changhai_T34', 'pam_tcga_scaleC_sig_prd', 'pdf')
	);
###
mydata.eb <- data.frame(scale(mydata.eb.raw))
genes <- intersect(rownames(mydata.eb), rownames(mydata.train$centroids));
groups.eb <- data.frame(cor = apply(mydata.eb[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'])));
groups.eb$group <- ifelse(groups.eb$cor>median(groups.eb$cor), 'high', 'low');
groups.eb$group <- factor(groups.eb$group, levels = c('low', 'high'));
groups.eb[, 3:4] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('bcr', 'time_to_bcr')];
groups.eb <- na.omit(groups.eb)
survobj <- Surv(groups.eb$time_to_bcr, groups.eb$bcr);
create.km.plot(survobj,
	groups.eb$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed', 'pam_tcga_scaleC_sig', 'pdf')
	);

groups.eb[, 5:6] <- myclin.eb[match(rownames(groups.eb), rownames(myclin.eb)), c('gs', 't')];
groups.eb.hi <- groups.eb[groups.eb$t%in%c("T3", "T4"), ];
groups.eb.hi$group <- ifelse(groups.eb.hi$cor>median(groups.eb.hi$cor), 'high', 'low');
groups.eb.hi$group <- factor(groups.eb.hi$group, levels = c('low', 'high'));
groups.eb.hi <- na.omit(groups.eb.hi)
survobj <- Surv(groups.eb.hi$time_to_bcr, groups.eb.hi$bcr);
create.km.plot(survobj,
	groups.eb.hi$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE,
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	filename = generate.filename('KM_EBioMed_T34', 'pam_tcga_scaleC_sig', 'pdf')
	);
