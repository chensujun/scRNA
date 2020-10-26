library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(Seurat);
library(pamr);
library(plyr);
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/sig/pam');
### using signature genes and scaled data
mydata.raw <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
colnames(mydata.raw) <- gsub('.[0-9]+$','', colnames(mydata.raw))
myclin <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2019-01-05_tcga_survival.rds');
myclin$time_to_bcr <- myclin$days_bcr/365;
rownames(myclin) <- gsub('.01', '', rownames(myclin))
basal <- read.table('../basal.in.all.cell.markers.posi.genes.txt');
mydata.scale <- data.frame(t(scale(t(mydata.raw))));
groups <- data.frame(score = colMeans(mydata.scale[rownames(mydata.scale)%in%basal$V1, ]));
groups$group <- ifelse(groups$score>median(groups$score), 'high', 'low');
groups$group <- factor(groups$group, levels = c('low', 'high'));
survobj <- Surv(groups$time_to_bcr, groups$bcr);
create.km.plot(survobj,
  groups$group,
  xlab.label = 'Time (Years)',
  show.risktable = TRUE,
  filename = generate.filename('KM_tcga', 'mean_basal', 'pdf')
  );

####
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

genes <- intersect(rownames(mydata.scale), rownames(mydata.train$centroids));
groups$cor <- apply(mydata.scale[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'], method = 'spearman'));
groups$cor.ps <- apply(mydata.scale[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high'], method = 'pearson'));
groups$value <- pamr.predict(mydata.train, mydata.scale, 0, 'posterior')[, 'high'];
groups$prd <- pamr.predict(mydata.train, mydata.scale, 0);
groups$group <- ifelse(groups$cor.ps>median(groups$cor.ps), 'high', 'low');
groups$group <- factor(groups$group, levels = c('low', 'high'))
groups[, 6:7] <- myclin[match(rownames(groups), rownames(myclin)), c('bcr', 'time_to_bcr')];
survobj <- Surv(groups$time_to_bcr, groups$bcr);
create.km.plot(survobj,
  groups$group,
  xlab.label = 'Time (Years)',
  show.risktable = TRUE,
  filename = generate.filename('KM_tcga', 'pam_tcga_scale_sig', 'pdf')
  );

clin.t <- data.frame(readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA.prad.clinical.rds'));
rownames(clin.t) <- gsub('-', '.', rownames(clin.t));
rownames(clin.t) <- gsub('.01', '', rownames(clin.t));
groups[, 8:9] <- clin.t[match(rownames(groups), rownames(clin.t)), c('gleasonscore', 'pathologyTstage')];
groups.hi <- groups[groups$gleasonscore%in%c(9, 10), ];
mydata <- data.frame(t(scale(t(mydata.raw[, rownames(groups.hi)]))));
groups.hi$cor <- apply(mydata[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high']));
groups.hi$group <- ifelse(groups.hi$cor>median(groups.hi$cor), 'high', 'low');
groups.hi$group <- factor(groups.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.hi$time_to_bcr, groups.hi$bcr);
create.km.plot(survobj,
  groups.hi$group,
  xlab.label = 'Time (Years)',
  show.risktable = TRUE,
  filename = generate.filename('KM_tcga_GS910', 'pam_tcga_scale_sig', 'pdf')
  );
####
groups.hi$value <- pamr.predict(mydata.train, mydata[, rownames(groups.hi)], 0, 'posterior')[, 'high']
groups.hi$group <- ifelse(groups.hi$value>median(groups.hi$value), 'high', 'low');
groups.hi$group <- factor(groups.hi$group, levels = c('low', 'high'));

survobj <- Surv(groups.hi$time_to_bcr, groups.hi$bcr);
create.km.plot(survobj,
  groups.hi$group,
  xlab.label = 'Time (Years)',
  show.risktable = TRUE,
  filename = generate.filename('KM_tcga_GS910', 'pam_tcga_scale_sig_prd', 'pdf')
  );
###
groups.hi$group <- pamr.predict(mydata.train, mydata[, rownames(groups.hi)], 0);
survobj <- Surv(groups.hi$time_to_bcr, groups.hi$bcr);
create.km.plot(survobj,
  groups.hi$group,
  xlab.label = 'Time (Years)',
  show.risktable = TRUE,
  filename = generate.filename('KM_tcga_GS910', 'pam_tcga_scale_sig_prd_class', 'pdf')
  );

####
mydata.ch <- mydata.ch.scale;
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
  filename = generate.filename('KM_changhai', 'pam_tcga_scale_sig', 'pdf')
  );
###
groups.ch[, 5:6] <- myclin.ch[match(rownames(groups.ch), myclin.ch$id), c('GS', 'Pathology.T')];
groups.ch.hi <- groups.ch[groups.ch$GS%in%c(9, 10), ];
mydata.ch <- data.frame(t(scale(t(mydata.ch.raw[, rownames(groups.ch.hi)]))));
groups.ch.hi$cor <- apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high']));
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
groups.ch.hi$cor <- apply(mydata.ch[genes, ], 2, function(x) cor(x, mydata.train$centroids[genes, 'high']));
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

