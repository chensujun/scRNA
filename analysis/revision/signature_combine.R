library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(Seurat);
library(pamr);
library(plyr);
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/sig/');
####
basal <- read.table('basal.in.all.cell.markers.posi.genes.txt');
cc <- read.table('cycle.in.all.cell.markers.posi.genes.txt');
basal91 <- read.table('data/survival.signatures/basal91.weight.txt');
pam.basal <- read.table('data/survival.signatures/pam50.basal.txt');
pam.a <- read.table('data/survival.signatures/pam50.lum.a.txt');
pam.b <- read.table('data/survival.signatures/pam50.lum.b.txt');
steml <- read.table('data/survival.signatures/stemlike.cancerres2017.txt');
pcs1 <- read.table('data/survival.signatures/pcs1.cancerres2016.txt');
### tcga
mydata.raw <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
colnames(mydata.raw) <- gsub('.[0-9]+$','', colnames(mydata.raw))
myclin <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2019-01-05_tcga_survival.rds');
myclin$time_to_bcr <- myclin$days_bcr/365;
rownames(myclin) <- gsub('.01', '', rownames(myclin))
###
mydata.ch.raw <- readRDS('pam/2020-08-19_filter_genes.FPKM_TCGA_match_T.rds');
myclin.ch <- read.csv('/cluster/projects/hansengroup/Public_Dataset_Hub/changhai/clinical/clinical.csv', header = TRUE);
myclin.ch$id <- paste0('T', myclin.ch$Patient.ID);
myclin.ch <- myclin.ch[, -16];
myclin.ch$time_to_bcr <- as.numeric(gsub('m', '', myclin.ch$month.to.surgery.BCR));
myclin.ch$bcr <- ifelse(myclin.ch$BCR=='YES', TRUE, FALSE);
rownames(myclin.ch) <- paste0('T', myclin.ch$Patient.ID)
###
mydata.eb <- readRDS('data/EBioMed/2019-02-10_EBioMedicine2015.rds');
mydata.eb <- 2^mydata.eb;
mydata.eb.raw <- mydata.eb
myclin.eb <- readRDS('data/EBioMed/2019-02-10_EBioMedicine2015_survival.rds');
myclin.gs <- readRDS('data/EBioMed/2019-05-02_EBioMedicine2015_gleason.rds');
myclin.t <- readRDS('data/EBioMed/2019-06-07_EBioMedicine2015_tstage.rds');
myclin.eb$gs <- myclin.gs[rownames(myclin.eb)];
myclin.eb$t <- myclin.t[rownames(myclin.eb)];
###
####
calc_score <- function(mytest, method.scale = 'mean', mymethod = 'mean', mygenes, clinical, label = '', cor.method = 'spearman'){
	if(method.scale == 'mean'){
		print(paste0('taking mean z-score, from ', nrow(mytest), ' genes to ', nrow(na.omit(mytest))))
		mytest <- na.omit(data.frame(t(scale(t(mytest)))));
	}else{
		print(paste0('taking mean z-score, from ', nrow(mytest), ' genes to ', nrow(na.omit(mytest))))		
		mytest <- na.omit(data.frame(t(scale(t(mytest), center = apply(mytest, 1, median)))))
	}

	mygenes <- mygenes[rownames(mytest)];
	if(mymethod == 'mean'){
		mytest <- colMeans(mytest)
	}else if (mymethod == 'median') {
		mytest <- apply(mytest, 2, median)
	}else if (mymethod == 'cor') {
		mytest.dat <- mytest;
		mytest <- data.frame(mytest = apply(mytest.dat, 2, function(x) cor(x, mygenes, method = cor.method)));
	}else if (mymethod == 'weighted') {
		mytest.dat <- data.frame(t(sapply(seq(length(mygenes)), function(x) unlist(mytest[x, ])*mygenes[x])));
		rownames(mytest.dat) <- names(mygenes);
		mytest <- colMeans(mytest.dat)
		} else {
			print('Please specify a method for signature calculation')
		};
	mytest <- data.frame(mytest);
	print(summary(mytest$mytest));
	mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
	mytest$group <- factor(mytest$group, levels = c('low', 'high'));
	print(table(mytest$group));
	mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
	colnames(mytest) <- paste0(colnames(mytest), label);
	return(mytest);
};

mytest <- mydata.raw[rownames(mydata.raw)%in%basal$V1, ];
mytest.tcga <- calc_score(mytest, 'mean', 'mean', basal$V1, myclin, '_basal');
mytest.eb <- calc_score(mydata.eb[rownames(mydata.eb)%in%basal$V1, ], 'mean', 'mean', basal$V1, myclin.eb, '_basal');
mytest.ch <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%basal$V1, ], 'mean', 'mean', basal$V1, myclin.ch, '_basal');
####
mytest <- mydata.raw[rownames(mydata.raw)%in%cc$V1, ];
mytest.tcga[, (ncol(mytest.tcga)+1):(ncol(mytest.tcga)+2)] <- calc_score(mytest, 'mean', 'mean', cc$V1, myclin, '_cycle')[, 1:2];
mytest.eb[, (ncol(mytest.eb)+1):(ncol(mytest.eb)+2)] <- calc_score(mydata.eb[rownames(mydata.eb)%in%cc$V1, ], 
	'mean', 'mean', cc$V1, myclin.eb, '_cycle')[, 1:2];
mytest.ch[, (ncol(mytest.ch)+1):(ncol(mytest.ch)+2)] <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%cc$V1, ], 
	'mean', 'mean', cc$V1, myclin.ch, '_cycle')[, 1:2];
###
mygenes <- basal91$V2;
names(mygenes) <- basal91$V1;
mytest <- mydata.raw[names(mygenes), ];
mytest.tcga[, (ncol(mytest.tcga)+1):(ncol(mytest.tcga)+2)] <- calc_score(mytest, 'mean', 'weighted', mygenes, myclin, '_basal91')[, 1:2];
mytest.eb[, (ncol(mytest.eb)+1):(ncol(mytest.eb)+2)] <- calc_score(mydata.eb[rownames(mydata.eb)%in%basal91$V1, ], 
	'mean', 'weighted', mygenes, myclin.eb, '_basal91')[, 1:2];
mytest.ch[, (ncol(mytest.ch)+1):(ncol(mytest.ch)+2)] <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%basal91$V1, ], 
	'mean', 'weighted', mygenes, myclin.ch, '_basal91')[, 1:2];
####
mytest <- mydata.raw[rownames(mydata.raw)%in%pcs1$V1, ];
mytest.tcga[, (ncol(mytest.tcga)+1):(ncol(mytest.tcga)+2)] <- calc_score(mytest, 'mean', 'mean', pcs1$V1, myclin, '_pcs1')[, 1:2];
mytest.eb[, (ncol(mytest.eb)+1):(ncol(mytest.eb)+2)] <- calc_score(mydata.eb[rownames(mydata.eb)%in%pcs1$V1, ], 
	'mean', 'mean', pcs1$V1, myclin.eb, '_pcs1')[, 1:2];
mytest.ch[, (ncol(mytest.ch)+1):(ncol(mytest.ch)+2)] <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%pcs1$V1, ], 
	'mean', 'mean', pcs1$V1, myclin.ch, '_pcs1')[, 1:2];
##
mygenes <- pam.basal$V2;
names(mygenes) <- pam.basal$V1;
mytest <- mydata.raw[rownames(mydata.raw)%in%pam.basal$V1, ];
mytest.tcga[, (ncol(mytest.tcga)+1):(ncol(mytest.tcga)+2)] <- calc_score(mytest, 'median', 'cor', mygenes, myclin, '_pamBasal')[, 1:2];
mytest.eb[, (ncol(mytest.eb)+1):(ncol(mytest.eb)+2)] <- calc_score(mydata.eb[rownames(mydata.eb)%in%pam.basal$V1, ], 
	'median', 'cor', mygenes, myclin.eb, '_pamBasal')[, 1:2];
mytest.ch[, (ncol(mytest.ch)+1):(ncol(mytest.ch)+2)] <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%pam.basal$V1, ], 
	'median', 'cor', mygenes, myclin.ch, '_pamBasal')[, 1:2];
##
mygenes <- pam.a$V2;
names(mygenes) <- pam.a$V1;
mytest <- mydata.raw[rownames(mydata.raw)%in%pam.a$V1, ];
mytest.tcga[, (ncol(mytest.tcga)+1):(ncol(mytest.tcga)+2)] <- calc_score(mytest, 'median', 'cor', mygenes, myclin, '_pamA')[, 1:2];
mytest.eb[, (ncol(mytest.eb)+1):(ncol(mytest.eb)+2)] <- calc_score(mydata.eb[rownames(mydata.eb)%in%pam.a$V1, ], 
	'median', 'cor', mygenes, myclin.eb, '_pamA')[, 1:2];
mytest.ch[, (ncol(mytest.ch)+1):(ncol(mytest.ch)+2)] <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%pam.a$V1, ], 
	'median', 'cor', mygenes, myclin.ch, '_pamA')[, 1:2];
##
mygenes <- pam.b$V2;
names(mygenes) <- pam.b$V1;
mytest <- mydata.raw[rownames(mydata.raw)%in%pam.basal$V1, ];
mytest.tcga[, (ncol(mytest.tcga)+1):(ncol(mytest.tcga)+2)] <- calc_score(mytest, 'median', 'cor', mygenes, myclin, '_pamB')[, 1:2];
mytest.eb[, (ncol(mytest.eb)+1):(ncol(mytest.eb)+2)] <- calc_score(mydata.eb[rownames(mydata.eb)%in%pam.b$V1, ], 
	'median', 'cor', mygenes, myclin.eb, '_pamB')[, 1:2];
mytest.ch[, (ncol(mytest.ch)+1):(ncol(mytest.ch)+2)] <- calc_score(mydata.ch.raw[rownames(mydata.ch.raw)%in%pam.b$V1, ], 
	'median', 'cor', mygenes, myclin.ch, '_pamB')[, 1:2];
####
mytest.tcga$group_basal_basal91 <- paste0(mytest.tcga$group_basal, '_', mytest.tcga$group_basal91);
mytest.tcga$group_basal_cycle <- paste0(mytest.tcga$group_basal, '_', mytest.tcga$group_cycle)
mytest.tcga$group_basal_pcs1 <- paste0(mytest.tcga$group_basal, '_', mytest.tcga$group_pcs1)
mytest.tcga$group_basal_pamBasal <- paste0(mytest.tcga$group_basal, '_', mytest.tcga$group_pamBasal)
mytest.tcga$group_basal_pamA <- paste0(mytest.tcga$group_basal, '_', mytest.tcga$group_pamA);
mytest.tcga$group_basal_pamB <- paste0(mytest.tcga$group_basal, '_', mytest.tcga$group_pamB);
#
mytest.tcga$group_cycle_basal91 <- paste0(mytest.tcga$group_cycle, '_', mytest.tcga$group_basal91);
mytest.tcga$group_cycle_pcs1 <- paste0(mytest.tcga$group_cycle, '_', mytest.tcga$group_pcs1)
mytest.tcga$group_cycle_pamBasal <- paste0(mytest.tcga$group_cycle, '_', mytest.tcga$group_pamBasal)
mytest.tcga$group_cycle_pamA <- paste0(mytest.tcga$group_cycle, '_', mytest.tcga$group_pamA);
mytest.tcga$group_cycle_pamB <- paste0(mytest.tcga$group_cycle, '_', mytest.tcga$group_pamB);
####
mytest.eb$group_basal_basal91 <- paste0(mytest.eb$group_basal, '_', mytest.eb$group_basal91);
mytest.eb$group_basal_cycle <- paste0(mytest.eb$group_basal, '_', mytest.eb$group_cycle)
mytest.eb$group_basal_pcs1 <- paste0(mytest.eb$group_basal, '_', mytest.eb$group_pcs1)
mytest.eb$group_basal_pamBasal <- paste0(mytest.eb$group_basal, '_', mytest.eb$group_pamBasal)
mytest.eb$group_basal_pamA <- paste0(mytest.eb$group_basal, '_', mytest.eb$group_pamA);
mytest.eb$group_basal_pamB <- paste0(mytest.eb$group_basal, '_', mytest.eb$group_pamB);
#
mytest.eb$group_cycle_basal91 <- paste0(mytest.eb$group_cycle, '_', mytest.eb$group_basal91);
mytest.eb$group_cycle_pcs1 <- paste0(mytest.eb$group_cycle, '_', mytest.eb$group_pcs1)
mytest.eb$group_cycle_pamBasal <- paste0(mytest.eb$group_cycle, '_', mytest.eb$group_pamBasal)
mytest.eb$group_cycle_pamA <- paste0(mytest.eb$group_cycle, '_', mytest.eb$group_pamA);
mytest.eb$group_cycle_pamB <- paste0(mytest.eb$group_cycle, '_', mytest.eb$group_pamB);
####
mytest.ch$group_basal_basal91 <- paste0(mytest.ch$group_basal, '_', mytest.ch$group_basal91);
mytest.ch$group_basal_cycle <- paste0(mytest.ch$group_basal, '_', mytest.ch$group_cycle)
mytest.ch$group_basal_pcs1 <- paste0(mytest.ch$group_basal, '_', mytest.ch$group_pcs1)
mytest.ch$group_basal_pamBasal <- paste0(mytest.ch$group_basal, '_', mytest.ch$group_pamBasal)
mytest.ch$group_basal_pamA <- paste0(mytest.ch$group_basal, '_', mytest.ch$group_pamA);
mytest.ch$group_basal_pamB <- paste0(mytest.ch$group_basal, '_', mytest.ch$group_pamB);
#
mytest.ch$group_cycle_basal91 <- paste0(mytest.ch$group_cycle, '_', mytest.ch$group_basal91);
mytest.ch$group_cycle_pcs1 <- paste0(mytest.ch$group_cycle, '_', mytest.ch$group_pcs1)
mytest.ch$group_cycle_pamBasal <- paste0(mytest.ch$group_cycle, '_', mytest.ch$group_pamBasal)
mytest.ch$group_cycle_pamA <- paste0(mytest.ch$group_cycle, '_', mytest.ch$group_pamA);
mytest.ch$group_cycle_pamB <- paste0(mytest.ch$group_cycle, '_', mytest.ch$group_pamB);
####
####
plot_km <- function(mytest, group, datatype, sigtype){
	survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
	myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
		mytest[, group],
		xlab.label = 'Time (Years)',
		show.risktable = TRUE,
		generate.filename(paste0('KM_', datatype), sigtype, 'pdf')
		);
};
plot_km(mytest.tcga, 'group_basal_basal91', 'tcga', 'basal_basal91');
plot_km(mytest.tcga, 'group_basal_cycle', 'tcga', 'basal_cycle');
plot_km(mytest.tcga, 'group_basal_pcs1', 'tcga', 'basal_pcs1');
plot_km(mytest.tcga, 'group_basal_pamBasal', 'tcga', 'basal_pamBasal');
plot_km(mytest.tcga, 'group_basal_pamA', 'tcga', 'basal_pamA');
plot_km(mytest.tcga, 'group_basal_pamB', 'tcga', 'basal_pamB');

plot_km(mytest.tcga, 'group_cycle_basal91', 'tcga', 'cycle_basal91');
plot_km(mytest.tcga, 'group_cycle_pcs1', 'tcga', 'cycle_pcs1');
plot_km(mytest.tcga, 'group_cycle_pamBasal', 'tcga', 'cycle_pamBasal');
plot_km(mytest.tcga, 'group_cycle_pamA', 'tcga', 'cycle_pamA');
plot_km(mytest.tcga, 'group_cycle_pamB', 'tcga', 'cycle_pamB');

plot_km(mytest.eb, 'group_basal_basal91', 'eb', 'basal_basal91');
plot_km(mytest.eb, 'group_basal_cycle', 'eb', 'basal_cycle');
plot_km(mytest.eb, 'group_basal_pcs1', 'eb', 'basal_pcs1');
plot_km(mytest.eb, 'group_basal_pamBasal', 'eb', 'basal_pamBasal');
plot_km(mytest.eb, 'group_basal_pamA', 'eb', 'basal_pamA');
plot_km(mytest.eb, 'group_basal_pamB', 'eb', 'basal_pamB');

plot_km(mytest.eb, 'group_cycle_basal91', 'eb', 'cycle_basal91');
plot_km(mytest.eb, 'group_cycle_pcs1', 'eb', 'cycle_pcs1');
plot_km(mytest.eb, 'group_cycle_pamBasal', 'eb', 'cycle_pamBasal');
plot_km(mytest.eb, 'group_cycle_pamA', 'eb', 'cycle_pamA');
plot_km(mytest.eb, 'group_cycle_pamB', 'eb', 'cycle_pamB');

plot_km(mytest.ch, 'group_basal_basal91', 'ch', 'basal_basal91');
plot_km(mytest.ch, 'group_basal_cycle', 'ch', 'basal_cycle');
plot_km(mytest.ch, 'group_basal_pcs1', 'ch', 'basal_pcs1');
plot_km(mytest.ch, 'group_basal_pamBasal', 'ch', 'basal_pamBasal');
plot_km(mytest.ch, 'group_basal_pamA', 'ch', 'basal_pamA');
plot_km(mytest.ch, 'group_basal_pamB', 'ch', 'basal_pamB');

plot_km(mytest.ch, 'group_cycle_basal91', 'ch', 'cycle_basal91');
plot_km(mytest.ch, 'group_cycle_pcs1', 'ch', 'cycle_pcs1');
plot_km(mytest.ch, 'group_cycle_pamBasal', 'ch', 'cycle_pamBasal');
plot_km(mytest.ch, 'group_cycle_pamA', 'ch', 'cycle_pamA');
plot_km(mytest.ch, 'group_cycle_pamB', 'ch', 'cycle_pamB');

