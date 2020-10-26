library(BoutrosLab.plotting.survival);
library(powerSurvEpi);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/F1')
clinical <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-01-05_tcga_survival.rds');
### The parameter is not suitable, just an example , do not use the following annotated parts 
#X1 <- rep(c(1, 0), each = nrow(clinical)/2);
#X2 <- sample(c(1, 0), length(X1), replace = TRUE);
#res <- numDEpi(X1, X2, power = 0.8, theta = 3.6, alpha = 0.05);
#res.power <- power.stratify(
#	n = nrow(clinical),
#	timeUnit = 1.2,
#	gVec = c(0.09, 0.91),
#	PVec = c(0.5, 0.5),
#	HR = 3,
#	lambda0Vec = c(2, 1), 
#	power.ini = 0.8,
#	power.low = 0.001,
#	power.upp = 0.999,
#	alpha = 0.05,
#	verbose = TRUE
#);
#mydata <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
#annot <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA.prad.clinical.rds');
#a1 <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA_PRAD_clinicalMatrix_Xena.txt', sep = '\t', header = TRUE);
#mysamp <- gsub('-01$', '', gsub('\\.', '-', colnames(mydata)));
#res.power <- power.stratify(
#	n = 87,
#	timeUnit = 1.2,
#	gVec = c(0.09, 0.91),
#	PVec = c(0.5, 0.5),
#	HR = 0.27,
#	lambda0Vec = c(2, 1), 
#	#power.ini = 0.8,
#	#power.low = 0.001,
#	#power.upp = 0.999,
#	alpha = 0.05,
#	verbose = TRUE
#);
#####
###
mydata <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
annot <- data.frame(readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA.prad.clinical.rds'));
a1 <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA_PRAD_clinicalMatrix_Xena.txt', sep = '\t', header = TRUE);
mysamp <- gsub('-01$', '', gsub('\\.', '-', colnames(mydata)));
a2 <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA_PRAD_survival_Xena.txt', header = TRUE);
a2 <- a2[grep('01', a2$sample), ];
####
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/cycle.in.all.cell.markers.posi.genes.txt';
clinical$time_to_bcr <- clinical$days_bcr/365;
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
###
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'cycle')), 'pdf')
        );
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 3, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 3.6, alpha = 0.05)$power;
#powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
#    nE = 134, nC = 136, RR = 3, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 3.6, alpha = 0.1)$power;

####
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));
mytest$gs <- annot[mytest$samp, ]$gleasonscore;
mytest$gs_group <- ifelse(mytest$gs%in%c(9, 10), 'h', 'l');
mytest$t <- annot[mytest$samp, ]$pathologyTstage;
mytest$t_group <- ifelse(grepl('t2', mytest$t), 'l', 'h');
###
survobj <- Surv(mytest[mytest$gs_group=='h', ]$time_to_bcr, mytest[mytest$gs_group=='h', ]$bcr);

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest[mytest$gs_group=='h', ]$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal')), 'pdf')
        );

survobj <- Surv(mytest[mytest$gs_group=='l', ]$time_to_bcr, mytest[mytest$gs_group=='l', ]$bcr);

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest[mytest$gs_group=='l', ]$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_low')), 'pdf')
        );

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_all')), 'pdf')
        );

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest[mytest$gs_group=='h', ]$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_flip')), 'pdf')
        );


mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;

mytest.dat <- mytest[mytest$gs_group=='l', ];
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.1)$power;

#powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
#    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 3, alpha = 0.05)$power;#

#powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
#    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.3, alpha = 0.05)$power;


mytest.dat <- mytest[mytest$gs_group=='h', ];
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.1)$power;

####

mytest.dat <- mytest[mytest$t_group=='l', ];
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.1)$power;



mytest.dat <- mytest[mytest$t_group=='h', ];
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.1)$power;

####

####
mypwr <- vector();
for(i in seq(1000)){
	set.seed(i);
	mytest$group <- sample(c('E', 'C'), nrow(mytest), replace = TRUE)
	mypwr[i] <- powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
	    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 3, alpha = 0.05)$power;
};

mypwr.h <- mypwr.l <- vector();
for(i in seq(1000)){
	set.seed(i);
	mytest.dat <- mytest[mytest$t_group=='l', ];
	mytest.dat$group <- sample(c('E', 'C'), nrow(mytest.dat), replace = TRUE)
	mypwr.l[i] <- powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
	    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
	
	###
	mytest.dat <- mytest[mytest$t_group=='h', ];
	mytest.dat$group <- sample(c('E', 'C'), nrow(mytest.dat), replace = TRUE)
	mypwr.h[i] <- powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest.dat, 
	    nE = nrow(mytest.dat[mytest.dat$group=='E', ]), nC = nrow(mytest.dat[mytest.dat$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
};
####
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
annot$samp <- paste0(gsub('-', '.', rownames(annot)), '.01');
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));
mytest$gs <- annot[mytest$samp, ]$gleasonscore;
mytest$gs_group <- ifelse(mytest$gs%in%c(9, 10), 'h', 'l');
mytest$t <- annot[mytest$samp, ]$pathologyTstage;
mytest$t_group <- ifelse(grepl('t2', mytest$t), 'l', 'h');
###
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_all')), 'pdf')
        );
mytest$PFI <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI;
mytest$PFI.time <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI.time/365;
survobj <- Surv(mytest$PFI.time, mytest$PFI);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean_PFS', '_', gsub('\\..*', '', 'basal_all')), 'pdf')
        );

mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
annot$samp <- paste0(gsub('-', '.', rownames(annot)), '.01');
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, colnames(mydata)%in%annot[annot$gleasonscore%in%c(9, 10), ]$samp];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));
mytest$gs <- annot[mytest$samp, ]$gleasonscore;
mytest$gs_group <- ifelse(mytest$gs%in%c(9, 10), 'h', 'l');
mytest$t <- annot[mytest$samp, ]$pathologyTstage;
mytest$t_group <- ifelse(grepl('t2', mytest$t), 'l', 'h');
###
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);

create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_G910')), 'pdf')
        );
mytest$PFI <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI;
mytest$PFI.time <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI.time/365;
survobj <- Surv(mytest$PFI.time, mytest$PFI);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean_PFS', '_', gsub('\\..*', '', 'basal_G910')), 'pdf')
        );

mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;

###
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
annot$samp <- paste0(gsub('-', '.', rownames(annot)), '.01');
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, colnames(mydata)%in%annot[annot$gleasonscore%in%c(6, 7, 8), ]$samp];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));
mytest$gs <- annot[mytest$samp, ]$gleasonscore;
mytest$gs_group <- ifelse(mytest$gs%in%c(9, 10), 'h', 'l');
mytest$t <- annot[mytest$samp, ]$pathologyTstage;
mytest$t_group <- ifelse(grepl('t2', mytest$t), 'l', 'h');
###
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_GS678')), 'pdf')
        );

mytest$PFI <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI;
mytest$PFI.time <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI.time/365;
survobj <- Surv(mytest$PFI.time, mytest$PFI);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean_PFS', '_', gsub('\\..*', '', 'basal_GS678')), 'pdf')
        );

mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
annot$samp <- paste0(gsub('-', '.', rownames(annot)), '.01');
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, colnames(mydata)%in%annot[grepl('t2', annot$pathologyTstage), ]$samp];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));
mytest$gs <- annot[mytest$samp, ]$gleasonscore;
mytest$gs_group <- ifelse(mytest$gs%in%c(9, 10), 'h', 'l');
mytest$t <- annot[mytest$samp, ]$pathologyTstage;
mytest$t_group <- ifelse(grepl('t2', mytest$t), 'l', 'h');
###
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_t2')), 'pdf')
        );
###
mytest$PFI <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI;
mytest$PFI.time <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI.time/365;
survobj <- Surv(mytest$PFI.time, mytest$PFI);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean_PFS', '_', gsub('\\..*', '', 'basal_t2')), 'pdf')
        );
####
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
###
mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
annot$samp <- paste0(gsub('-', '.', rownames(annot)), '.01');
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, colnames(mydata)%in%annot[!grepl('t2', annot$pathologyTstage), ]$samp];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));
mytest$gs <- annot[mytest$samp, ]$gleasonscore;
mytest$gs_group <- ifelse(mytest$gs%in%c(9, 10), 'h', 'l');
mytest$t <- annot[mytest$samp, ]$pathologyTstage;
mytest$t_group <- ifelse(grepl('t2', mytest$t), 'l', 'h');
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean', '_', gsub('\\..*', '', 'basal_t34')), 'pdf')
        );

mytest$PFI <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI;
mytest$PFI.time <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI.time/365;
survobj <- Surv(mytest$PFI.time, mytest$PFI);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean_PFS', '_', gsub('\\..*', '', 'basal_t34')), 'pdf')
        );

###
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;

####

mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/cycle.in.all.cell.markers.posi.genes.txt';
clinical$time_to_bcr <- clinical$days_bcr/365;
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
mytest$samp <- gsub('-01$', '', gsub('\\.', '-', rownames(mytest)));

###
mytest$PFI <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI;
mytest$PFI.time <- a2[match(mytest$samp, a2$X_PATIENT), ]$PFI.time/365;
survobj <- Surv(mytest$PFI.time, mytest$PFI);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_TCGA', paste0('mean_PFS', '_', gsub('\\..*', '', 'cycle')), 'pdf')
        );
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(PFI.time, PFI) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 3.6, alpha = 0.05)$power;
