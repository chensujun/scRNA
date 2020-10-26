library(BoutrosLab.plotting.survival);
library(powerSurvEpi);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/F1')

mylist <- '/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig/basal.in.all.cell.markers.posi.genes.txt';
a1 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-05-02_EBioMedicine2015_gleason.rds');
a2 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-02-10_EBioMedicine2015_survival.rds');
a3 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-06-07_EBioMedicine2015_tstage.rds');
mydata <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/2019-02-10_EBioMedicine2015.rds');
mydata <- 2^mydata;
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mygenes <- mygenes[mygenes$V1%in%rownames(mydata), , drop = FALSE];
mytest <- mydata[mygenes$V1, colnames(mydata)%in%names(na.omit(a1[a1%in%c(5,6,7,8)]))];
#mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
#mytest <- mytest[, colnames(mydata)%in%names(na.omit(a1[a1%in%c(5,6,7,8)]))];
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- a2[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_EBioM', paste0('mean_1', '_', gsub('\\..*', '', 'basal_GS5678')), 'pdf')
        );

###
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
mytest <- mydata[mygenes$V1, colnames(mydata)%in%names(na.omit(a1[!a1%in%c(5,6,7,8)]))];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- a2[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_EBioM', paste0('mean', '_', gsub('\\..*', '', 'basal_GS910')), 'pdf')
        );

###
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;

###
###
mytest <- mydata[mygenes$V1, colnames(mydata)%in%names(na.omit(a3[a3%in%c('T2')]))];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- a2[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_EBioM', paste0('mean', '_', gsub('\\..*', '', 'basal_T2')), 'pdf')
        );

###
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
####
###
mytest <- mydata[mygenes$V1, colnames(mydata)%in%names(na.omit(a3[a3%in%c('T3', 'T4')]))];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- a2[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_EBioM', paste0('mean', '_', gsub('\\..*', '', 'basal_T34')), 'pdf')
        );

###
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
###
mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
mytest <- na.omit(mytest);
mytest <- colMeans(mytest);
mytest <- data.frame(mytest);
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
mytest[, 3:4] <- a2[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE,
        filename = generate.filename('KM_EBioM', paste0('mean', '_', gsub('\\..*', '', 'basal_all')), 'pdf')
        );

###
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'E', 'C');
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.05)$power;
powerCT(formula = Surv(time_to_bcr, bcr) ~ group, dat = mytest, 
    nE = nrow(mytest[mytest$group=='E', ]), nC = nrow(mytest[mytest$group=='C', ]), RR = 0.27, alpha = 0.1)$power;
