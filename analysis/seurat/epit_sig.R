library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(Hmisc);
library(Seurat);
mydata <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-01-05_tcga_fpkm.rds');
sig.basal <- c('WFDC2', 'SLPI','SFN', 'GADD45A', 'RARRES1');
sig.cycle <- c('STMN1', 'UBE2C');

get_group <- function(mydata, mygene){
	mytest <- mydata[mygene, ];
	mytest <- data.frame(t(scale(t(mytest))));
	mytest <- colMeans(mytest)
	mytest <- data.frame(mean = mytest);
	#mytest$group <- ifelse(mytest$mean>median(mytest$mean), 'high', 'low');
	#mytest$group <- factor(mytest$group, levels = c('low', 'high'));
	mytest$group <- factor(as.numeric(cut2(mytest$mean, g = 3)));
	return(mytest)
}

mytest.basal <- get_group(mydata, sig.basal)
mytest.cycle <- get_group(mydata, sig.cycle)

mytest <- cbind(mytest.basal, mytest.cycle);
colnames(mytest)[2:3] <- c('group_basal', 'group_cycle');
#colnames(mytest) <- paste0(colnames(mytest), rep(c('_basal', '_cycle'), each = 2));
#mytest$group <- paste0(mytest$group_basal, '_', mytest$group_cycle);
#mytest$diff <- mytest$mean_basal - mytest$mean_cycle;
#mytest$group <- ifelse(mytest$diff>median(mytest$diff), 'high', 'low')
mytest$group <- paste0(mytest$group_basal, '_', mytest$group_cycle);
mytest[, (ncol(mytest)+1):(ncol(mytest)+2)] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE
        );
pdf(generate.filename('KM_TCGA', paste0('basal_cycle'), 'pdf'));
myplot;
dev.off();

###
mydata <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-04-03_MSKCC_medianZscore.rds');
clinical <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-02-09_MSKCC_survival.rds');
mydata <- mydata[, -1];
mytest.basal <- get_group(mydata, sig.basal)
mytest.cycle <- get_group(mydata, sig.cycle)

mytest <- cbind(mytest.basal, mytest.cycle[, 'group']);
colnames(mytest)[2:3] <- c('group_basal', 'group_cycle');
mytest$group <- paste0(mytest$group_basal, '_', mytest$group_cycle);
mytest[, 5:6] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];

survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
        mytest$group,
        xlab.label = 'Time (Years)',
        show.risktable = TRUE
        );
pdf(generate.filename('KM_MSKCC', paste0('basal_cycle'), 'pdf'));
myplot;
dev.off();
