library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
#tcga count table from TCGA-PRAD-genc19.rpkm.csv
#tcga survival data from tcga.pdata.Rdata
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mygene <- args[2];
mydata <- readRDS('~/survival/data/2019-01-05_tcga_fpkm.rds');
clinical <- readRDS('~/survival/data/2019-01-05_tcga_survival.rds');
clinical$time_to_bcr <- clinical$days_bcr/365;
if(mytype=='ens'){
	stop('Please convert to HUGO gene name')
}else{
	mytest <- mydata[mygene, ];
};
mytest <- data.frame(t(mytest));
colnames(mytest) <- 'exp'
print(summary(mytest$exp));
mytest$group <- ifelse(mytest$exp>median(mytest$exp), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
print(table(mytest$group));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);
pdf(generate.filename('KM_TCGA', mygene, 'pdf'));
myplot;
dev.off();
