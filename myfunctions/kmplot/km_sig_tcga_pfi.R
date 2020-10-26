library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
#tcga count table from TCGA-PRAD-genc19.rpkm.csv
#tcga survival data from tcga.pdata.Rdata
#both data courtesy of Musa
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mygene <- args[2];
mydata <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
clinical <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/data/2020-09-02_tcga_survival_PFI.rds');
clinical$time_to_pfi <- clinical$PFI.time/365;
clinical$time_to_bcr <- clinical$days_bcr/365;
clinical$time_to_pfs <- clinical$PFS.time/365;
if(mytype=='ens'){
	stop('Please convert to HUGO gene name')
}else{
	mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
};
mytest <- mydata[mygenes$V1, ];
mytest <- data.frame(t(scale(t(mytest))));
if(all(is.na(mytest))){
	print(paste0('taking median z-score of ', nrow(mytest), ' genes'))
}else{
	print(paste0('omiting NA values, from ', nrow(mytest), ' genes to ', nrow(na.omit(mytest))))
	mytest <- na.omit(mytest)
};
if(mymethod == 'mean'){
	mytest <- colMeans(mytest)
}else if(mymethod == 'harmonic'){
	mytest <- apply(mytest, 2, function(x) mean(x^(-1))^(-1))
}else if(mymethod == 'min'){
	mytest <- apply(mytest, 2, min)
}else{
	mytest <- apply(mytest, 2, median)
};
mytest <- data.frame(mytest);
print(summary(mytest$mytest));
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
print(table(mytest$group));
mytest[, 3:4] <- clinical[rownames(mytest), c('PFI', 'time_to_pfi')];

survobj <- Surv(mytest$time_to_pfi, mytest$PFI);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);

pdf(generate.filename('KM_TCGA_PFI', paste0(mymethod, '_', gsub('\\..*', '', mylist)), 'pdf'));
myplot;
dev.off();

mytest[, 5:6] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);

pdf(generate.filename('KM_TCGA_BCR', paste0(mymethod, '_', gsub('\\..*', '', mylist)), 'pdf'));
myplot;
dev.off();


mytest[, 7:8] <- clinical[rownames(mytest), c('PFS', 'time_to_pfs')];
survobj <- Surv(mytest$time_to_pfs, mytest$PFS);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);

pdf(generate.filename('KM_TCGA_PFS', paste0(mymethod, '_', gsub('\\..*', '', mylist)), 'pdf'));
myplot;
dev.off();
