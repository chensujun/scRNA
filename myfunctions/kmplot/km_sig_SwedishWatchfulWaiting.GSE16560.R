# command in linux 
# Rscript35 km_sig_tcga.R sb mylist.txt mean ; Rscript35 km_sig_cpc.R sb mylist.txt mean ; Rscript35 km_sig_MSKCC.Cancer.Cell.2010.R sb mylist.txt mean



library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mylist <- args[2];
mymethod <- args[3];
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
mydata <- readRDS('~/survival/data/2019-02-16_SwedishWatchfulWaiting.GSE16560.rds');
# str(mydata)
# 'data.frame':	60308 obs. of  144 variables:
#  $ CPCG0100: num  0 0.938 0 0 0 ...
#  $ CPCG0183: num  0 2.07 0 0 0 ...
#  $ CPCG0184: num  0 3.34 0 0 0 ...
# pdf("temp.hist.pdf")
# rownames(mydata)
# hist(mydata[,1])
# dev.off()
# ### it's not-log2ed exp matrix
clinical <- readRDS('~/survival/data/2019-02-16_SwedishWatchfulWaiting.GSE16560_survival.rds');
# str(clinical)
# 'data.frame':	144 obs. of  2 variables:
#  $ bcr        : logi  TRUE FALSE TRUE TRUE TRUE TRUE ...
#  $ time_to_bcr: num  0.337 9.155 0.942 2.686 0.649 ...
allEns <- readRDS('~/survival/data/2019-01-05_geneName_convert.rds');
rownames(clinical) <- gsub('-F1', '', rownames(clinical));

if(mytype=='sb'){
	mytest <- mydata[mygenes$V1, ];
}else{
	mytest <- mydata[rownames(allEns[allEns$gene%in%mygenes$V1, ]), ]
};
mytest <- data.frame(t(scale(t(mytest))));
if(all(is.na(mytest))){
	print(paste0('taking z-score of ', nrow(mytest), ' genes'))
}else{
	print(paste0('omiting NA values, from ', nrow(mytest), ' genes to ', nrow(na.omit(mytest))))
	mytest <- na.omit(mytest)
};
if(mymethod == 'mean'){
	mytest <- colMeans(mytest)
}else{
	mytest <- apply(mytest, 2, median)
};
mytest <- data.frame(mytest);
print(summary(mytest$mytest));
mytest$group <- ifelse(mytest$mytest>median(mytest$mytest), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
print(table(mytest$group));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);
pdf(generate.filename('KM_SwedishWatchfulWaiting.GSE16560', paste0(mymethod, '_', gsub('\\..*', '', mylist)), 'pdf'));
myplot;
dev.off();


