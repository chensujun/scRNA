library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
#tcga count table from TCGA-PRAD-genc19.rpkm.csv
#tcga survival data from tcga.pdata.Rdata
#both data courtesy of Musa
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mysig <- args[2];
mydata <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-01-05_tcga_fpkm.rds');
clinical <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-01-05_tcga_survival.rds');
clinical$time_to_bcr <- clinical$days_bcr/365;
if(mytype=='ens'){
	stop('Please convert to HUGO gene name')
}else{
	mygenes <- readRDS(mysig);
};
mygenes <- mygenes[names(mygenes)%in%rownames(mydata)];
mytest <- mydata[names(mygenes), ];
mytest <- data.frame(t(scale(t(mytest))));
if(all(is.na(mytest))){
	print(paste0('taking median z-score of ', nrow(mytest), ' genes'))
}else{
	print(paste0('omiting NA values, from ', nrow(mytest), ' genes to ', nrow(na.omit(mytest))))
	mytest <- na.omit(mytest)
};

mytest.dat <- mytest;
mytest <- data.frame(cor = apply(mytest.dat, 2, function(x) cor(x, mygenes)));
mytest$pval <- apply(mytest.dat, 2, function(x) cor.test(x, mygenes)$p.value);

mytest$group <- ifelse(mytest$cor>median(mytest$cor), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
print(table(mytest$group));
mytest[, 4:5] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
mytest <- na.omit(mytest);
saveRDS(mytest, file = generate.filename('KM_TCGA', paste0('cor_', gsub('.*sig.|.rds', '', mysig)), 'rds'));
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);
pdf(generate.filename('KM_TCGA', paste0('cor_', gsub('.*sig.|.rds', '', mysig)), 'pdf'));
myplot;
dev.off();
