library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mygene <- args[2];
mydata <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary//data/2019-01-05_cpc_fpkm.rds');
clinical <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary//data/2019-01-05_cpc_survival.rds');
allEns <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary//data/2019-01-05_geneName_convert.rds');
rownames(clinical) <- gsub('-F1', '', rownames(clinical));
if(mytype=='ens'){
	mytest <- mydata[mygene, ]	
}else{
	mytest <- mydata[rownames(allEns[allEns$gene==mygene, ]), ]
};
mytest <- data.frame(t(mytest));
colnames(mytest) <- 'exp'
print(summary(mytest$exp));
mytest$group <- ifelse(mytest$exp>median(mytest$exp), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
print(table(mytest$group));
mytest[, 3:4] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
pdf(generate.filename('KM_CPC', mygene, 'pdf'));
create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);
dev.off();
