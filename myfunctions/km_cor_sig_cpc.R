library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mysig <- args[2];
mydata <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-01-05_cpc_fpkm.rds');
clinical <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/2019-01-05_cpc_survival.rds');
rownames(clinical) <- gsub('-F1', '', rownames(clinical));
conf <- read.config.file('~/landscape/rnaseq_landscape/master_config_rnaseq.R');
load(conf$feature.file);
mygenes <- readRDS(mysig);

if(mytype=='ens'){
	mygenes <- mygenes[names(mygenes)%in%rownames(mydata)];
}else{
	mygenes <- mygenes[names(mygenes)%in%allEns$gene];
	names(mygenes) <- rownames(allEns[match(names(mygenes), allEns$gene), ]);
};
mytest <- mydata[match(names(mygenes), rownames(mydata)), ];
#####
mytest <- data.frame(t(scale(t(mytest))));
if(all(!is.na(mytest))){
	print(paste0('taking z-score of ', nrow(mytest), ' genes'))
}else{
	print(paste0('omiting NA values, from ', nrow(mytest), ' genes to ', nrow(na.omit(mytest))))
	mytest <- na.omit(mytest)
};
mygenes <- mygenes[names(mygenes)%in%rownames(mytest)];
mytest.dat <- mytest;
mytest <- data.frame(cor = apply(mytest.dat, 2, function(x) cor(x, mygenes)));
mytest$pval <- apply(mytest.dat, 2, function(x) cor.test(x, mygenes)$p.value);

mytest$group <- ifelse(mytest$cor>median(mytest$cor), 'high', 'low');
mytest$group <- factor(mytest$group, levels = c('low', 'high'));
print(table(mytest$group));
mytest[, 4:5] <- clinical[rownames(mytest), c('bcr', 'time_to_bcr')];
saveRDS(mytest, file = generate.filename('KM_CPC', paste0('cor_', gsub('.*sig.|.rds', '', mysig)), 'rds'));
survobj <- Surv(mytest$time_to_bcr, mytest$bcr);
myplot <- create.km.plot(survival.truncate(list(survobj), time = 10)[[1]],
	mytest$group,
	xlab.label = 'Time (Years)',
	show.risktable = TRUE
	);
pdf(generate.filename('KM_CPC', paste0('cor_', gsub('.*sig.|.rds', '', mysig)), 'pdf'));
myplot;
dev.off();

