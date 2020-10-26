# command in linux 
# Rscript35 km_sig_tcga.R sb mylist.txt mean ; Rscript35 km_sig_cpc.R sb mylist.txt mean ; Rscript35 km_sig_MSKCC.Cancer.Cell.2010.R sb mylist.txt mean


library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
args <- commandArgs(trailingOnly = TRUE);
mytype <- args[1];
mylist <- args[2];
mymethod <- args[3];
mygenes <- read.table(mylist, as.is = TRUE, header = FALSE);
##############################################################
# download the expression data with the input gene list
library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/")
a = getCancerStudies(mycgds)
a[,2]
dataset.2check=c("Prostate Adenocarcinoma (MSKCC, Cancer Cell 2010)", "Prostate Adenocarcinoma CNA study (MSKCC, PNAS 2014)", "Metastatic Prostate Adenocarcinoma (MCTP, Nature 2012)")
# dataset 2 "Prostate Adenocarcinoma (MSKCC, Cancer Cell 2010)"
dataset.2check.id<-a[,1][a[,2] %in% dataset.2check]
dataset.2check.all.descrip<-a[a[,2] %in% dataset.2check,]
dataset.2check.all.descrip[,2]
dataset.2check.id[2] 
dataset.temp=dataset.2check.id[2]
c = getCaseLists(mycgds,dataset.temp)
c[,1]
c[,3]
c.temp= "prad_mskcc_mrna_primary"

profiles.temp=getGeneticProfiles(mycgds,dataset.temp)
exp.temp=profiles.temp[2,1]
#   "log2 whole transcript mRNA expression values (Affymetrix Human Exon 1.0 ST arrays)"
exp.temp

genes2read=mygenes$V1
# genes2read= read.csv( file = "./data/genes.scale.data.all.primary.cluster.up.regulated.txt", stringsAsFactors=FALSE,header = FALSE)
# genes2read=genes2read$V1
exp.matrix.2<-getProfileData(mycgds, genes2read, exp.temp, c.temp)
exp.matrix.2=data.frame(t(exp.matrix.2))
# str(exp.matrix.2)
# 'data.frame':	168 obs. of  131 variables:
#  $ PCA0001: num  9.63 8.35 10.74 9.02 10.64 ...
 # $ PCA0002: num  10.71 8.13 8.78 9.32 10.6 ...
samples.1=colnames(exp.matrix.2)

clinical.2 = getClinicalData(mycgds,c.temp)
# str(clinical.2)
# 'data.frame':	131 obs. of  21 variables:
#  $ CANCER_TYPE            : chr  "Prostate Cancer" "Prostate Cancer" "Prostate Cancer" "Prostate Cancer" ...
#  $ CANCER_TYPE_DETAILED   : chr  "Prostate Adenocarcinoma" "Prostate Adenocarcinoma" "Prostate Adenocarcinoma" "Prostate Adenocarcinoma" ...
samples.2=rownames(clinical.2)
samples=intersect(samples.1,samples.2)
exp.matrix.2=exp.matrix.2[,samples]
clinical.2=clinical.2[samples,c("DFS_STATUS","DFS_MONTHS")]
colnames(clinical.2)=c('bcr', 'time_to_bcr')
clinical.2$'time_to_bcr'=clinical.2$'time_to_bcr'/12
clinical.2$'bcr'=(clinical.2$'bcr'=="Recurred")
# download the expression data with the input gene list
##############################################################
saveRDS(exp.matrix.2,'./data/2019-02-09_MSKCC_medianZscore.rds')
saveRDS(clinical.2,'./data/2019-02-09_MSKCC_survival.rds')

mydata <- readRDS('./data/2019-02-09_MSKCC_medianZscore.rds');
mydata=2^mydata
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
clinical <- readRDS('./data/2019-02-09_MSKCC_survival.rds');
# str(clinical)
# 'data.frame':	144 obs. of  2 variables:
#  $ bcr        : logi  TRUE FALSE TRUE TRUE TRUE TRUE ...
#  $ time_to_bcr: num  0.337 9.155 0.942 2.686 0.649 ...
allEns <- readRDS('./data/2019-01-05_geneName_convert.rds');
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
pdf(generate.filename('KM_MSKCC.cancer.cell.2010', paste0(mymethod, '_', gsub('\\..*', '', mylist)), 'pdf'));
myplot;
dev.off();









