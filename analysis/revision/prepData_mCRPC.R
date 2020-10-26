library(BoutrosLab.plotting.general);
library(Seurat);
library(readxl);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig');
exp.m <- read.table('data/mCRPC.felix/2018_04_15_matrix_rna_tpm.txt', header = TRUE, as.is = TRUE, row.names = 1);
colnames(exp.m) <- gsub('.BL.*|.PRO.*', '', colnames(exp.m));
rownames(clinic) <- gsub('-', '.', clinic$Patient.ID);
clinic <- data.frame(read_excel('data/mCRPC.felix/SupplementalTable4.xlsx'));
clinic$OS.mCRPC <- as.numeric(as.vector(clinic$OS.mCRPC));
clinic$OS.mCRPC <- clinic$OS.mCRPC/365;
saveRDS(exp.m, file = 'data/mCRPC.felix/2020-06-29_mCRPC.felix.exp.nonLog2ed.rds');
saveRDS(clinic, file = 'data/mCRPC.felix/2020-06-29_mCRPC.felix.survival.rds');
####
exp.m <- read.table('data/mCRPC.2019/data_mRNA_seq_fpkm_polya.txt', header = TRUE, as.is = TRUE);
sym=exp.m$Hugo_Symbol
sym[duplicated(sym)]
exp.m[sym=="MARCH1",1:5]
exp.m[sym=="MARCH2",1:5]
exp.m[sym=="VAMP7",1:5]
exp.m[sym=="DHRSX",1:5]
exp.m[sym=="GTPBP6",1:5]
exp.m=aggregate(exp.m, by = list(sym), FUN = mean, na.rm = TRUE)
rownames(exp.m)=exp.m$"Group.1"
exp.m[1:5, 1:5]
exp.m=exp.m[,-c(1,2)]
exp.m[1:5, 1:5]

map= read.csv("./data/mCRPC.2019/data_clinical_sample.txt", stringsAsFactors=FALSE,row.names=1, sep="\t",skip = 4)# mapping from patient to sample
map[1:5, 1:2]
map$PATIENT_ID[duplicated(map$PATIENT_ID)]
rownames(map)=gsub("-", ".", rownames(map))
map=map[rownames(map)%in%colnames(exp.m),]
map$PATIENT_ID[duplicated(map$PATIENT_ID)]
map$PATIENT_ID[!map$PATIENT_ID%in%rownames(clinical)]
map=map[map$PATIENT_ID%in%rownames(clinical),]
map.1=rownames(map)
names(map.1)=as.character(map$PATIENT_ID)
clinical=clinical[rownames(clinical)%in%names(map.1),]
rownames(clinical)=map.1[rownames(clinical)]
map <- map[rownames(map)!='SC_9086_T', ];
shared.names=intersect(colnames(exp.m), rownames(clinical))
clinical=clinical[shared.names,]
exp.m=exp.m[,shared.names]
saveRDS(exp.m,'./data/mCRPC.2019/2020-06-29_mCRPC.2019.exp.m.nonLog2ed.polyA.rds')
exp.m=log2(exp.m+1)
saveRDS(exp.m,'./data/mCRPC.2019/2020-06-29_mCRPC.2019.exp.m.Log2ed.polyA.rds')
saveRDS(clinical,'./data/mCRPC.2019/2020-06-29_mCRPC.2019.suvival.polyA.rds')
