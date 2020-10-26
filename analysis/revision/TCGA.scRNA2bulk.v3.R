require(scran)
require(Seurat)
require(limma)
require(TCGA2STAT)

### read TCGA data
# prad.RNASeq=readRDS("TCGA.prad.RNASeq.rds")# this includes samples with "daystopsa"
# prad.RNASeq <- getTCGA(disease="PRAD", data.type="RNASeq2", clinical=TRUE, cvars = "gleasonscore", type="RPKM")
# saveRDS(prad.RNASeq, "TCGA.prad.RNASeq.RPKM.with.Gleason.rds")
prad.RNASeq=readRDS("TCGA.prad.RNASeq.RPKM.with.Gleason.rds")
prad.norm <- TumorNormalMatch(prad.RNASeq$dat)
prad.norm = prad.norm$normal
colnames(prad.norm)=paste(colnames(prad.norm), "-N", sep="")
table(rownames(prad.norm)==colnames(prad.RNASeq$merged.dat)[-c(1,2)])
names.shared=intersect(rownames(prad.norm), colnames(prad.RNASeq$merged.dat))
prad.norm=prad.norm[names.shared, ]
prad.norm=t(prad.norm)
gleason.norm=rep("normal", nrow(prad.norm))
names(gleason.norm)=rownames(prad.norm)
# clinical
gleason= prad.RNASeq$merged.dat$GLEASONSCORE
names(gleason)=prad.RNASeq$merged.dat$bcr
exp.m = prad.RNASeq$merged.dat
rownames(exp.m)=prad.RNASeq$merged.dat$bcr
exp.m[1:5,1:5]
exp.m=exp.m[, -c(1,2)]
exp.m=exp.m[, names.shared]
exp.m[1:5,1:5]
str(exp.m[1:5, 1:5])
exp.m=as.matrix(exp.m)
exp.m=rbind(exp.m, prad.norm)
gleason=c(gleason, gleason.norm)
table(colnames(exp.m)==names(gleason))

summary(as.numeric(exp.m[1,]))
summary(as.numeric(exp.m[11,]))
summary(as.numeric(exp.m[111,]))
summary(as.numeric(exp.m[,1]))
summary(as.numeric(exp.m[,11]))
summary(as.numeric(exp.m[,111]))
hist(as.numeric(exp.m[1,]))
dev.off()
hist(log2(as.numeric(exp.m[1,])))
dev.off()

exp.m=log2(exp.m+1)
exp.m = t(scale(exp.m, scale = TRUE))
# exp.m = t( exp.m )
summary(apply(exp.m, 2, median))
table(names(gleason)==colnames(exp.m))


# calcualte cell abundance for every patient
# exp.m=readRDS('./TCGA.withGleason.log2.centered.NonScaled.sharedGeneWithCsRNAseq.rds')
# # puri.m=readRDS("TCGA.exp.m.purified.rds")
# gleason=readRDS('./TCGA.Gleason.withGleason.rds')
# # created by TCGA.mappingGleason2cells_v2_puriTCGA2subtypes.R

# exp.m <- readRDS('./survival/data/2019-01-05_tcga_fpkm.rds');
# clin.t <- data.frame(readRDS('./TCGA.prad.clinical.rds'));
# exp.m=as.matrix(exp.m)
# exp.m=log2(exp.m+1)
# colnames(exp.m) <- gsub('.[0-9]+$','', colnames(exp.m))
# gleason=as.character(clin.t$gleasonscore)
# names(gleason)=rownames(clin.t)
# names(gleason) <- gsub('-', '.', names(gleason))
# samples.shared=intersect(names(gleason), colnames(exp.m))
# exp.m=exp.m[,samples.shared]
# gleason=gleason[samples.shared]
# str(gleason)
# str(exp.m)

### read single-cell data
# epiTumor=readRDS("epiTumor.basal.batchEffectRemovedMatirx.seuratClustered.rds")
# epiTumor<-ScaleData(object = epiTumor)
# epiTumor@meta.data$gleason=rep("9", length(epiTumor@cell.names))
# epiTumor@meta.data$gleason[epiTumor@meta.data$fig.patient %in% c("JD1800159SL", "JD1800162SL", "JD1800174SL", "JD1800176SL")]="7"

# pros13.markers = readRDS("pros13.finerCellTypesFromSusu.finalVersion.Markers.batchEffectRemovedMatirx.seuratClustered.rds")
# pros13.celltype.markers = readRDS("pros13.finerCellTypesFromSusu.finalVersion.celltypeMarkers.batchEffectRemovedMatirx.seuratClustered.rds")
# pros13 = readRDS("pros13.finerCellTypesFromSusu.finalVersion.SeuratObj.batchEffectRemovedMatirx.seuratClustered.rds")
# pros13@meta.data = readRDS("pros13.finerCellTypesFromSusu.finalVersion.SeuratMetadata.batchEffectRemovedMatirx.seuratClustered.rds")
# pros13<-ScaleData(object = pros13)
# pros13@meta.data$gleason=rep("9", length(pros13@cell.names))
# pros13@meta.data$gleason[pros13@meta.data$fig.patient %in% c("JD1800159SL", "JD1800162SL", "JD1800174SL", "JD1800176SL")]="7"

genes.all.sc = intersect(rownames(exp.m), rownames(pros13@scale.data))
exp.m = exp.m[genes.all.sc, ]

### summarize profiles for TCGA gs7 and gs9's patients
patients=unique(pros13@meta.data$fig.patient)
# patients=setdiff(unique(pros13@meta.data$fig.patient), "JD1800154SL")
sum.profile=matrix(0, nrow=nrow(pros13@raw.data), ncol=length(patients))
colnames(sum.profile)=patients
rownames(sum.profile)=rownames(pros13@raw.data)
# epi.pct=c(0.6, 0.7, 0.8, 0.8, 0.8, 0.9, 0.1, 0.7, 0.9, 0.6, 0.8, 0.8, 0.5)
# names(epi.pct)=c("JD1800153SL", "JD1800154SL", "JD1800155SL", "JD1800156SL", "JD1800159SL", "JD1800162SL", "JD1800171SL", "JD1800172SL", "JD1800173SL", "JD1800174SL", "JD1800175SL", "JD1800176SL", "JD1800177SL")
for(i in patients)
{
	### as mean of scaled expression across cells 
	# temp.profile=rowMeans(epiTumor@scale.data[genes.all.sc, epiTumor@meta.data$fig.patient==i])
	### as loged summarized counts of all cells
	temp.profile=rowSums(pros13@raw.data[, pros13@meta.data$fig.patient==i])
	sum.profile[, i]=log2(temp.profile/sum(temp.profile)*1000000+1)
	# sum.profile[, i]=temp.profile/sum(temp.profile)*1000+1
	### as loged summarized counts of all cells, but weighted by epithelial percentage
	# temp.profile.epi=rowSums(pros13@raw.data[, (pros13@meta.data$fig.patient==i)&(pros13@meta.data$finer.cell.type%in%c("Luminal", "Basal/intermediate"))])
	# temp.profile.str=rowSums(pros13@raw.data[, (pros13@meta.data$fig.patient==i)&(!pros13@meta.data$finer.cell.type%in%c("Luminal", "Basal/intermediate"))])
	# temp.profile.epi=temp.profile.epi/sum(temp.profile.epi)*1000*epi.pct[i]
	# temp.profile.str=temp.profile.str/sum(temp.profile.str)*1000*(1-epi.pct[i])
	# sum.profile[, i]=log2(temp.profile.epi+temp.profile.str+1)
}
sd.sum.profile=apply(sum.profile, 1, sd)
sum.profile=sum.profile[sd.sum.profile>0,]
sum.profile = t(scale(t(sum.profile) , scale = TRUE))


# get DEGs between summarized profiles of GS 9 and GS 7
mapP2GS=rep("9", length(unique(pros13@meta.data$fig.patient)))
names(mapP2GS)=unique(pros13@meta.data$fig.patient)
mapP2GS[c("JD1800159SL", "JD1800162SL", "JD1800174SL", "JD1800176SL")]="7"
gs=mapP2GS[colnames(sum.profile)]
str(sum.profile)
gs
table(names(gs)==colnames(sum.profile))
genes=rownames(sum.profile)
deg.sum.p=matrix(0, nrow=length(genes), ncol=2)
colnames(deg.sum.p)=c("t.test.pvalue", "log2fc")
rownames(deg.sum.p)=genes
for(i in genes)
{
	e1 = sum.profile[i, gs=="9"]
	e2 = sum.profile[i, gs=="7"]
	if(sd(c(e1,e2))==0)
	{
		pvalue=1
	}
	else
	{
		pvalue=wilcox.test(e1, e2)$p.value
	}
	log2fc=mean(e1)-mean(e2)
	deg.sum.p[i,]=c(pvalue, log2fc)
}
# abs???
deg.sum.p=deg.sum.p[(deg.sum.p[,1]<0.05),] # & (deg.sum.p[,2])>log2(1.5)
str(deg.sum.p)
deg.sum.p=deg.sum.p[rownames(deg.sum.p)%in%genes.all.sc, ]
str(deg.sum.p)


### DEGs.9:7 - DEGs.7:6, simpleDEG, consider +/- normal samples, NON-Purified expression data, compare one gleason level to all others, All Cell Types
### ----------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
str(exp.m)
table(gleason)
gleason[gleason%in%c("6","7")]="GS 6, 7"
gleason[gleason%in%c("8","9","10")]="GS 8, 9, 10"

table(colnames(exp.m)==names(gleason))
table(gleason)
gleason = factor(gleason)


### check correlation between summarized profiles and TCGA samples'
gs
str(sum.profile)
str(deg.sum.p)
sum.profile=sum.profile[rownames(deg.sum.p), ]
str(sum.profile)
sum.profile.mean.as.gs=rowMeans(sum.profile[,gs=="9"])
exp.m.as.gs=exp.m[rownames(deg.sum.p), ]
corr.to.gs9<-apply(exp.m.as.gs, 2, function(x) cor(x, sum.profile.mean.as.gs))
pdf("temp.boxplot.gs.TCGA.corr.sumPrf.pdf", width=8, height=8, paper='special')
# sig=setdiff( rownames(deg.9_7), rownames(deg.7_6) )
for_temp_color=corr.to.gs9
temp.df=data.frame(gs.sig = for_temp_color, gs = gleason)
res <- wilcox.test(temp.df$gs.sig[temp.df$gs=="GS 6, 7"], temp.df$gs.sig[temp.df$gs=="GS 8, 9, 10"])
require(canprot)
ef.size<- CLES(temp.df$gs.sig[temp.df$gs=="GS 8, 9, 10"], temp.df$gs.sig[temp.df$gs=="GS 6, 7"])
boxplot(gs.sig~gs, data=temp.df, ylab="mean expression of up-regulated genes identified in single-cell data",main =paste("p-value: ", formatC(res$p.value, digits=3), ", Wilcoxon rank sum test\n", "CLES: ", formatC(ef.size, digits=3) , sep=""))
dev.off()











