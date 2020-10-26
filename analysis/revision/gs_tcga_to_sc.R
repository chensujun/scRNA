library(BoutrosLab.plotting.general);
library(readxl);
library(Seurat);
library(matrixStats);
library(reshape);
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/GS')
rpkm.t <- readRDS('~/chensj/scRNA/primary/data/2019-01-05_tcga_fpkm.rds');
clin.t <- data.frame(readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/data/TCGA.prad.clinical.rds'));
colnames(rpkm.t) <- gsub('.[0-9]+$','', colnames(rpkm.t));
rownames(clin.t) <- gsub('-', '.', rownames(clin.t));
clin.sub <- clin.t[clin.t$gleasonscore%in%c(10, 6), ];
rpkm.sub <- rpkm.t[, colnames(rpkm.t)%in%rownames(clin.sub)];
clin.sub <- clin.sub[colnames(rpkm.sub), ];
group <- clin.sub$gleasonscore;
pval <- data.frame(pval = apply(rpkm.sub, 1, function(x) wilcox.test(x~group)$p.value));
pval <- na.omit(pval);
pval$padj <- p.adjust(pval$pval);
saveRDS(pval, paste0(Sys.Date(), '_pval_GS10v6.rds'));
pval <- pval[pval$pval<0.01, ]
rpkm.scale <- data.frame(t(scale(t(rpkm.sub[rownames(pval), ]))))
pval$gs6 <- rowMeans(rpkm.scale[rownames(pval), group==6]);
pval$gs10 <- rowMeans(rpkm.scale[rownames(pval), group==10])
###
pval <- pval[rownames(pval)%in%rownames(seurat.all@data), ];
mycount <- seurat.all@data[rownames(pval), ];
mycount.scale <- data.frame(t(scale(t(mycount))));
mycount.scale <- na.omit(mycount.scale)
pval <- pval[rownames(mycount.scale), ];
mycor <- data.frame(cor = apply(mycount.scale, 2, function(x) cor(x, pval$gs6)));
mycor$samp <- gsub('SL', '', seurat.all@meta.data[gsub('\\.', '-', rownames(mycor)), ]$orig.ident)
annot <- data.frame(read_excel('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F1/Table_S1_SampleInfo.xlsx'));
mycor$gs <- annot[match(mycor$samp, annot$sample), ]$pGleason_score;
mycor$group <- ifelse(mycor$gs%in%c('4+5', '5+5'), 'G1', 'G2');
###
mycount.samp <- data.frame(sapply(unique(seurat.all@meta.data$orig.ident), function(x) rowMeans(mycount[, seurat.all@meta.data$orig.ident==x])))
mycount.samp <- data.frame(t(scale(t(mycount.samp))));
mycount.samp <- mycount.samp[rownames(pval), ];
mycor.samp <- data.frame(cor = apply(mycount.samp, 2, function(x) cor(x, pval$gs6)));
mycor.samp$gs <- annot[match(gsub('SL', '', rownames(mycor.samp)), annot$sample), ]$pGleason_score;
mycor.samp$group <- ifelse(mycor.samp$gs%in%c('4+5', '5+5'), 'G1', 'G2');
###### calc sig from averaged pseudo bulk from scRNA-seq
mycount.samp <- data.frame(sapply(unique(seurat.all@meta.data$orig.ident), function(x) rowMeans(seurat.all@data[, seurat.all@meta.data$orig.ident==x])))
all(colnames(mycount.samp)==rownames(mycor.samp));
pval.sc <- data.frame(pval = apply(mycount.samp, 1, function(x) wilcox.test(x~mycor.samp$group)$p.value));
pval.sc <- na.omit(pval.sc)
pval.sc$padj <- p.adjust(pval.sc$pval);
saveRDS(pval.sc, file = paste0(Sys.Date(), '_pval_GS9v7_sc.rds'));
pval.sc <- pval.sc[pval.sc$pval<0.01, ];
pval.sc$gs9 <- rowMeans(t(scale(t(mycount.samp)))[rownames(pval.sc), mycor.samp$group=='G1']);
rpkm.test <- rpkm.t[rownames(rpkm.t)%in%rownames(pval.sc), ];
pval.sc <- pval.sc[rownames(rpkm.test), ];
rpkm.test <- na.omit(data.frame(t(scale(t(rpkm.test)))));
pval.sc <- pval.sc[rownames(rpkm.test), ]
pval.t <- data.frame(cor = apply(rpkm.test, 2, function(x) cor(x, pval.sc$gs9)));
pval.t$gs <- clin.t[rownames(pval.t), ]$gleasonscore;
pval.t$group <- as.factor(ifelse(pval.t$gs%in%c(6, 7), '6+7', '8+9+10'));
saveRDS(pval.t, file =paste0(Sys.Date(), '_cor_sigSC.rds'));
mypval <- wilcox.test(pval.t$cor~pval.t$group)$p.value
create.boxplot(
	formula = cor~group,
	data = pval.t, 
	add.stripplot = TRUE,
	abline.h = 0,
	abline.lty = 2,
	xlab.label = 'TCGA GS group',
	ylab.label = 'Correlation with scRNA signature',
	add.text = TRUE,
	text.x = 1.5,
	text.y = 0.6,
	text.labels = scientific.notation(mypval),
	filename = generate.filename('correlation', 'sigSC', 'pdf'),
	style = 'Nature'
	);
####
pval <- readRDS('2020-03-05_pval_GS10v6.rds');
pval <- pval[rownames(pval)%in%rownames(seurat.all@data), ];
mycount <- seurat.all@data[rownames(pval), ];
mycount.samp <- data.frame(sapply(unique(seurat.all@meta.data$orig.ident), function(x) rowMeans(seurat.all@data[, seurat.all@meta.data$orig.ident==x])))
myraw <- seurat.all@raw.data[rownames(pval), ];
myraw.samp <- data.frame(sapply(unique(seurat.all@meta.data$orig.ident), function(x) rowSums(myraw[, seurat.all@meta.data$orig.ident==x])))
myraw.samp[, 1:13] <- sapply(seq(ncol(myraw.samp)), function(x) 1000000*myraw.samp[, x]/sum(myraw.samp[, x]));
myraw.samp <- log2(myraw.samp + 1);
saveRDS(myraw.samp, file = generate.filename('sum_profile', 'all', 'rds'))
#myraw.samp <- data.frame(t(scale(t(myraw.samp))));
#myraw.samp <- na.omit(myraw.samp);
annot <- data.frame(read_excel('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F1/Table_S1_SampleInfo.xlsx'));
group <- factor(ifelse(annot[match(gsub('SL', '', colnames(myraw.samp)), annot$sample), ]$pGleason_score%in%c('4+5', '5+5'), 'G1', 'G2'))
pval.sc <- data.frame(pval = apply(myraw.samp, 1, function(x) wilcox.test(x~group)$p.value));
pval.sc <- na.omit(pval.sc)
pval.sc$padj <- p.adjust(pval.sc$pval);
saveRDS(pval.sc, file = paste0(Sys.Date(), '_pval_GS9v7_sc_raw.rds'));
pval.sc <- pval.sc[pval.sc$pval<0.01, ];
pval.sc$gs9 <- rowMeans(t(scale(t(myraw.samp)))[rownames(pval.sc), group=='G1']);
#rpkm.test <- rpkm.t[rownames(rpkm.t)%in%rownames(pval.sc), ];
rpkm.test <- log2(rpkm.t[rownames(rpkm.t)%in%rownames(pval.sc), ] + 1);
pval.sc <- pval.sc[rownames(rpkm.test), ];
rpkm.test <- na.omit(data.frame(t(scale(t(rpkm.test)))));
pval.sc <- pval.sc[rownames(rpkm.test), ]
pval.t <- data.frame(cor = apply(rpkm.test, 2, function(x) cor(x, pval.sc$gs9)));
pval.t$gs <- clin.t[rownames(pval.t), ]$gleasonscore;
pval.t$group <- as.factor(ifelse(pval.t$gs%in%c(6, 7), '6+7', '8+9+10'));
saveRDS(pval.t, file =paste0(Sys.Date(), '_cor_sigSC.rds'));
mypval <- wilcox.test(pval.t$cor~pval.t$group)$p.value
create.boxplot(
	formula = cor~group,
	data = pval.t, 
	add.stripplot = TRUE,
	abline.h = 0,
	abline.lty = 2,
	xlab.label = 'TCGA GS group',
	ylab.label = 'Correlation with scRNA signature',
	add.text = TRUE,
	text.x = 1.5,
	text.y = 0.4,
	text.labels = scientific.notation(mypval),
	filename = generate.filename('correlation', 'sigSC_raw', 'pdf'),
	style = 'Nature'
	);
###
rv <- rowVars(as.matrix(rpkm.t));
ntop <- 500;
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))];
rpkm.test <- log2(rpkm.t[select, ]);
rpkm.test <- na.omit(data.frame(t(scale(t(rpkm.test)))));
rpkm.avg <- data.frame(all = rowMeans(rpkm.test), 
	gs9 = rowMeans(rpkm.test[, colnames(rpkm.test)%in%rownames(clin.t[clin.t$gleasonscore%in%c(8, 9,10), ])]));

pval.avg <- data.frame(cor = apply(rpkm.test, 2, function(x) cor(x, rpkm.avg$all)));
pval.avg$cor9 <- apply(rpkm.test, 2, function(x) cor(x, rpkm.avg$gs9))
pval.avg$gs <- clin.t[rownames(pval.t), ]$gleasonscore;
pval.avg$group <- as.factor(ifelse(pval.t$gs%in%c(6, 7), '6+7', '8+9+10'));
wilcox.test(pval.avg$cor~pval.avg$group)$p.value;
####
#### t test for differential genes
annot <- data.frame(read_excel('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F1/Table_S1_SampleInfo.xlsx'));
group <- factor(ifelse(annot[match(gsub('SL', '', colnames(myraw.samp)), annot$sample), ]$pGleason_score%in%c('4+5', '5+5'), 'G1', 'G2'))
pval.sc <- data.frame(pval = apply(myraw.samp, 1, function(x) t.test(x~group)$p.value));
pval.sc <- na.omit(pval.sc)
pval.sc$padj <- p.adjust(pval.sc$pval);
saveRDS(pval.sc, file = paste0(Sys.Date(), '_pval_GS9v7_sc_raw_ttest.rds'));
pval.sc <- pval.sc[pval.sc$pval<0.01, ];
pval.sc$gs9 <- rowMeans(t(scale(t(myraw.samp)))[rownames(pval.sc), group=='G1']);
#rpkm.test <- rpkm.t[rownames(rpkm.t)%in%rownames(pval.sc), ];
rpkm.test <- log2(rpkm.t[rownames(rpkm.t)%in%rownames(pval.sc), ] + 1);
pval.sc <- pval.sc[rownames(rpkm.test), ];
rpkm.test <- na.omit(data.frame(t(scale(t(rpkm.test)))));
pval.sc <- pval.sc[rownames(rpkm.test), ]
pval.t <- data.frame(cor = apply(rpkm.test, 2, function(x) cor(x, pval.sc$gs9)));
pval.t$gs <- clin.t[rownames(pval.t), ]$gleasonscore;
pval.t$group <- as.factor(ifelse(pval.t$gs%in%c(6, 7), '6+7', '8+9+10'));
saveRDS(pval.t, file =paste0(Sys.Date(), '_cor_sigSC.rds'));
mypval <- wilcox.test(pval.t$cor~pval.t$group)$p.value
create.boxplot(
	formula = cor~group,
	data = pval.t, 
	add.stripplot = TRUE,
	abline.h = 0,
	abline.lty = 2,
	xlab.label = 'TCGA GS group',
	ylab.label = 'Correlation with scRNA signature',
	add.text = TRUE,
	text.x = 1.5,
	text.y = 0.4,
	text.labels = scientific.notation(mypval),
	filename = generate.filename('correlation', 'sigSC_raw', 'pdf'),
	style = 'Nature'
	);
#### use most variable genes for reference profile
#### normalize correlation, similarity etc
#### microdessection score in scRNAseq data
#### use epithelia cells only to generate bulk profile
#### use changhai data?
#### table for epi% for each sample unresolved cells
#### check pathology for mal cells with Bobby
#### summarize % mal cells and compare with pathology
#### basal cell % table 
#### change circle size to abs cell number 
#### change % to number in rebuttal letter
#### remove RF3A 
#### add before correction KLK3 expression, calculate pvalue, no sigificant change 
#### RF5 add trend line >> change RF5 
#### RF5 label signifiant correlatted ones, add prad
#### Move figure S3K to rebutal letter 

#### epithelial cluster, malg%/non malignant distribution
#### epithelial cluster, colour by malignancy
#### malig vs non epithelial cells diff gene *
#### cellularity for sample 171
#### set 0 to 0 F3   change to barplot for F3 
#### #/% KLK3+ tcells, contingency table 
#### t cell function and marker gene expression correlation per patient 
#### compare within effector cells for KLK3 and function Correlation, check also treg activity
#### plot score distribution and color according to cnv results
########
#### color according to linear cnv malignancy for F1A
#### add % malig in F1D
########
#### test cell cycle/basal/aEC in PNAS paper, felix paper
#### 1. cnv in RF
#### 2. DPP4 to supp
#### 3. GS to supp
#### 4. change figure name, supp to extended data
#### check more data sample polyA
#### use the general cell type for GS, generate EF to supp
#### remove RF1B, change color of tsne
#### F3D, plot only for tumors, re-arrange fig order
name <- 'tcell_KLK3'
genes.marker <- c('KLK3');
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
pdf(generate.filename('plotgene', paste0(name, 0), 'pdf'), width = 5, height = 5);
FeaturePlot(seurat.t, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 4, nCol = 1, min.cutoff = 0,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
