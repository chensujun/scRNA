###using guassian model to select predict malignancy
library(BoutrosLab.plotting.general);
library(mclust);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/multi/cluster');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
annot <- seurat.all@meta.data;
rownames(annot) <- gsub('-', '.', rownames(annot));
to.plot <- readRDS('../2020-05-26_cnv_score_pri_0518.rds');
to.plot$type <- annot[rownames(to.plot), ]$type;
mydata <- to.plot[, 1:2];
pdf(generate.filename('distribution', 'cnv', 'pdf'));
plot(density(to.plot$score));
plot(density(to.plot$corr));
dev.off();
dens <- densityMclust(mydata);
saveRDS(dens, file = generate.filename('distribution', 'cnv', 'rds'));

ddat <- mydata$score
h1 <- hist(ddat);
x <- seq(min(ddat)-diff(range(ddat))/10, max(ddat)+diff(range(ddat))/10, length = 200);
cdens <- predict(dens, mydata, what = 'cdens');
cdens <- t(apply(cdens, 1, function(d) d*dens$parameters$pro));
col <- adjustcolor(mclust.options("classPlotColors")[1:2], alpha = 0.3)
plot(h1, xlab = "CNV Score", freq = FALSE, main = "", border = FALSE, col = col[1],
                   xlim = range(x), )
plot(h2, add = TRUE, freq = FALSE, border = FALSE, col = col[2])
matplot(x, cdens, type = "l", lwd = 1, add = TRUE, lty = 1:3, col = 1);
box();

out.s <- cutoff::em(ddat, 'normal', 'normal');
hist(ddat, 100, F, xlim = c(0, 0.1))
lines(out.s,lwd=1.5,col="red")

pdf(generate.filename('distribution', 'cnv', 'pdf'));
plotDensityMclustd(ddat);
plotDensityMclust1(ddat, data = NULL, hist.col = 'lightgrey', hist.border = 'white', breaks = 'Sturges');
dev.off();

pdf(generate.filename('pairwise', 'cluster', 'pdf'));
clPairs(mydata, to.plot$col);
clPairs(mydata, to.plot$type);
dev.off();
set.seed(100);
BIC <- mclustBIC(mydata);
saveRDS(BIC, file = generate.filename('mclust', 'BIC', 'rds'));
pdf(generate.filename('mclust', 'BIC', 'pdf'));
plot(BIC);
dev.off();
###
set.seed(0);
mod1 <- Mclust(mydata[, 1:2], x = BIC, G = 2);
boot1 <- MclustBootstrap(mod1, nboot = 1000, type = 'bs');
summary(boot1, what = 'ci');
pdf(generate.filename('mclust', 'ci', 'pdf'));
plot(mod1, what = 'density');
plot(mod1, what = 'uncertainty', ngrid = 200);
dev.off();
###
to.plot$class <- mod1$classification[rownames(mydata)];
to.plot$pval <- apply(mod1$z, 1, min);
mod2 <- Mclust(mydata[, 1:2], x = BIC, G = 3);
to.plot$class3 <- mod2$classification[rownames(to.plot)];
to.plot$pval3 <- mod2$z[, 1];
###
create.scatterplot(
	formula = corr~score,
	data = to.plot,
	col = to.plot$class,
	cex = 0.5,
	xaxis.cex = 2,
	yaxis.cex = 2,
	xlimits = c(0, 0.3),
	xlab.label = 'CNV Score', 
	ylab.label = 'CNV Correlation',
	filename = generate.filename('mclust', 'class', 'pdf'),
	style = 'Nature',
	width = 15,
	height = 15
	);
####
to.plot$group <- ifelse(to.plot$type%in%c('Basal/intermediate', 'Luminal'), 'epi', 'non');
#mod3 <- MclustDA(mydata[, 1:2], to.plot$group, modelNames = 'VVV');
###
to.plot$uncertainty <- mod1$uncertainty;
table(ifelse(to.plot$uncertainty>0.05, 3, to.plot$class), to.plot$type);
####
to.plot$class <- ifelse(to.plot$uncertainty>0.05, 3, to.plot$class)
create.scatterplot(
	formula = corr~score,
	data = to.plot,
	col = to.plot$class,
	cex = 0.5,
	xaxis.cex = 2,
	yaxis.cex = 2,
	xlimits = c(0, 0.3),
	xlab.label = 'CNV Score', 
	ylab.label = 'CNV Correlation',
	filename = generate.filename('mclust', 'class', 'pdf'),
	style = 'Nature',
	width = 15,
	height = 15
	);
#### summarize results
annot <- seurat.all@meta.data;
rownames(annot) <- gsub('-', '.', rownames(annot));
to.plot <- readRDS('../2020-05-26_cnv_score_pri_0518.rds');
to.plot$type <- annot[rownames(to.plot), ]$type;
to.plot$samp <- annot[rownames(to.plot), ]$orig.ident;
to.plot$type <- ifelse(annot[match(rownames(to.plot), rownames(annot)), ]$type%in%c('Basal/intermediate', 'Luminal'), 'epithelia', 'non');
mydata <- to.plot[, 1:2];
BIC <- readRDS('2020-05-28_mclust_BIC.rds');
set.seed(0);
mod1 <- Mclust(mydata[, 1:2], x = BIC, G = 2);
saveRDS(mod1, file = generate.filename('mclust', 'mod1', 'rds'));
to.plot$class <- mod1$classification[rownames(mydata)];
to.plot$uncertainty <- mod1$uncertainty;
to.plot$class <- ifelse(to.plot$uncertainty>0.05, 3, to.plot$class);
saveRDS(to.plot, file = generate.filename('class', 'gmm', 'rds'));
###
to.plot <- readRDS('2020-06-23_class_gmm.rds');
###
for(name in unique(to.plot$samp)){
	plot_distr_class(to.plot[to.plot$samp==name, ], paste0('multiB3_mclust', name), gsub('JD1800|SL', '', name), TRUE)
};
###
to.plot$type <- annot[rownames(to.plot), ]$type;
to.plot$type2 <- ifelse(to.plot$type%in%c('Basal/intermediate', 'Luminal'), 'epi', 'non');
pct <- data.frame(name = unique(to.plot$samp), pct = NA, pct.b = NA);

for(name in unique(to.plot$samp)){
	idat  <- to.plot[to.plot$samp==name, ];
	print(name)
	print(nrow(idat[idat$class==2, ])/nrow(idat[idat$type2=='epi', ]));
	pct[pct$name==name, ]$pct <- round(nrow(idat[idat$class==2, ])/nrow(idat[idat$type2=='epi', ])*100, 2)
	pct[pct$name==name, ]$pct.b <- round(nrow(idat[idat$class==2&idat$type=='Basal/intermediate', ])/nrow(idat[idat$type=='Basal/intermediate', ])*100, 2)
};
##### summarize results
mal.gmm <- readRDS('malignant.epithelial.meta.data.rds');
seurat.epi <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/epitumour/objects/2019-08-27_sample13_epi.rds');
mal.gmm$cluster.ori <- seurat.epi@meta.data[rownames(mal.gmm), ]$fig.cluster;
epi <- seurat.epi@meta.data;
rownames(epi) <- gsub('-', '.', rownames(epi));
to.plot <- to.plot[rownames(epi), ]
to.plot$cluster <- epi[rownames(to.plot), ]$fig.cluster;
write.csv(table(to.plot$cluster, to.plot$col), 'epi_cluster_linear.csv');
###
seurat.epi@meta.data$class <- to.plot[gsub('-', '.', rownames(seurat.epi@meta.data)), ]$class;
name <- 'epi'
pdf(generate.filename('plottsne', paste0(name, '_class'), 'pdf'));
DimPlot(seurat.epi, reduction.use = 'tsne', pt.size = 1, group.by = 'class', do.label = TRUE, vector.friendly = TRUE, label.size = 0, no.legend = FALSE, no.axes = FALSE) + 
labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
####
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/Figures/F1');
source('~/svn/singleCell/myfunctions/plot_distribution.R');
name <- 'all';
seurat.all@meta.data$type <- gsub('Macrophage', 'Monolytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
seurat.all <- SetAllIdent(seurat.all, id = 'type');
####
mywidth <- 16;
myheight <- 12;
myres <- 300;
#### F1A
to.plot <- readRDS('../2020-05-26_cnv_score_pri_0518.rds');
mytype <- c('malignant', 'non malignant', 'unresolved');
names(mytype) <- c('red', 'blue', 'black');
to.plot$class <- mytype[to.plot$col];
seurat.all@meta.data$class <- to.plot[gsub('-', '.', rownames(seurat.all@meta.data)), ]$class;
tiff(generate.filename('plottsne', paste0(name, '_cnv'), 'tiff'), width = mywidth, height = mywidth, res = myres, units = 'in')
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'class', do.label = FALSE, label.size = 24, 
	no.legend = TRUE, no.axes = TRUE);
dev.off();
###
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));

name <- 'all';
pdf(generate.filename('plottsne', paste0(name, '_cnv'), 'pdf'), width = mywidth + 3, height = mywidth)
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'class', do.label = FALSE, label.size = 24, 
	no.legend = FALSE, no.axes = FALSE, vector.friendly = TRUE, cols.use = names(mytype)) + 
	fontsize + labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();


tiff(generate.filename('plottsne', paste0(name, '_cnv'), 'tiff'), width = mywidth + 1, height = mywidth, res = myres, units = 'in')
DimPlot(seurat.all, reduction.use = 'tsne', pt.size = 1, group.by = 'class', do.label = FALSE, label.size = 24, 
	no.legend = TRUE, no.axes = FALSE, vector.friendly = TRUE, cols.use = names(mytype)) + 
	fontsize + labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
