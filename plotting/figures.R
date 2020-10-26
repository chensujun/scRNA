library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/Finalize');
### CNV, SF1f
cnv.gene <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/multi/2020-05-22_cnv_gene.rds');
ncnv <- data.frame(high = sapply(cnv.gene, function(x) length(x[x>0.04])));
ncnv$all <- sapply(cnv.gene, function(x) length(x));
ncnv$samp <- names(cnv.gene);
ncnv$pct <- 100*ncnv$high/ncnv$all;
ncnv$group <- factor(ifelse(grepl('JD', ncnv$samp), 'pca', 'mel'));
ncnv$h1 <- sapply(cnv.gene, function(x) length(x[x>0.1]));
ncnv$pct1 <- 100*ncnv$h1/ncnv$all;
pval <- scientific.notation(wilcox.test(ncnv$pct1~ncnv$group)$p.value)
create.boxplot(
	formula = pct1~group,
	data = ncnv,
	add.stripplot = TRUE,
	xaxis.lab = c('Mel', 'PCa'),
	xlab.label = "",
	ylab.label = '% of gene with strong CNV',
	add.text = TRUE,
	text.x = 1.5,
	text.y = 30.5,
	text.label = pval,
	filename = generate.filename('pct', 'altered01', 'pdf'),
	style = 'Nature',
	width = 4,
	height = 4
    );

### CAF, F5a
iseurat <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/fibroblast/objects/13_sample_fibroblast_20190718/GraphClust.seuset.rds');
mywidth <- 12;
myheight <- 12;
myres <- 300;
name <- 'caf';
mycol <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/data/2019-07-25_mycol_type.rds');
names(mycol) <- seq(0, 6);
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));
fontsize1 <- theme(axis.text=element_text(size=24), axis.title=element_text(size=24), legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1));
pdf(generate.filename('plottsne', paste0(name, '_cluster'), 'pdf'), width = mywidth, height = myheight)
DimPlot(iseurat, reduction.use = 'tsne', pt.size = 3, vector.friendly = TRUE, cols.use = mycol)+ fontsize +
	labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
## CAF, F5b
genes.marker <- c('FAP', 'S100A4', 'SPARC', 'ACTA2', 'PDGFRA', 'PDGFRB', 'CAV1', 'VIM');
pdf(generate.filename('plotgene', name, 'pdf'), width = 12, height = 12);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 4, nCol = 3, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
pdf(generate.filename('plotgene_leg', name, 'pdf'), width = 12, height = 12);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 4, nCol = 3, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE, no.legend = FALSE);
dev.off();

###CAF, F5c
to.plot <- read.csv('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/Figures/SF1/2020-02-27_plotdata_F_S1_reviewer.csv', row.names = 1);
create.scatterplot(
	formula = pCAF.2~emt.h1,
	data = to.plot,
	xlab.label = 'EMT score in epithelial cells',
	ylab.label = 'ACTA2+ CAF percentage',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	#xlimits = c(0, 2.5),
	xaxis.cex = 2,
	yaxis.cex = 2,
	filename = generate.filename('fib', 'emt', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$pCAF.2,
                y = to.plot$emt.h1,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 2
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    ),
    width = 5,
    height = 5,
);
####TAM, F3a
iseurat <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/macrophage/objects/Macrophage/1.GraphCluster/Manual.seuset.rds')
name <- 'tam'
pdf(generate.filename('plottsne', paste0(name, '_cluster'), 'pdf'), width = mywidth, height = myheight)
DimPlot(iseurat, reduction.use = 'tsne', pt.size = 3, vector.friendly = TRUE, cols.use = mycol)+ fontsize +
	labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();
##TAM, F3c,d
library(readxl);
sig.m1 <- data.frame(read_excel('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/macrophage/objects/1-s2.0-S0092867418307232-mmc4.xlsx', sheet = 2));
sig.m2s <- data.frame(read_excel('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/macrophage/objects/1-s2.0-S0092867418307232-mmc4.xlsx', sheet = 3));
gene.m1 <- gsub(' .*', '', sig.m1$Gene);
gene.m2s <- gsub(' .*', '', sig.m2s$Gene);
gene.m1 <- c(intersect(gene.m1, rownames(iseurat@data)), 'NOS2', 'NOS1', 'NOS3', 'FCGR1A', 'FCGR1B', 'IL12A', 'TNF');
gene.m2 <- c(intersect(gene.m2s, rownames(iseurat@data)), 'ARG1', 'ARG2', 'FCGR2A', 'FCGR2B', 'FCGR2C', 'FCER2', 'PDCD1LG2', 'CD274', 'MRC1', 'IL1RN', 
	'IL1R2', 'CTSA', 'CTSC', 'WNT7B', 'FASLG');
iseurat@meta.data$TAM_m1 <- colMeans(iseurat@data[rownames(iseurat@data)%in%gene.m1, ]);
iseurat@meta.data$TAM_m2 <- colMeans(iseurat@data[rownames(iseurat@data)%in%gene.m2s, ]);
###
pdf(generate.filename('plotviolin', paste0(name, '_geneset'), 'pdf'), width = 8, height = 10);
VlnPlot(iseurat, features.plot = gene.m1, x.lab.rot = TRUE, point.size.use = 0.1, nCol = 4);
dev.off();
###
to.plot <- iseurat@meta.data;
create.scatterplot(
	TAM_m1~TAM_m2,
	to.plot,
	groups = to.plot$subtype,
	col = default.colours(length(unique(to.plot$subtype))),
	cex = 0.5,
	add.xyline = TRUE,
	xyline.lty = 2,
	file = generate.filename('sig_tam', name, 'pdf'),
	style = 'Nature',
	legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$TAM_m1,
                y = to.plot$TAM_m2,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 2
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
	);




