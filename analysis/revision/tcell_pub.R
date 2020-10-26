library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape);
library(pheatmap);
library(plyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/tcell');
seurat.b <- readRDS('Breast/GraphClust.seuset-breast.rds');
name <- 'pub_breast';
mywidth <- 16;
myheight <- 12;
myres <- 300;
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));

pdf(generate.filename('plottsne', paste0(name, '_cluster'), 'pdf'), width = mywidth, height = myheight);
DimPlot(seurat.b, reduction.use = 'tsne', pt.size = 1, group.by = 'res.0.8', do.label = TRUE, vector.friendly = TRUE, label.size = 24, no.legend = TRUE, no.axes = FALSE) + 
fontsize + labs(x = 'Dimension 1', y = 'Dimension 2');
dev.off();

tiff(generate.filename('plottsne', paste0(name, '_cluster'), 'tiff'), width = mywidth, height = mywidth, res = myres, units = 'in')
DimPlot(seurat.b, reduction.use = 'tsne', pt.size = 1, group.by = 'type', do.label = TRUE, label.size = 24, no.legend = TRUE, no.axes = FALSE) + 
fontsize + labs(x = 'Dimension 1', y = 'Dimension 2')
dev.off();
####
mytype <- 'Breast';
seurat.b <- readRDS('Breast/GraphClust.seuset-breast.rds');
flist <- list.files(paste0(mytype, '/QuSAGE_GeneSet_result/GeneSet41_Gene_Function-change/'), 'txt', full.names = TRUE);
tmp <- read.table(flist[1], header = TRUE);
read_n_sort <- function(x, y){
	ifile <- read.table(x, header = TRUE);
	ifile <- ifile[match(y$GeneSet.name, ifile$GeneSet.name), ];
	return(ifile);
}
to.plot <- do.call(cbind, sapply(flist, function(x) read_n_sort(x, tmp)[, 2, drop = FALSE]));
colnames(to.plot) <- gsub('.*GeneSetInfo_|_GeneSet41.*', '', colnames(to.plot));
rownames(to.plot) <- tmp$GeneSet.name;
to.plot <- to.plot[grep('CD8|Pro|Treg|Cyto', rownames(to.plot)), ];
to.plot <- data.frame(t(to.plot));
myexp <- data.frame(t(seurat.b@data[c('KRT18', 'ITGA6', 'HIF1A', 'CEACAM1', 'CEACAM4'), ]));
myexp$cluster <- iseurat@meta.data$res.0.8;
myexp <- ddply(myexp, 'cluster', numcolwise(mean));
myexp$cluster <- paste0('Cluster_', '', myexp$cluster);
myexp <- myexp[match(rownames(to.plot), myexp$cluster), ];
to.plot[, 8:12] <- myexp[, -1];
create.scatterplot(
    formula = CD8_T_Cell_Activation_gene~KRT18,
    data = to.plot,
    xlab.label = 'mean KRT18',
    ylab.label = 'CD8_T_Cell_Activation_signature',
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    #xlimits = c(0, 2.5),
    filename = generate.filename('KLK3_CD8activation', 'BRCAtcell_KRT18', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$KRT18,
                y = to.plot$CD8_T_Cell_Activation_gene,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);

create.scatterplot(
	formula = CD8_T_Cell_Activation_gene~CEACAM4,
	data = to.plot,
	xlab.label = 'mean CEACAM4',
	ylab.label = 'CD8_T_Cell_Activation_signature',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	#xlimits = c(0, 2.5),
	filename = generate.filename('KLK3_CD8activation', 'BRCAtcell_CEACAM4', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$CEACAM4,
                y = to.plot$CD8_T_Cell_Activation_gene,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);

create.scatterplot(
	formula = CD8_T_Cell_Activation_gene~CEACAM1,
	data = to.plot,
	xlab.label = 'mean CEACAM1',
	ylab.label = 'CD8_T_Cell_Activation_signature',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	#xlimits = c(0, 2.5),
	filename = generate.filename('KLK3_CD8activation', 'BRCAtcell_CEACAM1', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$CEACAM1,
                y = to.plot$CD8_T_Cell_Activation_gene,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);

####
seurat.h <- readRDS('HCC/GraphClust.seuset-HCC.rds');
mytype <- 'HCC';
iseurat <- seurat.h;
flist <- list.files(paste0(mytype, '/QuSAGE_GeneSet_result/GeneSet41_Gene_Function-change/'), 'txt', full.names = TRUE);
tmp <- read.table(flist[1], header = TRUE);
read_n_sort <- function(x, y){
    ifile <- read.table(x, header = TRUE);
    ifile <- ifile[match(y$GeneSet.name, ifile$GeneSet.name), ];
    return(ifile);
}
to.plot <- do.call(cbind, sapply(flist, function(x) read_n_sort(x, tmp)[, 2, drop = FALSE]));
colnames(to.plot) <- gsub('.*GeneSetInfo_|_GeneSet41.*', '', colnames(to.plot));
rownames(to.plot) <- tmp$GeneSet.name;
to.plot <- to.plot[grep('CD8|Pro|Treg|Cyto', rownames(to.plot)), ];
to.plot <- data.frame(t(to.plot));
myexp <- data.frame(iseurat[['RNA']]@data);
myexp <- data.frame(t(myexp[c('ENO2', 'GOLM1'), ]))
myexp$cluster <- iseurat@meta.data$RNA_snn_res.0.8;
myexp <- ddply(myexp, 'cluster', numcolwise(mean));
myexp$cluster <- paste0('Cluster_', '', myexp$cluster);
myexp <- myexp[match(rownames(to.plot), myexp$cluster), ];
to.plot[, 8:10] <- myexp[, -1];
mygene <- 'ENO2';
saveRDS(myexp, file = paste0(Sys.Date(), '_myexp_hcc.rds'));
myexp <- readRDS('2020-03-30_myexp_hcc.rds');
to.plot[, 8:10] <- myexp[match(rownames(to.plot), myexp$cluster), -1];
create.scatterplot(
    formula = CD8_T_Cell_Activation_gene~ENO2,
    data = to.plot,
    xlab.label = 'mean KRT18',
    ylab.label = 'CD8_T_Cell_Activation_signature',
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    #xlimits = c(0, 2.5),
    filename = generate.filename('KLK3_CD8activation_tcell', paste0(mytype, '_', mygene), 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot[, mygene],
                y = to.plot$CD8_T_Cell_Activation_gene,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);

mygene <- 'GOLM1';
create.scatterplot(
    formula = CD8_T_Cell_Activation_gene~GOLM1,
    data = to.plot,
    xlab.label = 'mean CEACAM4',
    ylab.label = 'CD8_T_Cell_Activation_signature',
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    #xlimits = c(0, 2.5),
    filename = generate.filename('KLK3_CD8activation_tcell', paste0(mytype, '_', mygene), 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot[, mygene],
                y = to.plot$CD8_T_Cell_Activation_gene,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    )
);
#####
mytype <- 'CRC';
myexp <- read.table('data/CRC_RDS.ClusterMean.txt', header = TRUE, row.names = 1);
mypath <- read.table('data/CRC_Human_41GeneSets_190716.heatmap.matrix_CRC.tsv', header = TRUE, row.names = 1, sep = '\t');
colnames(mypath) <- gsub('luster\\.', '', colnames(mypath));
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
mypath <- mypath[grep('CD8|Pro|Treg|Cyto', rownames(mypath)), ];
#to.plot <- data.frame(cbind(t(mypath), t(myexp[c('ENO2', 'GOLM1'), ])));
to.plot <- data.frame(crc = apply(myexp, 1, function(x) cor(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman')));
to.plot.pval <- data.frame(crc = apply(myexp, 1, function(x) cor.test(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman')$p.value));
####
mytype <- 'NSCLC';
myexp <- read.table('data/NSCLC_RDS.ClusterMean.txt', header = TRUE, row.names = 1);
mypath <- read.table('data/NSCLC_QuSAGE_for_heatmap/QuSAGE_for_heatmap/Human_41GeneSets_190716_heatmap_input.txt', header = TRUE, row.names = 1, sep = '\t');
colnames(mypath) <- gsub('luster_', '', colnames(mypath));
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
mypath <- mypath[grep('CD8|Pro|Treg|Cyto', rownames(mypath)), ];
to.plot$nsclc <- apply(myexp[rownames(to.plot), ], 1, function(x) cor(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman'));
to.plot.pval$nsclc <- apply(myexp[rownames(to.plot), ], 1, function(x) cor.test(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman')$p.value);
####
mytype <- 'HCC';
myexp <- readRDS('data/HCC_2020-03-31_clusterMean.rds');
mypath <- read.table('HCC/QuSAGE_for_heatmap/GeneSet41_Gene_Function-change_heatmap_input.txt', header = TRUE, row.names = 1);
colnames(mypath) <- gsub('luster_', '', colnames(mypath));
to.plot$hcc <- apply(myexp[rownames(to.plot), ], 1, function(x) cor(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman'));
to.plot.pval$hcc <- apply(myexp[rownames(to.plot), ], 1, function(x) cor.test(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman')$p.value);
####
#mytype <- 'BRCA';
#myexp <- readRDS('data/BRCA_2020-03-31_clusterMean.rds');
#mypath <- read.table('Breast/QuSAGE_for_heatmap/Human_41GeneSets_190716_heatmap_input.txt', header = TRUE, row.names = 1);
#colnames(mypath) <- gsub('luster_', '', colnames(mypath));
#to.plot$brca <- apply(myexp[rownames(to.plot), ], 1, function(x) cor(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman'));
####
mypath <- readRDS('../../tcell/2019-08-01_qusage_tcell_curated_brca.rds');
mypath <- mypath$results.cor[mypath$results.cor=='CD8_T_Cell_Activation_gene', c(2,4)];
rownames(mypath) <- paste0('C', mypath$Cluster-1);
mypath <- data.frame(t(mypath[, -2, drop = FALSE]));
rownames(mypath) <- 'CD8_T_Cell_Activation_gene';
myexp <- readRDS('data/PRCA_2020-03-31_cluster_Mean.rds');
#seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
#myexp <- seurat.t@data[rownames(to.plot), ];
#myexp <- data.frame(t(myexp));
#myexp$cluster <- seurat.t@meta.data$res.0.8;
#myexp <- ddply(myexp, 'cluster', numcolwise(mean));
#rownames(myexp) <- paste0('C', myexp$cluster);
#myexp <- data.frame(t(myexp[, -1]));
#saveRDS(myexp, paste0('data/PRCA_', generate.filename('cluster', 'Mean', 'rds')));
to.plot$prad <- apply(myexp[rownames(to.plot), ], 1, function(x) cor(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman'));
to.plot.pval$prad <- apply(myexp[rownames(to.plot), ], 1, function(x) cor.test(x, unlist(mypath['CD8_T_Cell_Activation_gene', colnames(myexp)]), method = 'spearman')$p.value);
saveRDS(to.plot, file = generate.filename('correlation', 'all', 'rds'));
saveRDS(to.plot.pval, file = generate.filename('correlation', 'all_pval', 'rds'));
####
pdf(generate.filename('correlaiton', 'all', 'pdf'), height = 4);
pheatmap(t(to.plot), cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, 
	color = colorRampPalette(c('blue', 'white', 'red'))(100));
dev.off();
###

plot_gene_cor <- function(to.plot, gene, mytype, myheight = 6, mywidth = 6){
    to.plot$gene <- to.plot[, gene];
    xlab.max <- max(to.plot$gene);
    xlab.max <- ifelse(xlab.max<0.01, scientific.notation(xlab.max), round(xlab.max, 2))
create.scatterplot(
    formula = CD8_T_Cell_Activation_gene~gene,
    data = to.plot,
    xlab.label = '',
    ylab.label = '',
    xat = c(0, max(to.plot$gene)),
    xaxis.lab = c(0, xlab.max),
    xlimits = c(0, (max(to.plot$gene) + max(to.plot$gene)/10)),
    yat = c(-0.4, 0, 0.4),
    ylimits = c(-0.5, 0.5),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    #xlimits = c(0, 2.5),
    filename = generate.filename('CD8activation', paste0(mytype, '_', gene), 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot[, gene],
                y = to.plot$CD8_T_Cell_Activation_gene,
                label.items = c('spearman','spearman.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.04,
        y = 0.95,
        corner = c(0,1)
        )
    ),
    width = mywidth,
    height = myheight
    );
};

plot_gene_cor <- function(to.plot, gene, mytype, myheight = 6, mywidth = 6, plot.yaxis = FALSE){
    to.plot$gene <- to.plot[, gene];
    xlab.max <- max(to.plot$gene);
    legend.x <- -(xlab.max/10);
    xlab.max <- ifelse(xlab.max<0.01, scientific.notation(xlab.max), round(xlab.max, 2));
    fit <- lm(to.plot$CD8_T_Cell_Activation_gene~to.plot$gene);
    key = get.corr.key(
                    x = to.plot[, gene],
                    y = to.plot$CD8_T_Cell_Activation_gene,
                    label.items = c('spearman','spearman.p'),
                    alpha.background = 0,
                    key.cex = .5
                    );
    mycol <- ifelse(cor.test(to.plot$gene, to.plot$CD8_T_Cell_Activation_gene, method = 'spearman')$p.value<0.05, 'red', 'black')
    pdf(generate.filename('CD8activation', paste0(mytype, '_', gene), 'pdf'), height = myheight, width = mywidth);
    par(mar = rep(3, 4));
    plot(to.plot$gene, to.plot$CD8_T_Cell_Activation_gene, cex = .5, col = mycol,
        pch = 20, xlab = '', ylab = '', axes = FALSE, ylim = c(-0.4, 0.4), xlim = c(0, max(to.plot$gene)));
    axis(1, at = c(0, max(to.plot$gene)), labels = c(0, xlab.max));
    if(plot.yaxis){
    axis(2, las = 2, at = c(-0.4, 0, 0.4), );
    };
    box(bty = 'l');
    abline(fit, col = mycol);
    legend(legend.x, 0.45, legend = key$text$lab, bty = 'n', cex = 0.5);
    dev.off();

};
###
#mycor <- readRDS('2020-03-31_correlation_all.rds');
mycor <- readRDS('2020-06-21_correlation_all.rds');
####
mytype <- 'CRC';
myexp <- read.table('data/CRC_RDS.ClusterMean.txt', header = TRUE, row.names = 1);
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
genes <- rownames(mycor);
to.plot <- data.frame(t(myexp[genes, ]));
to.plot <- to.plot[, rownames(mycor)];
mygenes <- colnames(to.plot[, colSums(to.plot>0.05)>0]);
#mygenes <- c(colnames(to.plot[, colSums(to.plot)>0.5]));

mytype <- 'NSCLC';
myexp <- read.table('data/NSCLC_RDS.ClusterMean.txt', header = TRUE, row.names = 1);
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
genes <- rownames(mycor);
to.plot <- data.frame(t(myexp[genes, ]));
to.plot <- to.plot[, rownames(mycor)];
mygenes <- c(mygenes, colnames(to.plot[, colSums(to.plot>0.05)>0]));
#mygenes <- c(mygenes, colnames(to.plot[, colSums(to.plot)>0.5]));

####
mytype <- 'HCC';
myexp <- readRDS('data/HCC_2020-03-31_clusterMean.rds');
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
genes <- rownames(mycor);
to.plot <- data.frame(t(myexp[genes, ]));
to.plot <- to.plot[, rownames(mycor)];
mygenes <- c(mygenes, colnames(to.plot[, colSums(to.plot>0.05)>0]));
#mygenes <- c(mygenes, colnames(to.plot[, colSums(to.plot)>0.5]));

####
mytype <- 'PRAD';
myexp <- readRDS('data/PRCA_2020-03-31_cluster_Mean.rds');
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
genes <- rownames(mycor);
to.plot <- data.frame(t(myexp[genes, ]));
to.plot <- to.plot[, rownames(mycor)];
mygenes <- c(mygenes, colnames(to.plot[, colSums(to.plot>0.05)>0]));
#mygenes <- c(mygenes, colnames(to.plot[, colSums(to.plot)>0.5]));

####
mytype <- 'CRC';
myexp <- read.table('data/CRC_RDS.ClusterMean.txt', header = TRUE, row.names = 1);
mypath <- read.table('data/CRC_Human_41GeneSets_190716.heatmap.matrix_CRC.tsv', header = TRUE, row.names = 1, sep = '\t');
colnames(mypath) <- gsub('luster\\.', '', colnames(mypath));
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
mypath <- mypath[grep('CD8|Pro|Treg|Cyto', rownames(mypath)), ];
#genes <- rownames(mycor[rowSums(mycor[, 1:3]<0)>0, ]);
genes <- rownames(mycor);
genes <- unique(mygenes);
genes <- genes[order(genes)];
to.plot <- data.frame(cbind(t(mypath), t(myexp[genes, ])));
#plot_gene_cor(to.plot, 'ENO2', mytype);
#plot_gene_cor(to.plot, 'GOLM1', mytype);
plot_gene_cor(to.plot, genes[1], mytype, mywidth = 3, myheight = 3, plot.yaxis = TRUE);  
for(gene in genes[-1])
{
  plot_gene_cor(to.plot, gene, mytype, mywidth = 3, myheight = 3);  
};
####
mytype <- 'NSCLC';
myexp <- read.table('data/NSCLC_RDS.ClusterMean.txt', header = TRUE, row.names = 1);
mypath <- read.table('data/NSCLC_QuSAGE_for_heatmap/QuSAGE_for_heatmap/Human_41GeneSets_190716_heatmap_input.txt', header = TRUE, row.names = 1, sep = '\t');
colnames(mypath) <- gsub('luster_', '', colnames(mypath));
colnames(myexp) <- gsub('X', 'C', colnames(myexp));
mypath <- mypath[grep('CD8|Pro|Treg|Cyto', rownames(mypath)), ];
genes <- rownames(mycor);
genes <- unique(mygenes)
genes <- genes[order(genes)];
to.plot <- data.frame(cbind(t(mypath), t(myexp[genes, ])));
#plot_gene_cor(to.plot, 'ENO2', mytype);
#plot_gene_cor(to.plot, 'GOLM1', mytype);
plot_gene_cor(to.plot, genes[1], mytype, mywidth = 3, myheight = 3, plot.yaxis = TRUE);  
for(gene in genes[-1])
{
  plot_gene_cor(to.plot, gene, mytype, mywidth = 3, myheight = 3);  
}

####
mytype <- 'HCC';
myexp <- readRDS('data/HCC_2020-03-31_clusterMean.rds');
mypath <- read.table('HCC/QuSAGE_for_heatmap/GeneSet41_Gene_Function-change_heatmap_input.txt', header = TRUE, row.names = 1);
colnames(mypath) <- gsub('luster_', '', colnames(mypath));
genes <- rownames(mycor);
genes <- unique(mygenes)
genes <- genes[order(genes)];
to.plot <- data.frame(cbind(t(mypath), t(myexp[genes, ])));
#plot_gene_cor(to.plot, 'ENO2', mytype);
#plot_gene_cor(to.plot, 'GOLM1', mytype);
plot_gene_cor(to.plot, genes[1], mytype, mywidth = 3, myheight = 3, plot.yaxis = TRUE);  
for(gene in genes[-1])
{
  plot_gene_cor(to.plot, gene, mytype, mywidth = 3, myheight = 3);  
}
####
mytype <- 'PRAD';
myexp <- readRDS('data/PRCA_2020-03-31_cluster_Mean.rds');
mypath <- readRDS('../../tcell/2019-08-01_qusage_tcell_curated_brca.rds');
mypath <- mypath$results.cor[mypath$results.cor=='CD8_T_Cell_Activation_gene', c(2,4)];
rownames(mypath) <- paste0('C', mypath$Cluster-1);
mypath <- data.frame(t(mypath[, -2, drop = FALSE]));
rownames(mypath) <- 'CD8_T_Cell_Activation_gene';

genes <- rownames(mycor);
genes <- unique(mygenes)
genes <- genes[order(genes)];
to.plot <- data.frame(cbind(t(mypath), t(myexp[genes, ])));
#plot_gene_cor(to.plot, 'ENO2', mytype);
#plot_gene_cor(to.plot, 'GOLM1', mytype);
plot_gene_cor(to.plot, genes[1], mytype, mywidth = 3, myheight = 3, plot.yaxis = TRUE);  
for(gene in genes[-1])
{
  plot_gene_cor(to.plot, gene, mytype, mywidth = 3, myheight = 3);  
}

####
