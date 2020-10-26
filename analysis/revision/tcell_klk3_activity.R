library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape);
library(plyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
#load('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/2019-08-01_tcell_curated_brca_qusage_results_10_100.rda');
#myresult <- list(
#	results.cor = results.cor, 
#	results.clust.id = results.clust.id
#	);
#saveRDS(myresult, file = generate.filename('tcell', 'activity', 'rds'));
myresult <- readRDS('2020-03-13_tcell_activity.rds');
mypath <- myresult[[1]];
mypath <- mypath[grep('CD8|Pro|Treg|Cyto', mypath$pathway.name), ];
mypath.fc <- cast(mypath, pathway.name~Cluster, mean, value = 'log.fold.change');
rownames(mypath.fc) <- mypath.fc$pathway.name;
mypath.fc <- mypath.fc[, -1];

dat.klk <- data.frame(exp =seurat.t@data['KLK3', ], cluster = seurat.t@meta.data$res.0.8);
dat.klk <- ddply(dat.klk, 'cluster', numcolwise(median));
dat.klk[, 3:6] <- t(mypath.fc);
colnames(dat.klk)[3:6] <- c('CD8TAct', 'Cytolytics', 'ProInflam', 'Treg')

create.scatterplot(
	formula = CD8_T_Cell_Activation_gene~exp,
	data = dat.klk,
	xlab.label = 'mean KLK3',
	ylab.label = 'CD8_T_Cell_Activation_signature',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	filename = generate.filename('KLK3_CD8activation', 'tcell', 'pdf'),

);

create.scatterplot(
	formula = CD8_T_Cell_Activation_gene~exp,
	data = dat.klk,
	xlab.label = 'mean KLK3',
	ylab.label = 'CD8_T_Cell_Activation_signature',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	xlimits = c(0, 2.5),
	filename = generate.filename('KLK3_CD8activation', 'tcell', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = dat.klk$exp,
                y = dat.klk$CD8_T_Cell_Activation_gene,
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
	formula = Cytolytics_effector_pathway_gene~exp,
	data = dat.klk,
	xlab.label = 'mean KLK3',
	ylab.label = 'Cytolytics_effector_pathway_signature',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	xlimits = c(0, 2.5),
	filename = generate.filename('KLK3_Cytolytics', 'tcell', 'pdf'),
    legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = dat.klk$exp,
                y = dat.klk$Cytolytics_effector_pathway_gene,
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

);
####
dat.klk <- data.frame(exp =seurat.t@data['KLK3', ], cluster = seurat.t@meta.data$res.0.8);
dat.klk$group <- factor(ifelse(dat.klk$cluster%in%c(0,1,4,2 ), 'not', 'active'));
pval <- scientific.notation(wilcox.test(dat.klk$exp~dat.klk$group)$p.value);
create.violinplot(
    formula = exp~group,
    data = dat.klk,
    #add.stripplot = TRUE,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'group2', 'pdf'),
    #add.text = TRUE,
    #text.x = 1.5,
    #text.y = 7,
    #text.label = pval
    );
####
flist <- list.files('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/pathway_from_brca', full.name = TRUE);
gs <- list();
for(i in flist){
    ifile <- read.table(i, header = TRUE);
    gs[[gsub('.*/|\\..*', '', i)]] <- ifile[, 1]
};

dat.klk <- data.frame(exp =seurat.t@data['KLK3', ], cluster = seurat.t@meta.data$res.0.8);
dat.klk$t_act <- colMeans(seurat.t@data[rownames(seurat.t@data)%in%gs[['CD8_T_Cell_Activation_gene']], ])
dat.klk$pro_infla <- colMeans(seurat.t@data[rownames(seurat.t@data)%in%gs[['Pro_inflammatory_gene']], ])
dat.klk$cyto <- colMeans(seurat.t@data[rownames(seurat.t@data)%in%gs[['Cytolytics_effector_pathway_gene']], ])
create.scatterplot(
    formula = t_act~exp,
    data = dat.klk,
    xlab.label = 'KLK3 abundance',
    ylab.label = 'CD8_T_Cell_Activation_signature',
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('KLK3_CD8activation', 'tcell_cell', 'pdf')
);
####
to.plot <- dat.klk[dat.klk$cluster%in%c(2,3,5), ];
to.plot$group <- ifelse(to.plot$cluster==5, 'low', 'high');
pval <- scientific.notation(wilcox.test(to.plot$exp~to.plot$group)$p.value);
create.boxplot(
    formula = exp~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'CD8 T cell activation score',
    ylab.label = expression('KLK3 abundance'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'group_effector', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 7.3,
    text.label = pval
    );
####
plot_gene_cor <- function(to.plot, gene, mypath, myname, myheight = 6, mywidth = 6, plot.yaxis = FALSE){
    to.plot$gene <- to.plot[, gene];
    xlab.max <- max(to.plot$gene);
    legend.x <- -(xlab.max/10);
    xlab.max <- ifelse(xlab.max<0.01, scientific.notation(xlab.max), round(xlab.max, 2));
    to.plot$CD8_T_Cell_Activation_gene <- to.plot[, mypath]
    fit <- lm(to.plot$CD8_T_Cell_Activation_gene~to.plot$gene);
    key = get.corr.key(
                    x = to.plot[, gene],
                    y = to.plot$CD8_T_Cell_Activation_gene,
                    label.items = c('spearman','spearman.p'),
                    alpha.background = 0,
                    key.cex = .5
                    );
    mycol <- ifelse(cor.test(to.plot$gene, to.plot$CD8_T_Cell_Activation_gene, method = 'spearman')$p.value<0.05, 'red', 'black');
    xrange <- round(range(to.plot$CD8_T_Cell_Activation_gene), 2);
    pdf(generate.filename(mypath, paste0(gene, '_', myname), 'pdf'), height = myheight, width = mywidth);
    par(mar = rep(3, 4));
    plot(to.plot$gene, to.plot$CD8_T_Cell_Activation_gene, cex = .5, col = mycol,
        pch = 20, xlab = '', ylab = '', axes = FALSE, ylim = c(xrange[1], xrange[2]), xlim = c(0, max(to.plot$gene)));
    axis(1, at = c(0, max(to.plot$gene)), labels = c(0, xlab.max));
    if(plot.yaxis){
    axis(2, las = 2, at = c(xrange[1], xrange[2]), );
    };
    box(bty = 'l');
    abline(fit, col = mycol);
    legend(legend.x, 0.45, legend = key$text$lab, bty = 'n', cex = 0.5);
    dev.off();

};
to.plot <- data.frame(cd8 = colMeans(seurat.t@data[rownames(seurat.t@data)%in%gs[['CD8_T_Cell_Activation_gene']], ]));
to.plot$treg <- colMeans(seurat.t@data[rownames(seurat.t@data)%in%gs[['Treg_gene']], ]);
to.plot$cyto <- colMeans(seurat.t@data[rownames(seurat.t@data)%in%gs[['Cytolytics_effector_pathway_gene']], ]);
to.plot$klk3 <- seurat.t@data['KLK3', ]
to.plot$cluster <- seurat.t@meta.data$res.0.8;
for(i in c(2, 3, 5))
{
  plot_gene_cor(to.plot[to.plot$cluster==i, ], 'klk3', 'cd8', i, mywidth = 3, myheight = 3, plot.yaxis = TRUE); 
  plot_gene_cor(to.plot[to.plot$cluster==i, ], 'klk3', 'cyto', i, mywidth = 3, myheight = 3, plot.yaxis = TRUE);  
};
i <- 6;
plot_gene_cor(to.plot[to.plot$cluster==i, ], 'klk3', 'treg', i, mywidth = 3, myheight = 3, plot.yaxis = TRUE);  