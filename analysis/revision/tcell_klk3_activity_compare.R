library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape);
library(plyr);
library(Hmsic);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
seurat.t@meta.data$type <- 'CD8_effector';
seurat.t@meta.data[seurat.t@meta.data$res.0.8%in%c(1,4), ]$type <- 'CD8_naive';
seurat.t@meta.data[seurat.t@meta.data$res.0.8%in%c(0), ]$type <- 'CD4_conv'; 
seurat.t@meta.data[seurat.t@meta.data$res.0.8%in%c(6), ]$type <- 'CD4_Treg';
ref <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GeneSet41_Gene_Function.txt');

to.plot <- data.frame(cd8 = colMeans(seurat.t@data[rownames(seurat.t@data)%in%ref[ref$V2=='CD8_T_Cell_Activation_gene', ]$V1, ]));
to.plot$cytolytic <- colMeans(seurat.t@data[rownames(seurat.t@data)%in%ref[ref$V2=='Cytolytics_effector_pathway_gene', ]$V1, ]);
to.plot$klk3 <- seurat.t@data['KLK3', ];
to.plot$cluster <- seurat.t@meta.data$res.0.8;
to.plot$type <- seurat.t@meta.data$type;
to.plot$group <- ifelse(to.plot$klk3>0, 'yes', 'no');
to.plot$group1 <- as.numeric(as.factor(cut2(to.plot$klk3, g = 7)));

pval <- scientific.notation(wilcox.test(to.plot$cd8~to.plot$group)$p.value);
create.boxplot(
    formula = cd8~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('CD8 T cell activation score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3exp_cd8t', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 1.5,
    text.label = pval
    );
pval <- scientific.notation(wilcox.test(to.plot$cytolytic~to.plot$group)$p.value);
create.boxplot(
    formula = cytolytic~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('Cytolytic effector score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3exp_cytolytic', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 3,
    text.label = pval
    );
###
pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group1%in%c(1,2), ]$cytolytic~to.plot[to.plot$group1%in%c(1,2), ]$group1)$p.value);
pval2 <- scientific.notation(wilcox.test(to.plot[to.plot$group1%in%c(2,3), ]$cytolytic~to.plot[to.plot$group1%in%c(2,3), ]$group1)$p.value);
create.boxplot(
    formula = cytolytic~factor(group1),
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('Cytolytic effector score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3exp_cytolytic_g3', 'pdf'),
    #add.text = TRUE,
    text.x = c(1.5, 2.5),
    text.y = 3,
    text.label = c(pval1, pval2)
    );
###

###
myexp <- seurat.t@data[rownames(seurat.t@data)%in%ref[ref$V2%in%c('CD8_T_Cell_Activation_gene', 'Cytolytics_effector_pathway_gene'), ]$V1, ]
myexp <- data.frame(t(scale(t(myexp))));
to.plot <- data.frame(cd8 = colMeans(myexp[rownames(myexp)%in%ref[ref$V2=='CD8_T_Cell_Activation_gene', ]$V1, ]));
to.plot$cytolytic <- colMeans(myexp[rownames(myexp)%in%ref[ref$V2=='Cytolytics_effector_pathway_gene', ]$V1, ]);
to.plot$klk3 <- seurat.t@data['KLK3', ];
to.plot$cluster <- seurat.t@meta.data$res.0.8;
to.plot$type <- seurat.t@meta.data$type;
to.plot$group <- ifelse(to.plot$klk3>0, 'yes', 'no');
pval <- scientific.notation(wilcox.test(to.plot$cd8~to.plot$group)$p.value);
create.boxplot(
    formula = cd8~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('CD8 T cell activation score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3exp_cd8t_scale', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 1,
    text.label = pval
    );
pval <- scientific.notation(wilcox.test(to.plot$cytolytic~to.plot$group)$p.value);
create.boxplot(
    formula = cytolytic~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('Cytolytic effector score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3exp_cytolytic_scale', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 2,
    text.label = pval
    );
###
to.plot <- to.plot[to.plot$type%in%c('CD8_effector', 'CD8_naive'), ];
pval <- scientific.notation(wilcox.test(to.plot$cd8~to.plot$group)$p.value);
create.boxplot(
    formula = cd8~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('CD8 T cell activation score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell_cd8', 'klk3exp_cd8t_scale', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 0.8,
    text.label = pval
    );
pval <- scientific.notation(wilcox.test(to.plot$cytolytic~to.plot$group)$p.value);
create.boxplot(
    formula = cytolytic~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('Cytolytic effector score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell_cd8', 'klk3exp_cytolytic_scale', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 2,
    text.label = pval
    );
###
to.plot <- to.plot[to.plot$type=='CD8_effector', ];
pval <- scientific.notation(wilcox.test(to.plot$cd8~to.plot$group)$p.value);
create.boxplot(
    formula = cd8~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('CD8 T cell activation score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell_eff', 'klk3exp_cd8t_scale', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 0.8,
    text.label = pval
    );
pval <- scientific.notation(wilcox.test(to.plot$cytolytic~to.plot$group)$p.value);
create.boxplot(
    formula = cytolytic~group,
    data = to.plot,
    #add.stripplot = TRUE,
    xlab.label = 'KLK3 expressed',
    ylab.label = expression('Cytolytic effector score'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell_eff', 'klk3exp_cytolytic_scale', 'pdf'),
    add.text = TRUE,
    text.x = 1.5,
    text.y = 2,
    text.label = pval
    );
