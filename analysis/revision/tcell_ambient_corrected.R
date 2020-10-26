library(BoutrosLab.plotting.general);
library(SoupX);
library(ggplot2);
library(Matrix);
library(Seurat);
library(reshape);
library(reshape2);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/ambient');
name <- 'corrected';
dataDirs = c("SoupX-HBB-HBA/");
sc = Read10X(dataDirs);
seurat.c <- CreateSeuratObject(raw.data = sc, project = 'corrected');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
myexp <- data.frame(exp = seurat.all@data['KLK3', ]);
myexp$type <- seurat.all@meta.data$type;
myexp$corrected <- log(seurat.c@data['KLK3', rownames(myexp)]+1);
myexp$ori <- log(seurat.all@raw.data['KLK3', ]+1);
to.plot <- myexp[myexp$type=='T', c('corrected', 'ori')];
to.plot <- melt(to.plot);
create.violinplot(
    formula = value~variable,
    data = to.plot,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance'),
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_corrected', 'pdf'),
    #add.text = TRUE,
    #text.x = 1.5,
    #text.y = 7,
    #text.label = pval
    );
####
dataDirs3 = c("SoupX-HBB-HBA-top3/");
sc3 = Read10X(dataDirs);
seurat.c3 <- CreateSeuratObject(raw.data = sc3, project = 'corrected');
myexp$corrected3 <- log(seurat.c3@data['KLK3', rownames(myexp)]+1);
to.plot <- myexp[myexp$type=='T', c('corrected3', 'ori')];
to.plot <- melt(to.plot);
create.violinplot(
    formula = value~variable,
    data = to.plot,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance'),
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_corrected3', 'pdf'),
    #add.text = TRUE,
    #text.x = 1.5,
    #text.y = 7,
    #text.label = pval
    );
####
out <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/2020-03-09_soupx_corrected_ig.rds');
myexp$corrected_ig <- log(out['KLK3', rownames(myexp)] + 1);
to.plot <- myexp[myexp$type=='T', c('corrected_ig', 'ori')];
to.plot <- melt(to.plot);
create.violinplot(
    formula = value~variable,
    data = to.plot,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance'),
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_correctedIG', 'pdf'),
    #add.text = TRUE,
    #text.x = 1.5,
    #text.y = 7,
    #text.label = pval
    );
####
to.plot <- myexp[myexp$type=='T', c('corrected_ig', 'type'), drop = FALSE];
create.violinplot(
    formula = corrected_ig~type,
    data = to.plot,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance (log'[2]*'UMI)'),
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_correctedIG1', 'pdf')
	);

to.plot <- myexp[myexp$type=='T', c('corrected3', 'type'), drop = FALSE];
create.violinplot(
    formula = corrected3~type,
    data = to.plot,
    xlab.label = 'Group',
    ylab.label = expression('KLK3 abundance (log'[2]*'UMI)'),
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_corrected31', 'pdf')
	);
###
#dat <- read.table('~/plot/RDS.NormData.txt', header = TRUE);
#myexp$corrected_b <- dat[match(rownames(myexp), dat$GeneID), ]$KLK3;
#myexp$corrected_3r <- log(seurat.c3@raw.data[sKLK3', rownames(myexp)]+1);
seurat.c1 <- readRDS('GraphClust.seuset-1.rds');
myexp$corrected_b <- log(seurat.c1@raw.data['KLK3', rownames(myexp)]+1);
to.plot <- myexp[myexp$type=='T', c('corrected_b', 'type'), drop = FALSE];
create.violinplot(
    formula = corrected_b~type,
    data = to.plot,
    xlab.label = 'Group',
    #ylab.label = expression('KLK3 abundance (log'[2]*'UMI)'),
    ylab.label = 'log(UMI+1)',
    #xaxis.labels = c('Original', 'Corrected'),
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_correctedHBB_bob', 'pdf')
	);
myexp$raw <- log(seurat.all@raw.data['KLK3', rownames(myexp)]+1);
saveRDS(myexp, file = generate.filename('tcell_KLK3', 'ambient', 'rds'));
####
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
myexp <- readRDS('2020-06-23_tcell_KLK3_ambient.rds');
to.plot <- myexp[myexp$type=='T', c('raw', 'corrected_b', 'corrected_ig'), drop = FALSE];
to.plot <- to.plot[rownames(seurat.t@meta.data), ];
to.plot$exp_ig <- ifelse(to.plot$corrected_ig>0, 'Y', 'N');
to.plot$exp_raw <- ifelse(to.plot$raw>0, 'Y', 'N');
pval <- scientific.notation(wilcox.test(to.plot$raw, to.plot$corrected_b)$p.value);
to.plot <- melt(to.plot);
create.violinplot(
    formula = value~variable,
    data = to.plot,
    xlab.label = 'Group',
    ylab.label = 'log(UMI+1)',
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    style = 'Nature',
    filename = generate.filename('tcell', 'klk3_correctedHBB_bob_raw', 'pdf')
    );
###

to.plot <- data.frame(matrix(c(1545, 1, 0, 1570), nrow =2, byrow = T));
rownames(to.plot) <- colnames(to.plot) <- c('FALSE', 'TRUE');
pdf(generate.filename('tcell', 'numKLK3', 'pdf'))
pheatmap(to.plot, display_number = T, number_color = 'black', number_foramt = "%d", fontsize_number = 20,
    color = colorRampPalette(c('white', 'blue'))(10)[1:8], 
    legend = FALSE, cluster_rows = FALSE, cluster_cols = FALSE);
dev.off();

