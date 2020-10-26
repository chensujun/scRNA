library(BoutrosLab.plotting.general);
library(readxl);
library(Seurat);
library(plyr);
library(gtools);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/multi');
source('/cluster/home/sujunc/chensj/scRNA/script/myfunctions/plot_functions.R');
source('/cluster/projects/hansengroup/sujunc/scRNA/script/myfunctions/plotcnv_functions.R');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
seurat.epi <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/epitumour/objects/2019-08-27_sample13_epi.rds');
mygene <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt', as.is = TRUE);
annot <- seurat.all@meta.data;
rownames(annot) <- gsub('-', '.', rownames(annot));
annot$samp <- annot$orig.ident;
namestem <- '2020-05-18_multiB3_';
name <- unique(annot$samp)[1];
mycnv.all <- read.table(paste0(namestem, mytype, '_', name, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
for(mytype in c('Fibroblast', 'Myeloid', 'T', 'Endothelial', 'Mast')){
    mycnv.value <- read.table(paste0(namestem, mytype, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
    mycnv.ref <- read.table(paste0(namestem, mytype, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
    mycnv.tmp <- cbind(mycnv.value, mycnv.ref[rownames(mycnv.value), ]);
    mycnv.all <- cbind(mycnv.all, mycnv.tmp[rownames(mycnv.all), ]) 
    print(dim(mycnv.all));
};
mycnvi[is.na(mycnvi)] <- 1;
mycnv <- mycnvi;
for(name in unique(annot$samp)[-c(1)]){
    mytype <- 'Epi';
    mycnv.all <- read.table(paste0(namestem, mytype, '_', name, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
    for(mytype in c('Fibroblast', 'Myeloid', 'T', 'Endothelial', 'Mast')){
        mycnv.value <- read.table(paste0(namestem, mytype, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
        mycnv.ref <- read.table(paste0(namestem, mytype, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
        mycnv.tmp <- cbind(mycnv.value, mycnv.ref[rownames(mycnv.value), ]);
        mycnv.all <- cbind(mycnv.all, mycnv.tmp[rownames(mycnv.all), ]) 
        print(dim(mycnv.all));
    };
    mycnvi <- mycnv.all[, colnames(mycnv.all)%in%rownames(annot[annot$samp==name, ])];
    mycnvi[is.na(mycnvi)] <- 1;
    mynames <- unique(c(rownames(mycnv), rownames(mycnvi)));
    mycnv <- cbind(mycnv[mynames, ], mycnvi[mynames, ]);
    rownames(mycnv) <- mynames;
};
#saveRDS(mycnv, file = generate.filename('cnv', 'all', 'rds'));
colnames(mycnv) <- gsub('\\.', '-', colnames(mycnv));
mycnv[is.na(mycnv)] <- 1
saveRDS(mycnv, file = generate.filename('cnv', 'all_naomit_order', 'rds'));
#mycnv <- readRDS('2020-06-25_cnv_all_naomit_order.rds');
#cnv.chr <- data.frame(matrix(nrow = length(unique(mygene$V2)), ncol = ncol(mycnvi)));
#rownames(cnv.chr) <- unique(mygene$V2);
#colnames(cnv.chr) <- colnames(mycnvi)
#for(chr in unique(mygene$V2)){
#	icnv <- mycnvi[rownames(mycnvi)%in%mygene[mygene$V2==chr, ]$V1, ];
#	cnv.chr[chr, ] <- colMeans(icnv);
#};
#cnv.avg <- data.frame(matrix(nrow = nrow(cnv.chr), ncol = 16));
#rownames(cnv.avg) <- rownames(cnv.chr);
#colnames(cnv.avg) <- paste0('C', unique(seurat.epi@meta.data$fig.cluster));
#for(i in unique(seurat.epi@meta.data$fig.cluster)){
#	cnv.avg[, paste0('C', i)] <- rowMeans(cnv.chr[, rownames(seurat.epi@meta.data[seurat.epi@meta.data$fig.cluster==i, ])])
#};
#dbreaks <- seq(-1, 1, .1);
#col <- colorRampPalette(c('#1c2f5d','white', '#6e121f'))(length(dbreaks));
#pheatmap(cnv.avg-1, cluster_rows = FALSE, cluster_cols = FALSE, breaks = dbreaks, color = col);
#dev.off();
#####
mycnv <- readRDS('2020-06-26_cnv_all_naomit_order.rds')
mycnvi <- mycnv[, rownames(seurat.epi@meta.data)];
cnv.cl <- data.frame(matrix(nrow = nrow(mycnvi), ncol = 16));
rownames(cnv.cl) <- rownames(mycnvi);
colnames(cnv.cl) <- paste0('C', unique(seurat.epi@meta.data$fig.cluster));
for(i in unique(seurat.epi@meta.data$fig.cluster)){
	cnv.cl[, paste0('C', i)] <- rowMeans(mycnvi[, rownames(seurat.epi@meta.data[seurat.epi@meta.data$fig.cluster==i, ])])
};
mygene$pos <- 1;
for(chr in unique(mygene$V2)){
	num <- floor(nrow(mygene[mygene$V2==chr, ])/100)
	mygene[mygene$V2==chr, ]$pos <- c(rep(seq(num), each = 100), rep(num+1, nrow(mygene[mygene$V2==chr, ])-num*100));
};
mygene$group <- paste0(mygene$V2, '_', mygene$pos);
cnv.cl$group <- mygene[match(rownames(cnv.cl), mygene$V1), ]$group;
cnv.cl <- na.omit(cnv.cl);
cnv.avg100 <- ddply(cnv.cl, 'group', numcolwise(mean));
rownames(cnv.avg100) <- cnv.avg100$group;
cnv.avg100 <- cnv.avg100[, -1];
cnv.avg100 <- cnv.avg100[mixedorder(rownames(cnv.avg100)), ];
cnv.avg100 <- cnv.avg100[, paste0('C', seq(0, 15))]
mychr <- table(gsub('_.*', '', rownames(cnv.avg100)));
mychr <- mychr[mixedorder(names(mychr))];
gaps_row <- sapply(seq(length(mychr)), function(x) 260-sum(mychr[1:x]))
save(cnv.cl, cnv.avg100, file = generate.filename('cnv', 'avg_order', 'rda'));

dbreaks <- seq(-1, 1, length.out = 21)
col <- colorRampPalette(c('blue','white', 'red'))(21)[c(1,3,5,9, 10:12, 13, 17, 19, 21)];
pdf(generate.filename('cnv_avg', 'order', 'pdf'));
pheatmap(cnv.avg100-1, cluster_rows = FALSE, cluster_cols = FALSE, breaks = dbreaks, color = col,
	gaps_row = sapply(seq(length(gaps_row)), function(x) sum(gaps_row[1:x])), show_rownames = FALSE);
dev.off();
###
create.heatmap(
		x = cnv.avg100-1,
		filename = generate.filename('cnv_avg', 'order', 'pdf'),
		same.as.matrix = TRUE,
		colour.scheme = c('blue', 'white', 'red'),
		total.colours = length(dbreaks) + 1,
		at = dbreaks,
		cluster.dimensions = 'none',
		row.colour = 'black',
		#col.lines = gaps_col,
		col.lwd = 0.5,
		row.lwd = 0.5,
		row.lines = gaps_row,
		covariate.legend = NULL,
		grid.col = FALSE,
		force.grid.col = TRUE,
		grid.row = TRUE,
		force.grid.row = TRUE,
		print.colour.key = TRUE,
		colourkey.cex = 1,
		xaxis.lab = seq(0, 15),
		yaxis.lab = NULL,
		xaxis.rot = 0,
		xaxis.fontface = 'plain'
		);

###
mygene.rs <- mygene[, 1:4];
mygene.rs$V1 <- mygene$V1[order(mygene$V1)];
write.table(mygene.rs, 'gencode_hg38_gene_pos_replaced_sortedName_noHLA.txt', quote = F, row.names = F, col.names = F, sep = '\t');
### shuffling the genes
namestem <- '2020-06-25_multiB3NameSorted_';
name <- unique(annot$samp)[1];
mytype <- 'Epi';
mycnv.all <- read.table(paste0(namestem, mytype, '_', name, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
for(mytype in c('Fibroblast', 'Myeloid', 'T', 'Endothelial', 'Mast')){
    mycnv.value <- read.table(paste0(namestem, mytype, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
    mycnv.ref <- read.table(paste0(namestem, mytype, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
    mycnv.tmp <- cbind(mycnv.value, mycnv.ref[rownames(mycnv.value), ]);
    mycnv.all <- cbind(mycnv.all, mycnv.tmp[rownames(mycnv.all), ]) 
    print(dim(mycnv.all));
};
mycnvi <- mycnv.all[, colnames(mycnv.all)%in%rownames(annot[annot$samp==name, ])];
mycnvi[is.na(mycnvi)] <- 1;
mycnv <- mycnvi;
#for(name in unique(annot$samp)[-c(1)]){
#for(name in c('JD1800153SL', 'JD1800154SL', 'JD1800155SL', 'JD1800171SL', 'JD1800175SL', 'JD1800176SL', 'JD1800177SL')){
#for(name in c('JD1800162SL', 'JD1800172SL', 'JD1800173SL', 'JD1800175SL')){
for(name in unique(annot$samp)[!unique(annot$samp)=='JD1800156SL']){
   mytype <- 'Epi';
    mycnv.all <- read.table(paste0(namestem, mytype, '_', name, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
    for(mytype in c('Fibroblast', 'Myeloid', 'T', 'Endothelial', 'Mast')){
        mycnv.value <- read.table(paste0(namestem, mytype, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
        mycnv.ref <- read.table(paste0(namestem, mytype, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
        mycnv.tmp <- cbind(mycnv.value, mycnv.ref[rownames(mycnv.value), ]);
        mycnv.all <- cbind(mycnv.all, mycnv.tmp[rownames(mycnv.all), ]) 
        print(dim(mycnv.all));
    };
    mycnvi <- mycnv.all[, colnames(mycnv.all)%in%rownames(annot[annot$samp==name, ])];
    mycnvi[is.na(mycnvi)] <- 1;
    mynames <- unique(c(rownames(mycnv), rownames(mycnvi)));
    mycnv <- cbind(mycnv[mynames, ], mycnvi[mynames, ]);
    rownames(mycnv) <- mynames;
};
#saveRDS(mycnv, file = generate.filename('cnv', 'all', 'rds'));
###
colnames(mycnv) <- gsub('\\.', '-', colnames(mycnv));
mycnv[is.na(mycnv)] <- 1
#saveRDS(mycnv, file = generate.filename('cnv', 'all_naomit_1', 'rds'));
saveRDS(mycnv, file = generate.filename('cnv', 'all_naomit_NameSort', 'rds'));
####
mycnvi <- mycnv[, colnames(mycnv)%in%rownames(seurat.epi@meta.data)];
cnv.cl <- data.frame(matrix(nrow = nrow(mycnvi), ncol = 16));
rownames(cnv.cl) <- rownames(mycnvi);
colnames(cnv.cl) <- paste0('C', unique(seurat.epi@meta.data$fig.cluster));
for(i in unique(seurat.epi@meta.data$fig.cluster)){
	cnv.cl[, paste0('C', i)] <- rowMeans(mycnvi[, colnames(mycnvi)%in%rownames(seurat.epi@meta.data[seurat.epi@meta.data$fig.cluster==i, ])])
};
mygene$pos <- 1;
for(chr in unique(mygene$V2)){
	num <- floor(nrow(mygene[mygene$V2==chr, ])/100)
	mygene[mygene$V2==chr, ]$pos <- c(rep(seq(num), each = 100), rep(num+1, nrow(mygene[mygene$V2==chr, ])-num*100));
};
mygene$group <- paste0(mygene$V2, '_', mygene$pos);
cnv.cl$group <- mygene[match(rownames(cnv.cl), mygene$V1), ]$group;
cnv.cl <- na.omit(cnv.cl);
cnv.avg100 <- ddply(cnv.cl, 'group', numcolwise(mean));
rownames(cnv.avg100) <- cnv.avg100$group;
cnv.avg100 <- cnv.avg100[, -1];
cnv.avg100 <- cnv.avg100[mixedorder(rownames(cnv.avg100)), ];
cnv.avg100 <- cnv.avg100[, paste0('C', seq(0, 15))]
mychr <- table(gsub('_.*', '', rownames(cnv.avg100)));
mychr <- mychr[mixedorder(names(mychr))];
gaps_row <- sapply(seq(length(mychr)), function(x) 260-sum(mychr[1:x]))
save(cnv.cl, cnv.avg100, file = generate.filename('cnv', 'avg_NameSort', 'rda'));

#dbreaks <- unique(c(seq(-1, -0.5, 0.25), seq(-0.5, 0.5, 0.1), seq(0.5, 1, 0.25)));
dbreaks <- seq(-1, 1, length.out = 21)
#col <- colorRampPalette(c('#1c2f5d','white', '#6e121f'))(length(dbreaks));
col <- colorRampPalette(c('blue','white', 'red'))(21)[c(1,3,5,9, 10:12, 13, 17, 19, 21)];
pdf(generate.filename('cnv_avg', 'NameSort', 'pdf'));
pheatmap(cnv.avg100-1, cluster_rows = FALSE, cluster_cols = FALSE, breaks = dbreaks, color = col,
	gaps_row = sapply(seq(length(gaps_row)), function(x) sum(gaps_row[1:x])), show_rownames = FALSE);
dev.off();
###
create.heatmap(
		x = cnv.avg100-1,
		filename = generate.filename('cnv_avg', 'NameSort', 'pdf'),
		same.as.matrix = TRUE,
		colour.scheme = c('blue', 'white', 'red'),
		total.colours = length(dbreaks) + 1,
		at = dbreaks,
		cluster.dimensions = 'none',
		row.colour = 'black',
		#col.lines = gaps_col,
		col.lwd = 0.5,
		row.lwd = 0.5,
		row.lines = gaps_row,
		covariate.legend = NULL,
		grid.col = FALSE,
		force.grid.col = TRUE,
		grid.row = TRUE,
		force.grid.row = TRUE,
		print.colour.key = TRUE,
		colourkey.cex = 1,
		xaxis.lab = seq(0, 15),
		yaxis.lab = NULL,
		xaxis.rot = 0,
		xaxis.fontface = 'plain'
		);

p1 <- create.heatmap(
		x = cnv.avg100-1,
		filename = NULL,
		same.as.matrix = TRUE,
		colour.scheme = c('blue', 'white', 'red'),
		total.colours = length(dbreaks) + 1,
		at = dbreaks,
		cluster.dimensions = 'none',
		row.colour = 'black',
		#col.lines = gaps_col,
		col.lwd = 0.5,
		row.lwd = 0.5,
		row.lines = gaps_row,
		covariate.legend = NULL,
		grid.col = FALSE,
		force.grid.col = TRUE,
		grid.row = TRUE,
		force.grid.row = TRUE,
		print.colour.key = TRUE,
		colourkey.cex = 1,
		xaxis.lab = seq(0, 15),
		yaxis.lab = NULL,
		xaxis.rot = 0
		);
pdf(generate.filename('cnv_avg', 'NameSort', 'pdf'))
print(p1)
dev.off();
