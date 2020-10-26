library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/output');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
dlist <- dir('./', '^Cluster-CNV');
for(i in dlist){
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);
	name <-  gsub('Cluster-CNV-', '', i);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare_all', name, 'rds'));
	assign(paste0('cnv.', name), cnv.all)
};
####
for(i in dlist){
	name <- gsub('Cluster-CNV-', '', i);
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);	
	assign('cnv.all', get(paste0('cnv.', name)));
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type_scnorm);
	cnv.epi <- cnv.all[cnv.all$type=='epithelia', ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	#mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	#mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation', name, 'rds'))
	assign(paste0('mycor.', name), mycor.all);
};

for(i in dlist){
	name <- gsub('Cluster-CNV-', '', i);
	to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot$col <- 'black';
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(iseurat@meta.data)), ]$type_scnorm;
	to.plot$pch <- mypch[to.plot$type]
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score', name, 'pdf')
		)
};
#### use top x % malignant cell CNV and get the highest correlation 
for(i in dlist){
	name <- gsub('Cluster-CNV-', '', i);
	file.value <- paste0('2020-02-14_infercnv_meansquare_all_', name, '.rds');
	cnv.all <- readRDS(file.value)
	#assign('cnv.all', get(paste0('cnv.', name)));
	assign(paste0('cnv.', name), cnv.all);
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type_scnorm);
	cnv.epi <- cnv.all[cnv.all$type=='epithelia', ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- data.frame(cor.5 = c(mycor, mycor.ref));
	####
	cnv.top <- cnv.epi[((nrow(cnv.epi)/20) + 1):(nrow(cnv.epi)/10), ];
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all$cor.10 <- c(mycor, mycor.ref);
	cnv.top <- cnv.epi[((nrow(cnv.epi)/10)+1):(nrow(cnv.epi)/5), ];
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all$cor.20 <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation_p', name, 'rds'))
	assign(paste0('mycor.', name), mycor.all);
};
for(i in dlist){
	mypch <- c(3,16,17)
	names(mypch) <- unique(to.plot$type);
	name <- gsub('Cluster-CNV-', '', i);
	file.value <- paste0('2020-02-14_infercnv_meansquare_all_', name, '.rds');
	cnv.all <- readRDS(file.value);
	to.plot <- cbind(score = get(paste0('cnv.', name)), get(paste0('mycor.', name)));
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	to.plot$pch <- mypch[to.plot$type]
	to.plot$cor <- to.plot$cor.5;
	to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_pmax', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
			points = list(
				pch = mypch,
				col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		)
};
####
#### use cells with the least CNV aberration (ref cells) for correlation calculating 
for(i in dlist){
	name <- gsub('Cluster-CNV-', '', i);
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	avg.top <- rowMeans(mycnv.ref);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation_ref', name, 'rds'))
	assign(paste0('mycor.', name), mycor.all);
};
####
for(i in dlist){
	name <- gsub('Cluster-CNV-', '', i);
	file.value <- paste0('2020-02-14_infercnv_meansquare_all_', name, '.rds');
	cnv.all <- readRDS(file.value)
	assign(paste0('cnv.', name), cnv.all);
	to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot$col <- 'black';
	to.plot[to.plot$score>0.04&to.plot$cor<0, ]$col <- 'red';
	to.plot[to.plot$score<0.04&to.plot$cor>0, ]$col <- 'blue';
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	mypch <- c(3,16,24);
	names(mypch) <- unique(to.plot$type);
	to.plot$pch <- mypch[to.plot$type]
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_ref', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
		points = list(
			pch = mypch,
			col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		)
};
####
library(pheatmap);
library(infercnv);
library(gtools);
load(paste0(i, '/run.final.infercnv_obj'));
mygene <- read.table('../gencode_hg38_gene_pos_replaced_sorted_noHLA.txt', as.is = TRUE);
mygene <- mygene[mygene$V1%in%rownames(mycnv), ];
mygene <- mygene[mixedorder(mygene$V2), ];
mychr <- table(mygene$V2);
mychr <- mychr[mixedorder(names(mychr))];
gaps_col <- sapply(seq(length(mychr)), function(x) sum(mychr[1:x]));
mybreaks <- unique(c(seq(0, 0.5, 0.1), seq(0.5, 1.5, 0.01), seq(1.5, 2, 0.1)));
mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
myorder <- read.table(paste0(i, '/infercnv.observation_groupings.txt'), header = TRUE);
rownames(myorder) <- gsub('-', '.', rownames(myorder));
to.plot <- mycnv[match(mygene$V1, rownames(mycnv)), match(rownames(myorder), colnames(mycnv))];
png(generate.filename('plotcnv', i, 'png'), width = 10, height = 8.2, units = 'in', res = 300);
pheatmap(t(to.plot), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
	breaks = mybreaks, color = mycol, gaps_col = gaps_col);
dev.off();
mygrid <- rep(NA, nrow(to.plot));
mygrid[gaps_col] <- 'black';
p <- create.heatmap(
	x = to.plot,
	filename = NULL,
	colour.scheme = c('blue', 'white', 'red'),
	total.colours = length(mybreaks),
	at = mybreaks,
	cluster.dimensions = 'none',
	row.colour = 'black',
	col.lines = gaps_col,
	col.lwd = 2,
	grid.col = TRUE,
	force.grid.col = TRUE,
	);
png(generate.filename('plotcnv', i, 'png'), width = 10, height = 8.2, units = 'in', res = 300);
print(p);
dev.off();
###
i <- 'JD1800159SL';
	name <- i
	mycnv <- read.table(paste0('2020-02-20_infercnv_ref_cnv0_', i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0('2020-02-20_infercnv_ref_cnv0_', i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);
	assign(paste0('cnv.', name), cnv.all);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare_all', name, 'rds'));
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type_scnorm);
	cnv.epi <- cnv.all[cnv.all$type=='epithelia', ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation', name, 'rds'));
	assign(paste0('mycor.', name), mycor.all);

	#file.value <- paste0('2020-02-14_infercnv_meansquare_all_', name, '.rds');
	#cnv.all <- readRDS(file.value);
	to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	mypch <- c(3,16,17)
	names(mypch) <- unique(to.plot$type);
	to.plot$pch <- mypch[to.plot$type]
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_cnv0', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
			points = list(
				pch = mypch,
				col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		);



i <- 'JD1800159SL';
setwd('../')
for(i in unique(seurat.all@meta.data$orig.ident)){
	name <- i
	mycnv <- read.table(paste0('2020-02-19_infercnv_ref_stroma_unknown_', i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0('2020-02-19_infercnv_ref_stroma_unknown_', i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);
	assign(paste0('cnv.', name), cnv.all);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare_all_stroma_unknown', name, 'rds'));
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type_scnorm);
	cnv.epi <- cnv.all[cnv.all$type=='epithelia', ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation_stroma_unknown', name, 'rds'));
	assign(paste0('mycor.', name), mycor.all);

	name <- i;
	#file.value <- paste0('2020-02-14_infercnv_meansquare_all_', name, '.rds');
	#cnv.all <- readRDS(file.value);
	to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	mypch <- c(3,17)
	names(mypch) <- c("epithelia", "non");
	to.plot$pch <- mypch[to.plot$type]
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_stroma_unknown', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
			points = list(
				pch = mypch,
				col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		);
};
#####
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
dlist <- dir('./', '2020-03-05_infercnv_raw_ref_spike_JD18001');
for(i in dlist){
	name <- gsub('2020-03-05_infercnv_raw_ref_spike_', '', i);
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare_raw_ref', name, 'rds'));
	assign(paste0('cnv.', name), cnv.all);
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type);
	cnv.epi <- cnv.all[cnv.all$type%in%c('Basal/intermediate', 'Luminal'), ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	#mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	#mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation_raw_ref', name, 'rds'))
	assign(paste0('mycor.', name), mycor.all);
};
for(i in dlist){
	name <- gsub('2020-03-02_infercnv_raw_ref_', '', i);
	cnv.all <- readRDS(paste0('2020-03-02_infercnv_meansquare_raw_ref_', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-03-02_infercnv_correlation_raw_ref_', name, '.rds'));
	#to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot <- data.frame(score = cnv.all, cor = mycor.all);
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	mypch <- c(3,16,17)
	names(mypch) <- c("epithelia", "Unknown", "stroma");
	to.plot$pch <- mypch[to.plot$type]
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	if(nrow(to.plot[to.plot$score>0.04&to.plot$cor>0.4, ])>0){
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';		
	}
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	name <- gsub('2020-03-02_infercnv_raw_ref_', '', i);
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_raw_ref', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
			points = list(
				pch = mypch,
				col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		)
	};
###
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
dlist <- dir('./', 'infercnv_raw_ref_scale');
for(i in dlist){
	name <- gsub('.*infercnv_raw_ref_scale', '', i);
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare_raw_ref', name, 'rds'));
	assign(paste0('cnv.', name), cnv.all);
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type);
	cnv.epi <- cnv.all[cnv.all$type%in%c('Basal/intermediate', 'Luminal'), ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	#mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	#mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation_raw_ref', name, 'rds'))
	assign(paste0('mycor.', name), mycor.all);
};
for(i in dlist){
	name <- gsub('.*infercnv_raw_ref_scale', '', i);
	cnv.all <- readRDS(paste0('2020-03-03_infercnv_meansquare_raw_ref_', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-03-03_infercnv_correlation_raw_ref_', name, '.rds'));
	#to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot <- data.frame(score = cnv.all, cor = mycor.all);
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	mypch <- c(3,16,17)
	names(mypch) <- c("epithelia", "Unknown", "stroma");
	to.plot$pch <- mypch[to.plot$type]
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	if(nrow(to.plot[to.plot$score>0.04&to.plot$cor>0.4, ])>0){
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';		
	}
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_raw_ref_scale', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
			points = list(
				pch = mypch,
				col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		)
	};
###
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/v0.8.2');
dlist <- dir('./', 'infercnv_norm_');
for(i in dlist){
	name <- gsub('.*infercnv_norm_stroma_old_', '', i);
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.value <- colMeans((mycnv-1)^2);
	cnv.ref <- colMeans((mycnv.ref-1)^2);
	cnv.all <- c(cnv.value, cnv.ref);
	saveRDS(cnv.all, file = generate.filename('infercnv_meansquare_raw_ref', name, 'rds'));
	assign(paste0('cnv.', name), cnv.all);
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type);
	cnv.epi <- cnv.all[cnv.all$type%in%c('Basal/intermediate', 'Luminal'), ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	#mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	#mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
	mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
	mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
	mycor.all <- c(mycor, mycor.ref);
	saveRDS(mycor.all, file = generate.filename('infercnv_correlation_raw_ref', name, 'rds'))
	assign(paste0('mycor.', name), mycor.all);
};
for(i in dlist){
	name <- gsub('.*infercnv_norm_stroma_old_', '', i);
	cnv.all <- readRDS(paste0('2020-03-03_infercnv_meansquare_raw_ref_', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-03-03_infercnv_correlation_raw_ref_', name, '.rds'));
	#to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
	to.plot <- data.frame(score = cnv.all, cor = mycor.all);
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	mypch <- c(3,16,17)
	names(mypch) <- c("epithelia", "Unknown", "stroma");
	to.plot$pch <- mypch[to.plot$type]
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	if(nrow(to.plot[to.plot$score>0.04&to.plot$cor>0.4, ])>0){
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';		
	}
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.04,
		style = 'Nature',
		filename = generate.filename('cnv_cor_score_norm_stroma_old', name, 'pdf'),
		key = list(
			text = list(
				lab = names(mypch),
				col = 'grey'
				),
			points = list(
				pch = mypch,
				col = 'grey'
			),
		x = 0.8,
		y = 0.1
			)
		)
	};
