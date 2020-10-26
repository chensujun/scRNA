library(BoutrosLab.plotting.general);
library(Seurat);
library(matrixStats);
library(pheatmap);
library(gtools);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/output/');
source('/cluster/projects/hansengroup/sujunc/scRNA/script/myfunctions/plot_functions.R');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
mygene <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt', as.is = TRUE);
lookup <- read.table('lookup.txt', as.is = TRUE);
dlist <- dir('./', '^Cluster-CNV');
for(i in dlist){
	#i <- dlist[1];
	print(i);
	name <- lookup[lookup$V1==gsub('Cluster-CNV-Sample', '', i), ]$V2;
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mygene <- mygene[mygene$V1%in%rownames(mycnv), ];
	#mygene <- mygene[mixedorder(mygene$V2), ];
	mychr <- table(mygene$V2);
	mychr <- mychr[mixedorder(names(mychr))];
	gaps_col <- sapply(seq(length(mychr)), function(x) sum(mychr[1:x]));
	mybreaks <- unique(c(seq(0, 0.5, 0.1), seq(0.5, 1.5, 0.01), seq(1.5, 2, 0.1)));
	#mybreaks <- seq(0, 2, length.out = 16);
	#mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
	custom_pal <- color.palette(c("darkblue", "white", "darkred"),c(2, 2));
	mycol <- custom_pal(length(mybreaks) + 1);
	myorder <- read.table(paste0(i, '/infercnv.observation_groupings.txt'), header = TRUE);
	rownames(myorder) <- gsub('-', '.', rownames(myorder));
	to.plot <- mycnv[match(mygene$V1, rownames(mycnv)), match(rownames(myorder), colnames(mycnv))];
	mygrid <- rep(NA, nrow(to.plot));
	mygrid[gaps_col] <- 'black';
	### create top covariate for chromosome annotation
	col.chr <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(unique(names(mychr))));
	names(col.chr) <- unique(names(mychr));
	col.type <- c('#66C2A5', '#FFED6F');
	names(col.type) <- c('non', 'epithelia');
	if(name %in%c('JD1800153SL','JD1800172SL')){
		chr.covariate <- list(
			rect = list(
				col = rep(col.chr, times = mychr),
				fill = rep(col.chr, times = mychr)
				)
			)
		color.key <- FALSE
	}else if(name == 'JD1800177SL'){
		chr.covariate <- NULL
		color.key <- FALSE

	}else {
		chr.covariate <- NULL
		color.key <- FALSE
		cov.legend <- NULL
	};
	#mytype <- seurat.all@meta.data[gsub('\\.', '-', colnames(to.plot)), ]$type_scnorm;
	mytype <- ifelse(seurat.all@meta.data[gsub('\\.', '-', colnames(to.plot)), ]$type%in%c('Basal/intermediate', 'Luminal'), 'epithelia', 'non');
	to.plot <- cbind(to.plot[, mytype=='epithelia'], to.plot[, mytype=='non']);
	mytype <- ifelse(seurat.all@meta.data[gsub('\\.', '-', colnames(to.plot)), ]$type%in%c('Basal/intermediate', 'Luminal'), 'epithelia', 'non');
	gaps_row <- length(mytype[mytype==mytype[2]]);
	mygrid.row <- rep(NA, ncol(to.plot));
	mygrid.row[gaps_row] <- 'black';

	cell.covariate <- list(
		rect = list(
			col = col.type[mytype],
			fill = col.type[mytype]
			)
		);
	print('start ploting')
	p <- create.heatmap(
		x = to.plot,
		filename = NULL,
		colour.scheme = mycol,
		total.colours = length(mybreaks) + 1,
		at = mybreaks,
		cluster.dimensions = 'none',
		row.colour = 'black',
		col.lines = gaps_col,
		col.lwd = 0.5,
		row.lwd = 0.5,
		row.lines = gaps_row,
		covariates.top = chr.covariate,
		covariates = cell.covariate,
		covariate.legend = NULL,
		grid.col = TRUE,
		force.grid.col = TRUE,
		grid.row = TRUE,
		force.grid.row = TRUE,
		print.colour.key = color.key
		);
	png(generate.filename('plotcnv', name, 'png'), width = 10, height = 3.5, units = 'in', res = 300);
	print(p);
	dev.off();
	assign(paste0('p.', name), p);
	save(p, file = generate.filename('plotcnv', name, 'rda'));
};

cov.legend <- list(
	legend = list(
		colours = col.type,
		border = 'white',
		labels = c('unknown', 'epithelia'),
		title = 'Group'
		)
	);
color.key = TRUE;
p <- create.heatmap(
	x = to.plot[1:10, 1:10],
	filename = NULL,
	colour.scheme = mycol,
	total.colours = length(mybreaks) + 1,
	at = mybreaks,
	cluster.dimensions = 'none',
	row.colour = 'black',
	col.lines = gaps_col,
	col.lwd = 0.5,
	row.lwd = 0.5,
	row.lines = gaps_row,
	covariates.top = chr.covariate,
	covariates = cell.covariate,
	covariate.legend = cov.legend,
	grid.col = TRUE,
	force.grid.col = TRUE,
	grid.row = TRUE,
	force.grid.row = TRUE,
	print.colour.key = color.key,
	colourkey.cex = 1.5
	);
pdf(generate.filename('plotcnv', 'legend', 'pdf'), width = 10, height = 3.5);
print(p);
dev.off();

for(i in unique(seurat.all@meta.data$orig.ident)){
	name <- i;
	print(i)
	p <- get(paste0('p.', i));
	png(generate.filename('plotcnv', name, 'png'), width = 10, height = 3.5, units = 'in', res = 100);
	print(p);
	dev.off();
}

cov.legend <- list(
	legend = list(
		colours = col.type,
		labels = c('unknown', 'epithelia'),
		title = 'Group'
		)
	);
legend <- legend.grob(
	cov.legend
	)
pall <- create.multiplot(
    plot.objects = list(
		p.JD1800153SL,
		p.JD1800154SL,
		p.JD1800155SL,
		p.JD1800159SL
        ),
    plot.layout = c(1, 4),
    panel.heights = rep(1, 4), 
    panel.widths = 1,
    y.spacing = rep(-1, 2),
    x.spacing = c(1, 1),    
    x.relation = 'free',
    y.relation = 'free',
    yaxis.tck = 0, 
    xaxis.tck = 0,
    print.new.legend = TRUE,
    right.padding = 20,
    left.padding = 0.01,
    yaxis.cex = 1,
    legend = list(
        inside = list(
            x = 1.05,
            y = 1,
            fun = legend
            )
        ),
    style = 'Nature',
    height = 5,
    width = 8
    );

png(generate.filename('plotcnv', name, 'png'), width = 10, height = 8, units = 'in', res = 100);
pall;
dev.off();

####
dlist <- dir('./', '^Cluster-CNV');
for(i in dlist){
	name <- gsub('Cluster-CNV-Sample', '', i);
	cnv.all <- readRDS(paste0('2020-02-14_infercnv_meansquare_all_Sample', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-02-14_infercnv_correlation_Sample', name, '.rds'));
	name <- lookup[lookup$V1==name, ]$V2;
	to.plot <- data.frame(score = cnv.all, cor = mycor.all);
	rownames(to.plot) <- gsub('\\.', '-', rownames(to.plot));
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	to.plot$klk3 <- seurat.all@data['KLK3', rownames(to.plot)];
	to.plot$pmsa <- seurat.all@data['FOLH1', rownames(to.plot)];
	to.plot$erg <- seurat.all@data['ERG', rownames(to.plot)];
	dat.exp <- seurat.all@data[, colnames(seurat.all@data)%in%rownames(to.plot[to.plot$type=='epithelia', ])];
	rv <- rowVars(dat.exp);
	dat.exp <- dat.exp[order(rv, decreasing=TRUE)[seq_len(1000)], ];
	hc <- hclust(dist(t(dat.exp)));
	dat.exp <- data.frame(t(scale(t(dat.exp))));
	mybreaks <- unique(c(seq(-4, -1, 1), seq(-1, 1, .1), seq(1, 4, 1)));
	mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
	ann_col <- to.plot[, 1:2];
	rownames(ann_col) <- gsub('-', '.', rownames(ann_col));
	ann_col$score <- ifelse(ann_col$score>0.05, 'Y', 'N');
	ann_col$cor <- ifelse(ann_col$cor>0.4, 'Y', 'N');
	ann_col <- ann_col[colnames(dat.exp[, hc$order]), ]
	png(generate.filename('plotexp_cnv', i, 'png'), width = 10, height = 8.2, units = 'in', res = 300);
	pheatmap((dat.exp[, hc$order]), cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
		breaks = mybreaks, color = mycol, annotation_col = ann_col);
	dev.off();

};
###
for(i in dlist){
	name <- gsub('Cluster-CNV-Sample', '', i);
	cnv.all <- readRDS(paste0('2020-02-14_infercnv_meansquare_all_Sample', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-02-14_infercnv_correlation_Sample', name, '.rds'));
	name <- lookup[lookup$V1==name, ]$V2;
	to.plot <- data.frame(score = cnv.all, cor = mycor.all);
	rownames(to.plot) <- gsub('\\.', '-', rownames(to.plot));
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type;
	to.plot <- ddply(to.plot, 'type', numcolwise(mean));
	assign(paste0('to.plot.', name), to.plot);
	create.scatterplot(
		formula = cor~score,
		data = to.plot,
		xlab.label = '',
		ylab.label = '',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.05,
		style = 'Nature',
		filename = generate.filename('cnv_cor_exp_type', name, 'pdf'),
		)


}
###
for(i in dlist){
	name <- gsub('Cluster-CNV-Sample', '', i);
	mycnv <- read.table(paste0(i, '/infercnv.observations.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	mycnv.ref <- read.table(paste0(i, '/infercnv.references.txt'), row.names = 1, as.is = TRUE, header=TRUE);
	cnv.all <- readRDS(paste0('2020-02-14_infercnv_meansquare_all_Sample', name, '.rds'));	
	mycor.all <- readRDS(paste0('2020-02-14_infercnv_correlation_Sample', name, '.rds'));
	name <- lookup[lookup$V1==name, ]$V2;	
	assign(paste0('cnv.', name), cnv.all);
	assign(paste0('mycor.', name), mycor.all);
	cnv.all <- data.frame(cnv = cnv.all, type = seurat.all@meta.data[gsub('\\.', '-', names(cnv.all)), ]$type);
	cnv.epi <- cnv.all[cnv.all$type%in%c('Basal/intermediate', 'Luminal'), ];
	cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
	cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
	iexp <- seurat.all@data[, gsub('\\.', '-', c(colnames(mycnv), colnames(mycnv.ref)))];
	colnames(iexp) <- gsub('-', '.', colnames(iexp));
	avg.top <- rowMeans(iexp[, rownames(cnv.top)]);
	mycor.exp <- apply(iexp, 2, function(x) cor(x, avg.top));
	saveRDS(mycor.exp, file = generate.filename('infercnv_correlation_exp', name, 'rds'))
	assign(paste0('mycor.exp.', name), mycor.all);
};


for(i in dlist){
	name <- gsub('Cluster-CNV-Sample', '', i);
	cnv.all <- readRDS(paste0('2020-02-14_infercnv_meansquare_all_Sample', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-02-14_infercnv_correlation_Sample', name, '.rds'));
	name <- lookup[lookup$V1==name, ]$V2;
	mycor.exp <- readRDS(paste0('2020-03-03_infercnv_correlation_exp_', name, '.rds'));
	to.plot <- data.frame(score = cnv.all, cor = mycor.all, exp.cor = mycor.exp);
	rownames(to.plot) <- gsub('\\.', '-', rownames(to.plot));
	#to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	to.plot$type_all <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type;
	to.plot$type <- ifelse(to.plot$type_all%in%c('Basal/intermediate', 'Luminal'),  'epithelia', 'non');
	#mypch <- c(3,16,17)
	#names(mypch) <- c("epithelia", "Unknown", "stroma");
	mypch <- c(3,16);
	names(mypch) <- c("epithelia", "non");
	to.plot$pch <- mypch[to.plot$type]
	#to.plot$cor <- to.plot$cor.5;
	#to.plot[to.plot$type=='epithelia', ]$cor <- apply(to.plot[to.plot$type=='epithelia', 2:4], 1, max);
	to.plot$col <- 'black';
	if(nrow(to.plot[to.plot$exp.cor>0.7&to.plot$cor>0.4, ])>0){
	to.plot[to.plot$exp.cor>0.7&to.plot$cor>0.4, ]$col <- 'red';		
	}
	to.plot[to.plot$exp.cor<0.7&to.plot$cor<0.4, ]$col <- 'blue';
	create.scatterplot(
		formula = cor~exp.cor,
		data = to.plot,
		col = to.plot$col,
		pch = to.plot$pch,
		xlab.label = 'Expression correlation',
		ylab.label = 'CNV correlation',
		xaxis.fontface = 'plain', 
		yaxis.fontface = 'plain',
		abline.h = 0.4,
		abline.v = 0.7,
		style = 'Nature',
		filename = generate.filename('cnv_cor_exp', name, 'pdf'),
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
