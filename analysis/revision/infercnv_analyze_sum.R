library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/output');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
dlist <- dir('./', '^Cluster-CNV');
dlist <- gsub('.*-', '', dlist);
for(i in dlist){
	name <- i;
	cnv.all <- readRDS(paste0('2020-02-14_infercnv_meansquare_all_', i, '.rds'));
	mycor.all <- readRDS(paste0('2020-02-14_infercnv_correlation_', i, '.rds'));
	to.plot <- data.frame(score = cnv.all, cor = mycor.all);
	to.plot$col <- 'black';
	to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';
	to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
	to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
	to.plot$type_all <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type;
	mypch <- c(3,16,17)
	names(mypch) <- unique(to.plot$type);
	to.plot$pch <- mypch[to.plot$type];
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
		xlimits = c(0, 0.1),
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
		);
};