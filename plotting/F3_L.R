library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape);
library(plyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/tcell');
#seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
#name <- 'all';
#seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monocytic', seurat.all@meta.data$type);
#seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
#seurat.all <- SetAllIdent(seurat.all, id = 'type');
#myexp <- data.frame(exp = seurat.all@data['KLK3', ]);
#myexp$type <- seurat.all@meta.data$type;
#myexp <- myexp[myexp$exp>0, ]
#myexp$group <- 'epithelia';
#myexp[myexp$type%in%c('Endothelia', 'Fibroblast'), ]$group <- 'stroma';
#myexp[!myexp$type%in%c('Endothelia', 'Fibroblast', 'Basal/intermediate', 'Luminal'), ]$group <- 'immune';
#myexp$group <- factor(myexp$group);#

####
#to.plot <- data.frame(exp = seurat.all@data['KLK3', ], type = seurat.all@meta.data$type, 
#	samp = 'sample13');
#to.plot <- to.plot[to.plot$type!='Epithelia', ]
#to.plot.m <- data.frame(reshape::cast(to.plot, type~samp, mean, value = 'exp', fill = 0));
#to.plot.f <- data.frame(reshape::cast(to.plot[to.plot$exp>0, ], type~samp, length, value = 'exp', fill = 0));
#to.plot.n <- data.frame(reshape::cast(to.plot, type~samp, length, value = 'exp', fill = 0));
##
#to.plot.n$group <- 'all';
#to.plot.f$group <- 'pos';
#to.plot <- rbind(to.plot.f, to.plot.n);
#to.plot.f[, 2] <- round(100*to.plot.f[, 2]/to.plot.n[, 2]);
#to.plot.n$group <- 'all';
#to.plot.f$group <- 'pos';
#to.plot <- rbind(to.plot.f, to.plot.n);
#to.plot.f[, 2:4] <- sapply(seq(3), function(x) 100*to.plot.f[, (x+1)]/to.plot.n[, (x+1)]);
#to.plot.f[, 2:4] <- round(to.plot.f[, 2:4])
#to.plot.f[6, 3] <- 'NA';#

#create.barplot(
#	data = to.plot[, c('type', 'S5.LEFT.2', 'group')],
#	formula = S5.LEFT.2~type,
#	groups = to.plot$group,
#	col = default.colours(2),
#	style = 'Nature',
#	#xaxis.lab = c('R', 'L', 'T'),
#	xlab.label = '',
#	ylab.labe = '# cells',
#	xaxis.cex = 2,
#	yaxis.cex = 2,
#	xaxis.rot = 45,
#	ylimits = c(0, 9000),
#	height = 4,
#	width = 8,
#	add.text = TRUE,
#	text.x = seq(10),
#	text.y = to.plot[11:20, ]$S5.LEFT.2 + 500,
#	text.labels = paste0(to.plot.f$S5.LEFT.2, '%'),
#	text.cex = 1.5,
#	text.fontface = 'plain',
#	key = list(
#	text = list(
#		lab = c('All', 'KLK3 positive'),
#		col = default.colours(2)
#		),
#	points = list(
#		pch = 22,
#		fill = default.colours(2),
#		col = 'black'
#	),
#	x = 0.02,
#	y = 0.9
#	),
#	filename = generate.filename('number_cells', 'samp3p1_left', 'pdf')
#	);
###
iseurat <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/samp3p1/objects/AllSample_GraphClust.seuset.rds');
mytype <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran//samp3p1/objects/CellType_Rename.txt', header = TRUE);
iseurat@meta.data$type <- mytype[match(rownames(iseurat@meta.data), mytype$Cell), ]$Cluster;
iseurat@meta.data$type <- gsub('Cell$', '', iseurat@meta.data$type);
iseurat@meta.data$type <- gsub('Macrophage|^DC', 'Monocytic', iseurat@meta.data$type);
imyexp <- data.frame(exp = iseurat@data['KLK3', ]);
imyexp$type <- iseurat@meta.data$type;
imyexp <- imyexp[grepl('LEFT', rownames(imyexp))&imyexp$exp>0, ];
imyexp$group <- 'epithelia';
imyexp[imyexp$type%in%c('Endothelia', 'Fibroblast'), ]$group <- 'stroma';
imyexp[!imyexp$type%in%c('Endothelia', 'Fibroblast', 'Epithelia'), ]$group <- 'immune';
#imyexp <- ddply(imyexp, 'group', numcolwise(mean));
pval1 <- scientific.notation(wilcox.test(imyexp[imyexp$group=='epithelia', ]$exp, imyexp[imyexp$group=='immune', ]$exp)$p.value);
pval1 <- ifelse(pval1==0, expression('< 2.2'%*%10^-16, ));
pval2 <- scientific.notation(wilcox.test(imyexp[imyexp$group=='stroma', ]$exp, imyexp[imyexp$group=='immune', ]$exp)$p.value)
create.boxplot(
	formula = exp~group,
	data = imyexp, 
	#add.stripplot = TRUE,
	xlab.label = 'Type',
	ylab.label = 'KLK3 abundance',
	add.text = TRUE,
	text.x = c(1.5, 2.5),
	text.y = 6.5,
	text.labels = c(pval1, pval2),
	text.cex = 1.5,
	xaxis.cex = 2,
	yaxis.cex =2,
	filename = generate.filename('exp_klk3p', 'type_left', 'pdf'),
	style = 'Nature'
	);
####

####
imyexp <- data.frame(exp = iseurat@data['KLK3', ]);
imyexp$type <- iseurat@meta.data$type;
imyexp <- imyexp[grepl('TUMOR', rownames(imyexp))&imyexp$exp>0, ];
imyexp$group <- 'epithelia';
imyexp[imyexp$type%in%c('Endothelia', 'Fibroblast'), ]$group <- 'stroma';
imyexp[!imyexp$type%in%c('Endothelia', 'Fibroblast', 'Epithelia'), ]$group <- 'immune';
#imyexp <- ddply(imyexp, 'group', numcolwise(mean));
pval1 <- scientific.notation(wilcox.test(imyexp[imyexp$group=='epithelia', ]$exp, imyexp[imyexp$group=='immune', ]$exp)$p.value);
pval1 <- ifelse(pval1==0, expression('< 2.2'%*%10^-16, ));
pval2 <- scientific.notation(wilcox.test(imyexp[imyexp$group=='stroma', ]$exp, imyexp[imyexp$group=='immune', ]$exp)$p.value)

create.boxplot(
	formula = exp~group,
	data = imyexp, 
	#add.stripplot = TRUE,
	xlab.label = 'Type',
	ylab.label = 'KLK3 abundance',
	add.text = TRUE,
	text.x = c(1.5, 2.5),
	text.y = 6.2,
	text.labels = c(pval1, pval2),
	text.cex = 1.5,
	xaxis.cex = 2,
	yaxis.cex =2,
	filename = generate.filename('exp_klk3p', 'type_tumor', 'pdf'),
	style = 'Nature',
	#width = 4.5,
	#height = 4.5
	);

###
summary(data.frame(t(iseurat@data[c('TP63', 'KRT14', 'KRT5'), iseurat@meta.data$orig.ident=='S6.TUMOR.M'&iseurat@meta.data$type=='Epithelia'])))
imyexp <- data.frame(t(as.matrix(iseurat@data[c('TP63', 'KRT14', 'KRT5'), ])));
imyexp$type <- iseurat@meta.data$type;
imyexp$orig.ident <- iseurat@meta.data$orig.ident;
imyexp <- imyexp[imyexp$type=='Epithelia', ];

to.plot <- data.frame(exp = iseurat@data['KLK3', ], type = iseurat@meta.data$type, 
	samp = iseurat@meta.data$orig.ident);
to.plot <- to.plot[to.plot$type!='Epithelia', ]
to.plot.m <- data.frame(reshape::cast(to.plot, type~samp, mean, value = 'exp', fill = 0));
to.plot.f <- data.frame(reshape::cast(to.plot[to.plot$exp>0, ], type~samp, length, value = 'exp', fill = 0));
to.plot.n <- data.frame(reshape::cast(to.plot, type~samp, length, value = 'exp', fill = 0));
#
to.plot.n$group <- 'all';
to.plot.f$group <- 'pos';
to.plot <- rbind(to.plot.f, to.plot.n);
to.plot.f[, 2:4] <- sapply(seq(3), function(x) 100*to.plot.f[, (x+1)]/to.plot.n[, (x+1)]);
to.plot.f[, 2:4] <- round(to.plot.f[, 2:4])
to.plot.f[6, 3] <- 'NA';

create.barplot(
	data = to.plot[, c('type', 'S5.LEFT.2', 'group')],
	formula = S5.LEFT.2~type,
	groups = to.plot$group,
	col = default.colours(2),
	style = 'Nature',
	#xaxis.lab = c('R', 'L', 'T'),
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 9000),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot[11:20, ]$S5.LEFT.2 + 500,
	text.labels = paste0(to.plot.f$S5.LEFT.2, '%'),
	text.cex = 1.5,
	text.fontface = 'plain',
	key = list(
	text = list(
		lab = c('All', 'KLK3 positive'),
		col = default.colours(2)
		),
	points = list(
		pch = 22,
		fill = default.colours(2),
		col = 'black'
	),
	x = 0.02,
	y = 0.9
	),
	filename = generate.filename('number_cells', 'samp3p1_left', 'pdf')
	);

create.barplot(
	data = to.plot[, c('type', 'S6.TUMOR.2', 'group')],
	formula = S6.TUMOR.M~type,
	groups = to.plot$group,
	col = default.colours(2),
	style = 'Nature',
	#xaxis.lab = c('R', 'L', 'T'),
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 9000),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot[11:20, ]$S6.TUMOR.M + 500,
	text.labels = paste0(to.plot$S6.TUMOR.M),
	text.cex = 1.5,
	text.fontface = 'plain',
	key = list(
	text = list(
		lab = c('All', 'KLK3 positive'),
		col = default.colours(2)
		),
	points = list(
		pch = 22,
		fill = default.colours(2),
		col = 'black'
	),
	x = 0.02,
	y = 0.9
	),
	filename = generate.filename('number_cells', 'samp3p1_right', 'pdf')
	);
####
create.barplot(
	data = to.plot[, c('type', 'S5.LEFT.2', 'group')],
	formula = S5.LEFT.2~type,
	groups = to.plot$group,
	col = default.colours(2),
	style = 'Nature',
	#xaxis.lab = c('R', 'L', 'T'),
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 9000),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot[11:20, ]$S5.LEFT.2 + 500,
	text.labels = paste0(to.plot$S5.LEFT.2),
	text.cex = 1.5,
	text.fontface = 'plain',
	key = list(
	text = list(
		lab = c('All', 'KLK3 positive'),
		col = default.colours(2)
		),
	points = list(
		pch = 22,
		fill = default.colours(2),
		col = 'black'
	),
	x = 0.02,
	y = 0.9
	),
	filename = generate.filename('number_cells', 'samp3p1_left', 'pdf')
	);

create.barplot(
	data = to.plot[, c('type', 'S6.TUMOR.M', 'group')],
	formula = S6.TUMOR.M~type,
	groups = to.plot$group,
	col = default.colours(2),
	style = 'Nature',
	#xaxis.lab = c('R', 'L', 'T'),
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 9000),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot[11:20, ]$S6.TUMOR.M + 500,
	text.labels = paste0(to.plot$S6.TUMOR.M),
	text.cex = 1.5,
	text.fontface = 'plain',
	key = list(
	text = list(
		lab = c('All', 'KLK3 positive'),
		col = default.colours(2)
		),
	points = list(
		pch = 22,
		fill = default.colours(2),
		col = 'black'
	),
	x = 0.02,
	y = 0.9
	),
	filename = generate.filename('number_cells', 'samp3p1_tumor', 'pdf')
	);
####
to.plot <- to.plot[to.plot$group=='pos', ];
create.barplot(
	data = to.plot[, c('type', 'S6.TUMOR.M', 'group')],
	formula = S6.TUMOR.M~type,
	style = 'Nature',
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 3000),
	yat = c(0, 1000, 2000, 3000),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot$S6.TUMOR.M + 200,
	text.labels = paste0(to.plot$S6.TUMOR.M),
	text.cex = 1.5,
	text.fontface = 'plain',
	filename = generate.filename('number_cells', 'samp3p1_tumor', 'pdf')
	);
####
create.barplot(
	data = to.plot[, c('type', 'S5.LEFT.2', 'group')],
	formula = S5.LEFT.2~type,
	style = 'Nature',
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 1500),
	yat = c(0, 500, 1000, 1500),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot$S5.LEFT.2 + 150,
	text.labels = paste0(to.plot$S5.LEFT.2),
	text.cex = 1.5,
	text.fontface = 'plain',
	filename = generate.filename('number_cells', 'samp3p1_left', 'pdf')
	);

create.barplot(
	data = to.plot[, c('type', 'S4.RIGHT.2', 'group')],
	formula = S4.RIGHT.2~type,
	style = 'Nature',
	xlab.label = '',
	ylab.labe = '# cells',
	xaxis.cex = 2,
	yaxis.cex = 2,
	xaxis.rot = 45,
	ylimits = c(0, 3),
	yat = c(0, 2),
	yaxis.lab = c(0, 2000),
	height = 4,
	width = 8,
	add.text = TRUE,
	text.x = seq(10),
	text.y = to.plot$S4.RIGHT.2 + .3,
	text.labels = paste0(to.plot$S4.RIGHT.2),
	text.cex = 1.5,
	text.fontface = 'plain',
	filename = generate.filename('number_cells', 'samp3p1_right', 'pdf')
	);
