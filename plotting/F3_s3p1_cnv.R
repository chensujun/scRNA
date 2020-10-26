library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F3');
iseurat <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/samp3p1/objects/AllSample_GraphClust.seuset.rds');
dat <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F3/CNV_cells.txt');
mytype <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran//samp3p1/objects/CellType_Rename.txt', header = TRUE);
iseurat@meta.data$type <- mytype[match(rownames(iseurat@meta.data), mytype$Cell), ]$Cluster;
iseurat@meta.data$type <- gsub('Cell$', '', iseurat@meta.data$type);
iseurat@meta.data$type <- gsub('Macrophage|^DC', 'Monocytic', iseurat@meta.data$type);
dat$type <- iseurat@meta.data[rownames(dat), ]$type;
dat$col <- 'black';
dat[dat$type=='Epithelia'&grepl('S6.TUMOR', rownames(dat)), ]$col <- 'red'
dat[dat$type=='Epithelia'&grepl('S5.LEFT', rownames(dat)), ]$col <- 'blue'
dat <- dat[order(dat$col), ]
create.scatterplot(
	formula = cnv_meanSquare ~ cnv_corr,
	data = dat,
	col = dat$col,
	xlab.label = 'Correlation',
	ylab.label = 'Ave.deviation',
	abline.h = 0.05,
	abline.v = 0.5,
	abline.lty = 2,
	xaxis.cex = 2,
	yaxis.cex = 2,
	key = list(
	text = list(
		lab = c('Epi-Tumor', 'Epi-LN.L', 'Others'),
		cex = 2,
		col = c('red', 'blue', 'black')
		),
	points = list(
		pch = 19,
		col = c('red', 'blue', 'black'),
		cex = 1.5
		),
	x = 0.01,
	y = 0.99
	),
	style = 'Nature',
	filename = generate.filename('cnv', 'samp3p1', 'pdf'),
	height = 4.5
	);

