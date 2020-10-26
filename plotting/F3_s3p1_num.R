library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(readxl);
library(plyr);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/Figures/F3');
iseurat <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/samp3p1/objects/AllSample_GraphClust.seuset.rds');
mytype <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran//samp3p1/objects/CellType_Rename.txt', header = TRUE);
iseurat@meta.data$type <- mytype[match(rownames(iseurat@meta.data), mytype$Cell), ]$Cluster;
iseurat@meta.data$type <- gsub('Cell$', '', iseurat@meta.data$type);
iseurat@meta.data$type <- gsub('Macrophage|^DC', 'Monocytic', iseurat@meta.data$type);
name <- 'samp3p1'
to.plot <- data.frame(exp = iseurat@data['KLK3', ], type = iseurat@meta.data$type, 
	samp = iseurat@meta.data$orig.ident);
to.plot <- to.plot[to.plot$type!='Epithelia', ]
to.plot.m <- data.frame(reshape::cast(to.plot, type~samp, mean, value = 'exp', fill = 0));
to.plot.f <- data.frame(reshape::cast(to.plot[to.plot$exp>0, ], type~samp, length, value = 'exp', fill = 0));
to.plot.n <- data.frame(reshape::cast(to.plot, type~samp, length, value = 'exp', fill = 0));
rownames(to.plot.m) <- to.plot.m$type;
to.plot.m <- to.plot.m[, -1];
rownames(to.plot.n) <- to.plot.n$type;
to.plot.n <- to.plot.n[, -1];
rownames(to.plot.f) <- to.plot.f$type;
to.plot.f <- to.plot.f[, -1];
to.plot <- data.frame(type = rep(seq(nrow(to.plot.m)), times = 3), samp = rep(seq(3), each = nrow(to.plot.m)));
colours <- colorRampPalette(c('blue', 'red'))(11);

create.scatterplot(
	formula = type ~ samp,
	data = to.plot,
	cex = unlist(to.plot.f) / 500 + 0.01,
	col = colours[round(10*unlist(to.plot.m))+1],
	xlab.label = '',
	ylab.label = '',
	xaxis.lab = c('R', 'L', 'T'),
	yaxis.lab = rownames(to.plot.m),
	xat = seq(3),
	yat = seq(nrow(to.plot.m)),
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	xlimits = c(0.5, 3.5),
	width = 4,
	style = 'BoutrosLab',
	filename = generate.filename('numbmean_klk3', name, 'pdf')
	);
