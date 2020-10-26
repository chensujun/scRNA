library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(reshape);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/overview');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
name <- 'all';
annot <- readRDS(conf$annot);
seurat.all@meta.data$status <- ifelse(seurat.all@meta.data$orig.ident == 'JD1800172SL', 'LN', 'PRI');
seurat.all@meta.data$idc <- ifelse(seurat.all@meta.data$orig.ident%in%annot[annot$IDCP == 0, ]$sample, FALSE, TRUE);
seurat.all@meta.data$GS <- annot[match(seurat.all@meta.data$orig.ident, annot$sample), ]$pathological_gleason_grade;
seurat.all@meta.data$pT <- annot[match(seurat.all@meta.data$orig.ident, annot$sample), ]$pathological_t;
seurat.all@meta.data$pN <- annot[match(seurat.all@meta.data$orig.ident, annot$sample), ]$pathological_n;
seurat.all@meta.data$pM <- annot[match(seurat.all@meta.data$orig.ident, annot$sample), ]$pathological_m;

genes.marker <- c('status', 'idc', 'GS', 'pT', 'pN', 'pM');
plot.list <- list();
for(i in genes.marker){
	plot.list[[i]] <- DimPlot(seurat.all, reduction.use = 'tsne', group.by = i, pt.size = .5, 
		no.axes = TRUE, vector.friendly = TRUE);
};
pdf(generate.filename('plotgene_clinical', name, 'pdf'), width = 15, height = 20);
plot_grid(plotlist = plot.list, ncol = 2);
dev.off();
####
to.plot <- (cor(t(table(seurat.all@meta.data$type, seurat.all@meta.data$orig.ident))));
pdf(generate.filename('correlation', 'number_all', 'pdf'));
pheatmap(to.plot, display_numbers = TRUE, color = colorRampPalette(c('blue', 'white', 'red'))(100),
	breaks = seq(-1, 1, length.out = 100));
dev.off();

mytype <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data/objects/forCellPhone-Novel-final-version.txt', header = TRUE);
seurat.all@meta.data$type.cluster <- mytype[match(rownames(seurat.all@meta.data), mytype$nGene), ]$type_epi;
to.plot <- (cor(t(table(seurat.all@meta.data$type.cluster, seurat.all@meta.data$orig.ident))));
pdf(generate.filename('correlation', 'number_all_subtype', 'pdf'), width = 7.5);
pheatmap(to.plot, display_numbers = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(100),
	breaks = seq(-1, 1, length.out = 100));
dev.off();

seurat.all@meta.data$typec <- paste0(seurat.all@meta.data$type, '_', seurat.all@meta.data$cluster)
to.plot <- (cor(t(table(seurat.all@meta.data$typec, seurat.all@meta.data$orig.ident))));
pdf(generate.filename('correlation', 'number_all_typecluter', 'pdf'), width = 7.5);
pheatmap(to.plot, display_numbers = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(100),
	breaks = seq(-1, 1, length.out = 100));
dev.off();
#####
comm <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pvalues_0.05_communication.txt', header = TRUE, sep = '\t', as.is = TRUE);
comm.t <- comm[grepl('Basal', comm$a)&grepl('TAM|^CD', comm$b), ];
comm.t$type <- 'luminal';
comm.t[grepl('^CD', comm.t$b), ]$type <- 'tcell';
comm.t[grepl('TAM', comm.t$b), ]$type <- 'myeloid';
comm.e <- comm[grepl('luminal', comm$a)&grepl('TAM|^CD', comm$b), ];
comm.e$type <- 'luminal';
comm.e[grepl('^CD', comm.e$b), ]$type <- 'tcell';
comm.e[grepl('TAM', comm.e$b), ]$type <- 'myeloid';
####
comm$fdr <- p.adjust(comm$P.Value)
comm.t <- comm[comm$fdr<0.05&comm$significant_mean>1, ];
comm.t$type <- 'luminal';
comm.t[grepl('^CD', comm.t$b), ]$type <- 'tcell';
comm.t[grepl('Fibro', comm.t$b), ]$type <- 'fibro';
comm.t[grepl('TAM|Mono|DC', comm.t$b), ]$type <- 'myeloid';
comm.t[grepl('vCAF|Endo', comm.t$b), ]$type <- 'endo';
comm.t[grepl('Basal', comm.t$b), ]$type <- 'basal';
comm.t[grepl('Mast', comm.t$b), ]$type <- 'mast';

comm.t$type.a <- 'luminal';
comm.t[grepl('^CD', comm.t$a), ]$type.a <- 'tcell';
comm.t[grepl('Fibro', comm.t$a), ]$type.a <- 'fibro';
comm.t[grepl('TAM|Mono|DC', comm.t$a), ]$type.a <- 'myeloid';
comm.t[grepl('vCAF|Endo', comm.t$a), ]$type.a <- 'endo';
comm.t[grepl('Basal', comm.t$a), ]$type.a <- 'basal';
comm.t[grepl('Mast', comm.t$a), ]$type.a <- 'mast';
comm.t <- comm.t[comm.t$type%in%c('luminal', 'basal')&!comm.t$type.a%in%c('luminal', 'basal'), ];


to.plot.p <- cast(comm.t, interacting_pair~type.a+type, mean, value = 'P.Value');
to.plot.m <- cast(comm.t, interacting_pair~type.a+type, mean, value = 'significant_mean');
rownames(to.plot.p) <- rownames(to.plot.m) <- to.plot.m$interacting_pair
to.plot.p <- -log10(to.plot.p[, -1]);
to.plot.m <- to.plot.m[, -1];
to.plot.p[sapply(to.plot.p, is.infinite)] <- 22;
to.plot.m <- to.plot.m[, !colnames(to.plot.m)%in%paste0(unique(comm.t$type), '_', unique(comm.t$type))];
to.plot.m$max <- apply(to.plot.m[, 1:10], 1, function(x) max(na.omit(x)));
to.plot.m$number <- apply(to.plot.m[, 1:10], 1, function(x) length(na.omit(x)));
####
npair <- apply(to.plot.m, 2, function(x) length(na.omit(x)));
names(npair) <- colnames(to.plot.m);
npair.u <- vector();
type <- names(table(comm.t$type));
for(j in seq(6)){
	for(i in seq(7-j)){
		npair.u[paste0(type[i], '_', type[i+j])] <- npair[paste0(type[i], '_', type[i+j])] + npair[paste0(type[i+j], '_', type[i])];
	}
}
####
#iseurat <- SubsetData(seurat.all, cells.use = rownames(seurat.all@meta.data[!seurat.all@meta.data$type%in%c('Luminal', 'Basal/intermediate'), ]));
to.plot.m <- to.plot.m[to.plot.m$number>1, colnames(to.plot.m)%in%colnames(to.plot.p)];
to.plot.p <- to.plot.p[rownames(to.plot.m), colnames(to.plot.m)];

dotkeymax <- 3;
colourkey <- FALSE;
cut1 <- 0.25;
xrot <- 90;
plot_comm <- function(to.plot.m, to.plot.p, name, width = 5, height = 20,
  dotkeymax = 3, colourkey = FALSE, cut1 = 0.25, xrot = 90,
  spot.colour.function = function(x) {colours = rep('black', length(unlist(to.plot.m))); return(colours)},
  spot.size.function =  function(x) {abs(x)}){
  #spot.colour.function <- function(x) {
  #    mycol <- colorRampPalette(c('black', 'yellow', 'red'))(length(unique(unlist(to.plot.m))))
  #    names(mycol) <- unique(unlist(to.plot.m))
  #    colours <- mycol[x]
  #    return(colours);
  #};
  dot.key <- list(
      # indicate which side of the plot the legend will appear
      space = "right",
      between.rows = 10,
      points = list(
              cex = spot.size.function(c(1, dotkeymax)),
              col = spot.colour.function(c(1, dotkeymax)),
              pch = 19
              ),
      # dot labels
      text = list(
              lab = as.character(c(1, dotkeymax)),
              cex = 1.5,
              adj = 1.0,
              fontface = "bold"
              )
      );

  
  create.dotmap(
        file = generate.filename('dotmap_comm', name, 'pdf'),
        x = to.plot.m,
        xaxis.cex = 1,
        yaxis.cex = 1.2,
        left.padding = 0,
        bottom.padding = 4,
        # use specified spot size and colour functions
        spot.size.function = spot.size.function,
        spot.colour.function = spot.colour.function,
        # create a legend matching the dot sizes
        key = dot.key,
        key.top = 1,
        xaxis.lab = gsub('JD1800|SL', '', colnames(to.plot.m)),
        yaxis.lab = rownames(to.plot.m),
        xaxis.rot = xrot,
        pch = 21,
        pch.border.col = 'transparent',
        # add the background
        bg.data = to.plot.p,
        # add a colourkey
        colourkey = colourkey,
        colour.scheme = c("black", "white"),
        total.colour = 5,
        bg.alpha = 1,
        at = c(0, -log10(cut1), -log10(0.01), 5),
        colourkey.labels.at = c(0, -log10(cut1), 2, 5),
        colourkey.labels = c(1, cut1, expression(10^-2), expression(''<=10^-5)),
        width = width,
        height = height,
        na.spot.size = 3,
        add.grid = TRUE,
        col.lwd = 1,
        style = 'Nature',
        col.colour = 'black',
        row.colour = 'black',
        );
};

plot_comm(to.plot.m, to.plot.p, name, width = 10, height = 20)