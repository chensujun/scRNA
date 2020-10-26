library(BoutrosLab.plotting.general);
library(viridis);
library(qusage);
library(gdata);
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');

find_cluster_seurat <- function(seurat.all, out, name, perplexity = 30, pcs.compute = 20, nPC = 8){
	seurat.all@meta.data$orig.cluster <- seurat.all@meta.data$cluster;
	cor.exp <- tcrossprod(out$rotation, out$corrected[colnames(seurat.all@data), ]);
	cor.exp <- cor.exp[rownames(seurat.all@data), colnames(seurat.all@data)];
	seurat.all@scale.data <- cor.exp;
	tmp <- seurat.all;
	seurat.all@meta.data$orig.type <- seurat.all@meta.data$type;
	### scale and center the mnn data
	seurat.all <-ScaleData(object = seurat.all, data.use=seurat.all@scale.data);
	### assign scaled data to data slot
	seurat.all@data <- seurat.all@scale.data;
	seurat.all <- FindVariableGenes(object= seurat.all,  mean.function = ExpMean, dispersion.function = LogVMR, 
		x.low.cutoff = -1, x.high.cutoff = Inf, y.cutoff = 1, do.plot = TRUE);
	#seurat.all <- RunPCA(object =seurat.all, pc.genes = names(gene.feature), pc.compute = 20,
	#	do.print = TRUE, pcs.print = 1:5, genes.print = 5, seed.use = 42);
	seurat.all <- RunPCA(object =seurat.all, pcs.compute = pcs.compute,
		do.print = TRUE, pcs.print = 1:5, genes.print = 5, seed.use = 42);
	#use JackStraw to find number of PCs
	#seurat.all <- JackStraw(seurat.all);
	#seurat.all <- JackStrawPlot(seurat.all, PCs = 1:20);
	#### grab significant PCs
	#pval.pc <- as.data.frame(seurat.all@dr$pca@jackstraw@overall.p.values);
	#pval.pc$padj <- p.adjust(pval.pc$Score);
	#nPC <- max(pval.pc[pval.pc$padj<0.01, ]$PC);
	pdf(generate.filename('Elbow', name, 'pdf'))
	PCElbowPlot(seurat.all, num.pc = pcs.compute);
	dev.off();
	seurat.all <- ProjectPCA(object = seurat.all, do.print = FALSE);
	seurat.all <- FindClusters(object = seurat.all, reduction.type = "pca", 
		dims.use = 1:nPC, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE);
	seurat.all <- RunTSNE(object = seurat.all, reduction.use = "pca",  dims.use = 1:nPC, perplexity = perplexity,
		do.fast = TRUE);
	seurat.all@meta.data$cluster <- as.numeric(as.vector(seurat.all@ident)) + 1;
	seurat.all@data <- tmp@data;
	#####
	markers <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/annotate_type/2019-05-15_markers_celltype.rds');
	my.seed <- 100;
	myresult <- run_qusage(seurat.all, nm = 'manual', gs = markers, my.seed = my.seed);
	mytype <- myresult[[2]]
	#seurat.all@meta.data$type_ntop <- droplevels(mytype[match(seurat.all@meta.data$cluster, mytype$Cluster), ]$pathway.name);
	seurat.all@meta.data$type <- droplevels(mytype[match(seurat.all@meta.data$cluster, mytype$Cluster), ]$pathway.name);
    seurat.all@meta.data$type <- gsub('Dendritic_cells_cell201801', 'Dendritic', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Endothelial_cells_cell201801', 'Endothelial', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Fibroblasts_cell201801', 'Fibroblast', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Mast_cells_cell20180.*', 'Mast_cell', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('_cell', '', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('basals_profRen', 'Basal', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('basal', 'Basal', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Basal$', 'Basal_intermediate', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('B_Plasmas201801', 'B_Plasma', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Epithelial', 'Luminal', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Macrophage.*', 'Macrophage', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('^B$', 'B/Plasma', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('_', '/', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Basal/full', 'Basal/intermediate', seurat.all@meta.data$type);
    seurat.all@meta.data$type <- gsub('Myocytes201801', 'Myocytes', seurat.all@meta.data$type);

	p1 <- DimPlot(seurat.all, reduction.use = 'tsne', group.by = 'type', pt.size = 1, no.axes = TRUE);
	p2 <- DimPlot(seurat.all, reduction.use = 'tsne', group.by = 'cluster', pt.size = 1, no.axes = TRUE);
	pdf(generate.filename('typec', name, 'pdf'), width = 15);
	plot_grid(p1, p2, rel_widths = c(1.15, 1));
	dev.off();
	####
	genes.major <- c('CD3D', 'CD79A', 'C1QA', 'CD163', 'EPCAM', 'PECAM1', 'VIM', 'DCN', 'S100A9');
	pdf(generate.filename('marker_major', name, 'pdf'));
	FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = genes.major, pt.size = 2,
		cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE)
	dev.off();

	####
	seurat.all <- SetAllIdent(object = seurat.all, id = 'type');
	return(seurat.all);
};


run_tsne_seurat <- function(seurat.all, out, name, perplexity = 30, pcs.compute = 20, nPC = 8){
	seurat.all@meta.data$orig.cluster <- seurat.all@meta.data$cluster;
	cor.exp <- tcrossprod(out$rotation, out$corrected[colnames(seurat.all@data), ]);
	cor.exp <- cor.exp[rownames(seurat.all@data), colnames(seurat.all@data)];
	seurat.all@scale.data <- cor.exp;
	tmp <- seurat.all;
	seurat.all@meta.data$orig.type <- seurat.all@meta.data$type;
	### scale and center the mnn data
	seurat.all <-ScaleData(object = seurat.all, data.use=seurat.all@scale.data);
	### assign scaled data to data slot
	seurat.all@data <- seurat.all@scale.data;
	seurat.all <- FindVariableGenes(object= seurat.all,  mean.function = ExpMean, dispersion.function = LogVMR, 
		x.low.cutoff = -1, x.high.cutoff = Inf, y.cutoff = 1, do.plot = TRUE);
	#seurat.all <- RunPCA(object =seurat.all, pc.genes = names(gene.feature), pc.compute = 20,
	#	do.print = TRUE, pcs.print = 1:5, genes.print = 5, seed.use = 42);
	seurat.all <- RunPCA(object =seurat.all, pcs.compute = pcs.compute,
		do.print = TRUE, pcs.print = 1:5, genes.print = 5, seed.use = 42);
	seurat.all <- ProjectPCA(object = seurat.all, do.print = FALSE);
	seurat.all <- RunTSNE(object = seurat.all, reduction.use = "pca",  dims.use = 1:nPC, perplexity = perplexity,
		do.fast = TRUE);
	seurat.all@meta.data$cluster <- as.numeric(as.vector(seurat.all@ident)) + 1;
	seurat.all@data <- tmp@data;
	#####
	return(seurat.all);
};

