library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(readxl);
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/macrophage');
source('/u/schen/svn/singleCell/myfunctions/grep_term.R');
source('~/svn/singleCell/myfunctions/plot_enrich_list.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
seurat.all <- readRDS(conf$sseurat_all);
iseurat <- readRDS(conf$seurat_macro);
#name <- 'myeloid';
name <- 'monocytic';
annot <- readRDS(conf$annot);
geneset <- read.table(conf$marker_geneset)
iseurat@meta.data$orig.ident <- seurat.all@meta.data[rownames(iseurat@meta.data), ]$orig.ident
iseurat@meta.data$status <- ifelse(iseurat@meta.data$orig.ident == 'JD1800172SL', 'LN', 'PRI');
iseurat@meta.data$idc <- ifelse(iseurat@meta.data$orig.ident%in%annot[annot$IDCP == 0, ]$sample, FALSE, TRUE);
iseurat@meta.data$GS <- annot[match(iseurat@meta.data$orig.ident, annot$sample), ]$pathological_gleason_grade;
iseurat@meta.data$pT <- annot[match(iseurat@meta.data$orig.ident, annot$sample), ]$pathological_t;
iseurat@meta.data$pN <- annot[match(iseurat@meta.data$orig.ident, annot$sample), ]$pathological_n;
iseurat@meta.data$pM <- annot[match(iseurat@meta.data$orig.ident, annot$sample), ]$pathological_m;
genes.marker <- c('status', 'idc', 'GS', 'pT', 'pN', 'pM');
plot.list <- list();
for(i in genes.marker){
	plot.list[[i]] <- DimPlot(iseurat, reduction.use = 'tsne', group.by = i, pt.size = 2, 
		no.axes = TRUE, vector.friendly = TRUE);
};
pdf(generate.filename('plotgene_clinical', name, 'pdf'), width = 15, height = 20);
plot_grid(plotlist = plot.list, ncol = 2);
dev.off();

#genes.marker <- c('HLA-DRA', 'HLA-DRB5', 'HLA-DRB1', 'FCGR3A', 'CD68', 'ANPEP', 'ITGAX', 'CD14', 'ITGAM', 'CD33')
#genes.marker <- c('HLA-DRA', 'FCGR3A', 'CD68', 'C1QA', 'ITGAX', 'CD14', 'ITGAM', 'CD163')
genes.marker <- c('HLA-DRA', 'FCGR3A', 'CD68', 'C1QA', 'CD86', 'MRC1', 'MRC2', 'CD163')
pdf(generate.filename('plottsne', paste0(name, '_marker'), 'pdf'), width = 8, height = 4.5);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 5, vector.friendly = TRUE, no.axes = TRUE, no.legend = TRUE, nCol = 4);
dev.off();
mytype <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data/objects/forCellPhone-Novel-final-version.txt', header = TRUE);
iseurat@meta.data$subtype <- droplevels(mytype[match(rownames(iseurat@meta.data), mytype$nGene), ]$type_epi);
iseurat@meta.data$subtype <- gsub('_[0-9]', '', iseurat@meta.data$subtype);
iseurat@meta.data$TAM1 <- colMeans(iseurat@data[rownames(iseurat@data)%in%geneset[geneset$V2=='TAM_M1_gene', ]$V1, ]);
iseurat@meta.data$TAM2 <- colMeans(iseurat@data[rownames(iseurat@data)%in%geneset[geneset$V2=='TAM_M2_gene', ]$V1, ]);
####
sig.m1 <- data.frame(read_excel('objects/1-s2.0-S0092867418307232-mmc4.xlsx', sheet = 2));
sig.m2s <- data.frame(read_excel('objects/1-s2.0-S0092867418307232-mmc4.xlsx', sheet = 3));
gene.m1 <- gsub(' .*', '', sig.m1$Gene);
gene.m2s <- gsub(' .*', '', sig.m2s$Gene);
gene.m1 <- c(intersect(gene.m1, rownames(iseurat@data)), 'NOS2', 'NOS1', 'NOS3', 'FCGR1A', 'FCGR1B', 'IL12A', 'TNF');
gene.m2 <- c(intersect(gene.m2s, rownames(iseurat@data)), 'ARG1', 'ARG2', 'FCGR2A', 'FCGR2B', 'FCGR2C', 'FCER2', 'PDCD1LG2', 'CD274', 'MRC1', 'IL1RN', 
	'IL1R2', 'CTSA', 'CTSC', 'WNT7B', 'FASLG');
iseurat@meta.data$TAM_m1 <- colMeans(iseurat@data[rownames(iseurat@data)%in%gene.m1, ]);
iseurat@meta.data$TAM_m2 <- colMeans(iseurat@data[rownames(iseurat@data)%in%gene.m2s, ]);


genes.marker <- c('TAM1', 'TAM2', 'TAM_m1', 'TAM_m2');
pdf(generate.filename('plottsne', paste0(name, '_geneset'), 'pdf'), width = 8, height = 2);
FeaturePlot(iseurat, reduction.use = 'tsne', features.plot = genes.marker, pt.size = 5, vector.friendly = TRUE, no.axes = TRUE, no.legend = TRUE, nCol = 4);
dev.off();
####
pdf(generate.filename('plotviolin', paste0(name, '_geneset'), 'pdf'), width = 8, height = 10);
VlnPlot(iseurat, features.plot = gene.m1, x.lab.rot = TRUE, point.size.use = 0.1, nCol = 4);
dev.off();

pdf(generate.filename('plotviolin', paste0(name, '_biocarta_NOS1'), 'pdf'), width = 8, height = 10);
VlnPlot(iseurat, features.plot = intersect(rownames(iseurat@data), a1$V1), x.lab.rot = TRUE, point.size.use = 0.1, nCol = 4);
dev.off();

#### genesets
to.plot <- iseurat@meta.data;
create.scatterplot(
	TAM1~TAM2,
	to.plot,
	groups = to.plot$subtype,
	col = default.colours(length(unique(to.plot$subtype))),
	cex = 0.5,
	add.xyline = TRUE,
	xyline.lty = 2,
	file = generate.filename('sig_tam', name, 'pdf'),
	style = 'Nature'
	);
create.scatterplot(
	TAM_m1~TAM_m2,
	to.plot,
	groups = to.plot$subtype,
	col = default.colours(length(unique(to.plot$subtype))),
	cex = 0.5,
	add.xyline = TRUE,
	xyline.lty = 2,
	file = generate.filename('sig_tam', name, 'pdf'),
	style = 'Nature'
	);

#### pathway
flist <- list.files('objects/Macrophage/2.MarkerGeneGOPath/GOAnalysis_result/', pattern = 'xlsx', full.names = TRUE);
for(i in grep('BP', flist)){
	ifile.b <- data.frame(read_excel(flist[i]));
	ifile.c <- data.frame(read_excel(flist[(i+1)]));
	ifile.m <- data.frame(read_excel(flist[(i+2)]));
	ifile <- rbind(cbind(ifile.b, type = 'BP'),
		cbind(ifile.c, type = 'CC'), cbind(ifile.m, type = 'MF'))	
	assign(paste0('mygo.', gsub('\\..*', '', gsub('.*/', '', flist[i]))), ifile)
};
####
flist <- list.files('objects/Macrophage/2.MarkerGeneGOPath/PathWayAnalysis_result/', pattern = 'xlsx', full.names = TRUE);
for(i in seq(length(flist))){
	ifile <- data.frame(read_excel(flist[i]));
	assign(paste0('mykegg.', gsub('\\..*', '', gsub('.*/', '', flist[i]))), ifile)
};
####
myego <- list(mykegg.0, mykegg.1, mykegg.2, mykegg.3, mykegg.4, mykegg.5, mykegg.6);
names(myego) <- paste0('C', seq(0, 6));
to.plot <- grep_data_keg(myego, cut.fdr = 0.05);
saveRDS(to.plot, generate.filename('enrich_cluster', name, 'rds'));
##
plot_enrich(to.plot, name, width = 8);

myego.go <- list(mygo.0, mygo.1, mygo.2, mygo.3, mygo.4, mygo.5, mygo.6);
names(myego.go) <- paste0('C', seq(0, 6));
to.plot.go <- grep_data_go(myego.go, cut.fdr = 0.05, Nmin = 10, Nmax = 500);
plot_enrich(to.plot.go, paste0(name, '_go'), width = 12, height = 12, spot.size.function = function(x) {ifelse(x>10, 10/3, abs(x)/3)});

