library(BoutrosLab.plotting.general);
library(viridis);
library(qusage);
library(gdata);
source('~/svn/singleCell/myfunctions/run_qusage.R');

annotype <- function(csce, name){
	#### using signatures from normal prostate scRNA-seq to differentiate stroma from epithelia
	myref <- read.xls('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/annotate_type/data/1-s2.0-S2211124718318771-mmc2.xlsx', as.is = TRUE);
	geneset.epi <- myref[, 2][sapply(myref[, 2], function(x) nchar(x))>0];
	geneset.stm <- myref[, 1][sapply(myref[, 1], function(x) nchar(x))>0];
	mygs <- list(stroma = geneset.stm, epithelia = geneset.epi)
	myresult <- run_qusage(csce, nm = 'scnorm', gs = mygs);
	mytype <- myresult[[2]];
	csce$type_scnorm <- mytype[match(csce$cluster, mytype$Cluster), ]$pathway.name;
	saveRDS(myresult, file = generate.filename(paste0('annotate_', name), 'result_scnorm', 'rds'));
	write.csv(data.frame(csce@colData), generate.filename(paste0('annotate_', name), 'colData_scnorm', 'csv'));
	write.table(data.frame(csce@colData)[, 'type_scnorm', drop = F], 'celltype_scran_scnorm.txt', quote = FALSE, row.names = TRUE, col.names = FALSE, sep = '\t');
	##### using signatures from normal prostate scRNA-seq to further subtype epithelia to luminal and basal
	mygs.epi <- list(BE = myref[, 3][sapply(myref[, 3], function(x) nchar(x))>0],
	        LE = myref[, 4][sapply(myref[, 4], function(x) nchar(x))>0],
	        OE = myref[, 5][sapply(myref[, 5], function(x) nchar(x))>0]
	        );
	##### 
	myresult <- run_qusage(csce, nm = 'scnorm_epi', gs = mygs.epi);
	mytype <- myresult[[2]];
	csce$type_epi <- mytype[match(csce$cluster, mytype$Cluster), ]$pathway.name;
	save(myresult, file = generate.filename(paste0('annotate_', name), 'result_scnorm_epi', 'rds'));
	#####
	# create a list object to store the markers
	# 
	#markers <- readRDS('2019-03-20_markers_celltype.rds');
	markers <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/annotate_type/2019-05-15_markers_celltype.rds');
	myresult <- run_qusage(csce, nm = 'manual', gs = markers, my.seed = 100);
	mytype <- myresult[[2]];
	csce$type <- mytype[match(csce$cluster, mytype$Cluster), ]$pathway.name;
	saveRDS(myresult, file = generate.filename('annotate_type', 'result_mannual', 'rds'));
	csce$type <- droplevels(csce$type);
	csce$type <- gsub('Dendritic_cells_cell201801', 'Dendritic', csce$type);
	csce$type <- gsub('Endothelial_cells_cell201801', 'Endothelial', csce$type);
	csce$type <- gsub('Fibroblasts_cell201801', 'Fibroblast', csce$type);
	csce$type <- gsub('Mast_cells_cell201801', 'Mast_cell', csce$type);
	csce$type <- gsub('_cell', '', csce$type);
	csce$type <- gsub('basals_profRen', 'Basal', csce$type);
	csce$type <- gsub('basal', 'Basal', csce$type);
	csce$type <- gsub('Basal$', 'Basal_intermediate', csce$type);
	csce$type <- gsub('B_Plasmas201801', 'B_Plasma', csce$type);
	csce$type <- gsub('Epithelial', 'Luminal', csce$type);
	csce$type <- gsub('Macrophage.*', 'Macrophage', csce$type);
	csce$type <- gsub('^B$', 'B/Plasma', csce$type);
	csce$type <- gsub('_', '/', csce$type);
	write.csv(data.frame(csce@colData), generate.filename(paste0('annotate_', name), 'colData_manual', 'csv'));
	saveRDS(csce, file = paste0('objects/', generate.filename(name, 'annotated', 'rds')));
	newcolData <- csce@colData;
	saveRDS(newcolData, file = generate.filename(paste0('annotate_', name), 'colData_manual', 'rds'));

	####
	fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
	        plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1));
	p1 <- plotTSNE(csce, colour_by="type_scnorm") + ggtitle("Stroma vs epithelial") +  
		geom_rug(colour = 'white', size = 2) + fontsize;
	p2 <- plotTSNE(csce, colour_by="type_epi") + ggtitle("Epithelial types") + 
		geom_rug(colour = 'white', size = 2) + fontsize;
	p3 <- plotTSNE(csce[, order(csce$type)], colour_by="type") + ggtitle("Type") + 
		geom_rug(colour = 'white', size = 2) + fontsize;

	pdf(generate.filename('plottsne', paste0(name, '_types'), 'pdf'), width = 18, height = 5)
	multiplot(p1, p2, p3, cols = 3)
	dev.off();
	#### return sce with annotation
	return(csce);
};
