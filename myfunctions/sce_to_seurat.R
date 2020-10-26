sce_to_seurat <- function(csce, name, newcolData){
	colnames(csce@reducedDims$MNN) <- paste0('dim', seq(50))
	colnames(csce@reducedDims$TSNE) <- paste0('tsne', seq(2))
	seurat.bak <- as.seurat(csce);
	###
	seurat.all <- CreateSeuratObject(raw.data = counts(csce));
	seurat.all <- NormalizeData(object = seurat.all, normalization.method = 'LogNormalize', scale.factor = 10000);
	seurat.all@data <- log(2^logcounts(csce));
	seurat.all@dr <- seurat.bak@dr;
	seurat.all@meta.data$orig.ident <- csce$batch;
	seurat.all@meta.data$type <- csce$type;
	seurat.all@ident <- as.factor(seurat.all@meta.data$type);
	names(seurat.all@ident) <- rownames(seurat.all@meta.data);
	newcolData <- readRDS(newcolData);
	seurat.all@meta.data$type <- newcolData$type;
	seurat.all@meta.data$type_scnorm <- newcolData$type_scnorm;
	seurat.all@meta.data$type_epi <- newcolData$type_epi;
	seurat.all@meta.data$cluster <- newcolData$cluster;
	seurat.all <- SetAllIdent(object = seurat.all, id = 'type');
	seurat.all <- ScaleData(seurat.all);
	saveRDS(seurat.all, generate.filename('seurat', name, 'rds'));
	return(seurat.all)
};