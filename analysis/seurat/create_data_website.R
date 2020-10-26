library(Seurat);
seurat.bak <- readRDS('2019-07-25_seurat_manual_all.rds');
perplexity = 30;
nPC = 8;
seurat.all <- RunTSNE(object = seurat.bak, reduction.use = "pca",  dims.use = 1:nPC, perplexity = perplexity,
		do.fast = TRUE, dim.embed = 3);
seurat.all@dr$tsne2 <- seurat.bak@dr$tsne;
seurat.all@scale.data <- NULL;
saveRDS(seurat.all, file = paste0(Sys.Date(), '_seurat_manual_all_reduced.rds'))
jm.out <- FindAllMarkers(seurat.all);