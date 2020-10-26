library(Seurat);
setwd('/home/zhugh/sujun');
iseurat <- readRDS('~/180381/processing/epiTumor.basal.batchEffectRemovedMatirx.seuratClustered.rds');
iseurat@scale.data <- NULL;
saveRDS(iseurat, file = paste0(Sys.Date(), '_sample13_epi.rds'));
iseurat <- readRDS('~/180381/processing/crpc.4samples.epiTumor.basal.batchEffectRemovedMatirx.seuratClustered.rds');
iseurat@scale.data <- NULL;
saveRDS(iseurat, file = paste0(Sys.Date(), '_sample4crpc_epi.rds'));
#cp ~/180381/processing/epiTumor.basal.markers.3types.batchEffectRemovedMatirx.seuratClustered.rds ~/sujun/epi_markers_sample13.rds
#cp ~/180381/processing/crpc.4samples.epiTumor.basal.markers.batchEffectRemovedMatirx.seuratClustered.rds ~/sujun/epi_markers_crpc4.rds
