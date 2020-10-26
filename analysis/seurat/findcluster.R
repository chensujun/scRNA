library(BoutrosLab.plotting.general);
library(Seurat);
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
#seurat.all <- readRDS(conf$sseurat_all);
#out <- readRDS(conf$results_mnn);
#name <- 'all_val'
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data');
seurat.all <- readRDS(conf$sseurat_all);
out <- readRDS(conf$results_mnn);
name <- 'all'
seurat.all <- find_cluster_seurat(seurat.all, out, name);
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat', name, 'rds')));
saveRDS(seurat.all@meta.data, file = generate.filename('seurat_coldata', name, 'rds'));
seurat.all@scale.data <- NULL;
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat_reduce', name, 'rds')));
####
pdf(generate.filename('marker_myocyte', name, 'pdf'));
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = markers[['Myocytes_cell201801']], pt.size = 1,
		cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE)
dev.off();

pdf(generate.filename('marker_bcell', name, 'pdf'));
FeaturePlot(seurat.all, reduction.use = 'tsne', pt.size = 1, cols.use = c('grey', 'red'), vector.friendly = TRUE, 
		features.plot = intersect(markers[['B_cell_Lymph_node']], rownames(seurat.all@data)))
dev.off();

pdf(generate.filename('marker_luminal', name, 'pdf'));
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = markers[['Luminal_cell']], pt.size = 1,
		cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE)
dev.off();

#### the memory couldn't handle 40k+ cells at once, need run the find_cluster_seurat function line by line
####### save the files after scale, then re-start R and read back in the object to FindVariableGenes 
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm');
seurat.all <- readRDS(conf$seurat_all11);
out <- readRDS(conf$results_mnn11);
name <- 'all11'
seurat.all <- find_cluster_seurat(seurat.all, out, name);
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat', name, 'rds')));
seurat.all@scale.data <- NULL;
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat_reduce', name, 'rds')));
saveRDS(seurat.all@meta.data, file = generate.filename('seurat_coldata', name, 'rds'));
saveRDS(seurat.all@dr, generate.filename('seurat_dr', name, 'rds'));
####


setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm');
seurat.all <- readRDS(conf$seurat_crpc);
out <- readRDS(conf$results_mnn11);
name <- 'crpc'
seurat.all <- find_cluster_seurat(seurat.all, out, name);
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat', name, 'rds')));
seurat.all@scale.data <- NULL;
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat_reduce', name, 'rds')));
saveRDS(seurat.all@meta.data, file = generate.filename('seurat_coldata', name, 'rds'));
saveRDS(seurat.all@dr, generate.filename('seurat_dr', name, 'rds'));
####
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm');
seurat.all <- readRDS(conf$seurat_norm);
out <- readRDS(conf$results_mnn11);
name <- 'norm'
seurat.all <- find_cluster_seurat(seurat.all, out, name);
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat', name, 'rds')));
seurat.all@scale.data <- NULL;
saveRDS(seurat.all, paste0('objects/', generate.filename('seurat_reduce', name, 'rds')));
saveRDS(seurat.all@meta.data, file = generate.filename('seurat_coldata', name, 'rds'));
saveRDS(seurat.all@dr, generate.filename('seurat_dr', name, 'rds'));
####
