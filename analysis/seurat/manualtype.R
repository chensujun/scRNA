library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
source('~/svn/singleCell/myfunctions/cluster_annot_seurat.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
#seurat.all <- readRDS(conf$sseurat_all);
#out <- readRDS(conf$results_mnn);
#name <- 'all_val'
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data');
name <- 'all';
seurat.all <- readRDS(conf$sseurat_all);
newcoldata <- read.csv('2019-07-24_annotate_type_colData_scnorm.csv', row.names = 1, as.is = TRUE);
seurat.all@meta.data <- newcoldata;
#mytype <- read.table('objects/forCellPhone-Novel-final-version.txt', header = TRUE);
#mytype$type <- gsub('[0-9]+$|Cluster', '', mytype$type_epi);
#mytype$type <- gsub('vCAF', 'Endothelial', gsub('CellCycle|luminal', 'Luminal', gsub('DC.*|TAM.*|Mono.*', 'Myeloid', gsub('CD.*', 'T', mytype$type))));
##write.csv(annot.go[annot.go$name%in%c('androgen biosynthetic process', 'androgen metabolic process'), ], generate.filename('androgen', 'metabolic_synthetic', 'csv'))
#iseurat <- SubsetData(seurat.all, cells.use = mytype$nGene);
#iseurat@meta.data$type_manual <- droplevels(mytype[match(rownames(iseurat@meta.data), mytype$nGene), ]$type);
#iseurat <- SetAllIdent(iseurat, id = 'type_manual')
#seurat.all <- iseurat;
#saveRDS(seurat.all, file = paste0('objects/', generate.filename('seurat_manual', name, 'rds')))
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', gsub('Macrophage', 'Myeloid', seurat.all@meta.data$type));
seurat.all <- SetAllIdent(seurat.all, id = 'type');
p1 <- DimPlot(seurat.all, reduction.use = 'tsne', group.by = 'type', pt.size = 1, no.axes = TRUE);
p2 <- DimPlot(seurat.all, reduction.use = 'tsne', group.by = 'cluster', pt.size = 1, no.axes = TRUE);
pdf(generate.filename('typec', name, 'pdf'), width = 15);
plot_grid(p1, p2, rel_widths = c(1.15, 1));
dev.off();
saveRDS(seurat.all, file = paste0('objects/', generate.filename('seurat_manual', name, 'rds')))
