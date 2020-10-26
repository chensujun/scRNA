library(BoutrosLab.plotting.general);
library(Seurat);
library(matrixStats);
library(pheatmap);
library(gtools);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
source('/cluster/projects/hansengroup/sujunc/scRNA/script/myfunctions/plot_functions.R');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
mygene <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt', as.is = TRUE);
dlist <- dir('./', '2020-03-09_infercnv_raw_nref0_spike_DE_g1_sd2');
seurat.all@meta.data$cnv_sep <- seurat.all@meta.data$cnv_cor <- NA;
seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monolytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
seurat.all <- SetAllIdent(seurat.all, id = 'type');
for(i in dlist){
	print(i);
	name <- gsub('2020-03-09_infercnv_raw_nref0_spike_DE_g1_sd2_', '', i);
	cnv.all <- readRDS(paste0('2020-03-09_infercnv_meansquare_raw_nref0_DE_', name, '.rds'));
	mycor.all <- readRDS(paste0('2020-03-09_infercnv_correlation_raw_nref0_DE_', name, '.rds'));
	seurat.all@meta.data[match(gsub('\\.', '-', names(cnv.all)), rownames(seurat.all@meta.data)), ]$cnv_sep <- cnv.all;
	seurat.all@meta.data[match(gsub('\\.', '-', names(mycor.all)), rownames(seurat.all@meta.data)), ]$cnv_cor <- mycor.all;
};
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/final');
saveRDS(seurat.all@meta.data[, c('cnv_sep', 'cnv_cor')], file = generate.filename('all', 'cnv_sep', 'rds'));
###
mywidth <- 16;
myheight <- 12;
myres <- 100;

fontsize2 <- theme(axis.text=element_text(size=48), axis.title.y=element_text(size=48, colour = 'black'), axis.title.x = element_text(size = 0),
	plot.title=element_text(size = 0),
    plot.margin=unit(c(.5,.5,1,2.5),"cm"), axis.text.x = element_text(angle = 45, hjust = 1, size = 48));

seurat.all@meta.data$type <- gsub('Basal/intermediate', 'Basal/int.', seurat.all@meta.data$type);
tiff(generate.filename('plotviolin', paste0(name, '_cnv'), 'tiff'), width = mywidth , height = myheight - 2, res = myres, units = 'in')
VlnPlot(seurat.all, features.plot = 'cnv_sep', x.lab.rot = TRUE, point.size.use = 0, group.by = 'type') + fontsize2 + ggtitle('') +
	#labs(x = '', y = expression(paste('(CNV-1)'^'2'))) +
	labs(x = '', y = expression('Ave. deviation')) +
	stat_summary(fun.y = median, geom = 'point', size = 25, colour = 'black', shape = 95) + 
	geom_hline(yintercept = 0.05, linetype = 'dashed', size = 2);
dev.off();
