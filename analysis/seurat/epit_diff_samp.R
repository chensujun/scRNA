library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(readxl);
library(reshape);
library(viridis);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/grep_term.R');
source('~/svn/singleCell/myfunctions/plot_enrich_list.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/epitumour");
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
name <- 'epi';
iseurat <- readRDS(conf$seurat_epi13);
iseurat <- SetIdent(iseurat, ident.use = iseurat@meta.data$fig.patient);
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
jm.out <- FindAllMarkers(iseurat, logfc.threshold = 0, only.pos = TRUE);
saveRDS(jm.out, file = generate.filename('diff_sample', name, 'rds'));
####
all.keg <- all.go <- list();
jm.out <- jm.out[jm.out$p_val_adj<0.05&abs(jm.out$avg_logFC)>log(1.5), ];
for(i in unique(jm.out$cluster)){
	idiff <- jm.out[jm.out$cluster==i, ];
	all.keg[[i]] <- test_enrich(annot.keg, idiff$gene)
	all.go[[i]] <- test_enrich(annot.go, idiff$gene)
};
save(all.keg, all.go, file = generate.filename('diff_sample_enrich', name, 'rda'));
to.plot <- grep_data(all.keg, Nmin = 10, Nmax = 500);
plot_enrich(to.plot, name, 12, 12, xrot = 45, cut1 = 0.25, dot.key = TRUE, colourkey = TRUE);

to.plot <- grep_data(all.go, Nmin = 10, Nmax = 500);
plot_enrich(to.plot, paste0(name, '_go'), width = 24, height = 24, xrot = 45, cut1 = 0.25, dot.key = TRUE, colourkey = TRUE);
####
iseurat <- ScaleData(iseurat);
plot.dat <- data.frame(
	"J159" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800159SL']),
	"J162" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800162SL']),
	"J174" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800174SL']),
	"J172" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800172SL']),
	"J173" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800173SL']),
	"J171" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800171SL']),
	"J175" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800175SL']),
	"J177" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800177SL']),
	"J176" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800176SL']),
	"J154" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800154SL']),
	"J155" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800155SL']),
	"J156" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800156SL']),
	"J153" = rowMeans(iseurat@data[, iseurat@meta.data$fig.patient=='JD1800153SL'])
	);
to.plot <- data.frame(t(scale(t(plot.dat[rownames(plot.dat)%in%jm.out$gene, ]))))
pdf(generate.filename('avgexp_diff', name, 'pdf'));
pheatmap(to.plot, show_rownames = FALSE, show_colnames = TRUE,
	color = colorRampPalette(c('blue', 'white', 'red'))(100)
	);
dev.off();