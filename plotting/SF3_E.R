library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape2);
library(pheatmap);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/SF3/');
name <- 'monolytic';
path.q <- read.table('~/chensj/scRNA/primary/scran/macrophage/objects/Macrophage/3.Qusage/QuSAGE_for_heatmap/KEGG_human_20190613_heatmap_input.txt');
path.q$max <- apply(path.q[, 1:7], 1, max);
path.q <- path.q[path.q$max>0.2|rowSums(path.q[, 1:7])>0.2, ];
to.plot <- data.frame(t(scale(t(path.q[, c(2, 4, 1, 3, 5, 6, 7)]))));
colnames(to.plot) <- gsub('Cluster_', 'C', colnames(to.plot));
#### Fig. S3E
pdf(generate.filename('kegg0.1', name, 'pdf'), height = 4, width = 6);
pheatmap(to.plot, cluster_cols = TRUE, border_color = NA, color = colorRampPalette(c('blue', 'white', 'red'))(100),
        breaks = seq(-max(to.plot), max(to.plot), length.out = 100));
dev.off();
###
write.csv(to.plot, generate.filename('plotdata', 'F_S3E', 'csv'));