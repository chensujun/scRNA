library(BoutrosLab.plotting.general);
library(mclust);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/multi/cluster');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/epitumour/objects/2019-08-27_sample13_epi.rds')
to.plot <- readRDS('../2020-05-26_cnv_score_pri_0518.rds');
mal.gmm <- readRDS('malignant.epithelial.meta.data.rds');
mal.gmm$cluster.ori <- seurat.all@meta.data[rownames(mal.gmm), ]$fig.cluster;
seurat.all@meta.data$col <- to.plot[gsub('-', '.', rownames(seurat.all@meta.data)), ]$col;
