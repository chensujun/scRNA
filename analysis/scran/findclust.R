library(BoutrosLab.plotting.general);
library(scran);
library(scater);
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
source('~/svn/singleCell/myfunctions/find_clusters_snn.R');
####
setwd("/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/normalize_data");
csce <- readRDS(conf$scran_mnn);
csce <- run_snn(csce, 'all');
csce.20k <- run_plot_snn(csce, 'all', 20);

#mv 2019-03-19_all_mnn_snn.rds objects/
#mv 2019-03-21_all_mnn_snn_20.rds objects/
###
save.session.profile(generate.filename('session', 'find_cluster_snn', 'txt'));
