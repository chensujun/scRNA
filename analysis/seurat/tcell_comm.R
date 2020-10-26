library(BoutrosLab.plotting.general);
library(Seurat);
library(gdata);
library(VennDiagram);
library(dendextend);
library(igraph);
library(qusage);
library(UpSetR);
setwd('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell');
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/run_qusage_seurat.R');
conf <- read.config.file('~/svn/singleCell/master_config_scRNA.R');
comm <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/tcell/raw_data/pvalues_0.05_communication.txt', header = TRUE, sep = '\t', as.is = TRUE);
comm.t <- comm[grepl('^CD', comm$a)&!comm$b%in%c('Endothelial6', 'Endothelial7'), ];
comm.t$type <- 'luminal';
comm.t[grepl('^CD', comm.t$b), ]$type <- 'tcell';
comm.t[grepl('Fibro', comm.t$b), ]$type <- 'fibro';
comm.t[grepl('TAM|Mono|DC', comm.t$b), ]$type <- 'myeloid';
comm.t[grepl('vCAF|Endo', comm.t$b), ]$type <- 'endo';

pair.f <- unique(comm.t[grepl('fibro', comm.t$type), ]$interacting_pair);
pair.e <- unique(comm.t[grepl('endo', comm.t$type), ]$interacting_pair);
pair.t <- unique(comm.t[grepl('tcell', comm.t$type), ]$interacting_pair);
pair.m <- unique(comm.t[grepl('myeloid', comm.t$type), ]$interacting_pair);
pair.l <- unique(comm.t[grepl('luminal', comm.t$type), ]$interacting_pair);

comm.t <- comm[grepl('^CD', comm$b)&!comm$a%in%c('Endothelial6', 'Endothelial7'), ];
comm.t$type <- 'luminal';
comm.t[grepl('^CD', comm.t$a), ]$type <- 'tcell';
comm.t[grepl('Fibro', comm.t$a), ]$type <- 'fibro';
comm.t[grepl('TAM|Mono|DC', comm.t$a), ]$type <- 'myeloid';
comm.t[grepl('vCAF|Endo', comm.t$a), ]$type <- 'endo';

pair.f <- unique(comm.t[grepl('fibro', comm.t$type), ]$interacting_pair);
pair.e <- unique(comm.t[grepl('endo', comm.t$type), ]$interacting_pair);
pair.t <- unique(comm.t[grepl('tcell', comm.t$type), ]$interacting_pair);
pair.m <- unique(comm.t[grepl('myeloid', comm.t$type), ]$interacting_pair);
pair.l <- unique(comm.t[grepl('luminal', comm.t$type), ]$interacting_pair);
