library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape2);
library(pheatmap);
library(readxl);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/Figures/SF3/');
name <- 'tcell';
fdir <- '~/chensj/scRNA/primary/scran/tcell/raw_data/TCell/5.Qusage/KEGG_human_20190613/';
flist <- list.files(fdir, 'xlsx');
tmp <- data.frame(read_excel(paste0(fdir, flist[1])));
read_n_sort <- function(x, y){
	ifile <- data.frame(read_excel(x));
	ifile <- ifile[match(y$GeneSet.name, ifile$GeneSet.name), ];
	return(ifile);
}
to.plot <- do.call(cbind, sapply(flist, function(x) read_n_sort(paste0(fdir, x), tmp)[, 2, drop = FALSE]));
to.plot <- data.frame(to.plot);
colnames(to.plot) <- gsub('GeneSetsInfo_|_KEGG_human_20190613.xlsx|luster_', '', flist);
rownames(to.plot) <- tmp$GeneSet.name;

dat.fdr <- do.call(cbind, sapply(flist, function(x) read_n_sort(paste0(fdir, x), tmp)[, 4, drop = FALSE]));
dat.fdr <- data.frame(dat.fdr);
colnames(dat.fdr) <- gsub('GeneSetsInfo_|_KEGG_human_20190613.xlsx|luster_', '', flist);
rownames(dat.fdr) <- tmp$GeneSet.name;

annot.keg <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
to.plot <- to.plot[rownames(to.plot)%in%annot.keg[grep('PATH:00|PATH:01', annot.keg$ont), ]$name, ];

####
mypath <- read.table('pathway_met.txt', sep = '\t');
to.plot <- to.plot[rownames(to.plot)%in%mypath$V1|grepl('D-Glutamine|Ubiquinone|2-Oxocarboxylic|Glycosaminoglycan biosynthesis - chondroitin sulfate|Glycosaminoglycan biosynthesis - keratan sulfate|Glycosphingolipid biosynthesis - lacto and neolacto', rownames(to.plot)), ]
plot.dat <- na.omit(data.frame(t(scale(t(to.plot)))));
plot.dat <- plot.dat[rownames(plot.dat)%in%rownames(dat.fdr[rowSums(dat.fdr<0.05)>0, ]), ];

pdf(generate.filename('kegg_met', name, 'pdf'), height = 10, width = 8);
pheatmap(plot.dat, cluster_cols = TRUE, border_color = NA, color = colorRampPalette(c('blue', 'white', 'red'))(100),
        breaks = seq(-max(plot.dat), max(plot.dat), length.out = 100));
dev.off();
###
write.csv(plot.dat, generate.filename('plotdata', 'F_S3H', 'csv'));