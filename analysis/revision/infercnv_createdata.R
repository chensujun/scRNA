library(BoutrosLab.plotting.general);
library(Seurat);
library(reshape2);
library(infercnv)
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv');
source('~/chensj/scRNA/script/myfunctions/get_dotplot_data.R');
source('~/chensj/scRNA/script/myfunctions/DotPlot_flip.R');
source('~/chensj/scRNA/script/myfunctions/utilities_internal.R');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
name <- 'all';
seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monolytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
seurat.all <- SetAllIdent(seurat.all, id = 'type');
mycount <- exp(seurat.all@data) - 1;
mycount <- mycount[!grepl('^HLA', rownames(mycount)), ];
for(i in unique(seurat.all@meta.data$orig.ident)){
	#assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
	mycounti <- mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])];
	write.table(mycounti, generate.filename('infercnv_count_norm', i, 'txt'), quote = FALSE, col.names = NA)
};
###
gene.pos <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/gencode_hg38_gene_pos_replaced_sorted.txt');
gene.pos <- gene.pos[gene.pos$V1%in%rownames(mycount), ];
write.table(gene.pos, 'gencode_hg38_gene_pos_replaced_sorted_noHLA.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
for(i in unique(seurat.all@meta.data$orig.ident)){
	itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, 'type_scnorm', drop = FALSE];
	rownames(itype) <- gsub('-', '.', rownames(itype))
	write.table(itype, generate.filename('infercnv_type', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE);
};
###
mycounti <- read.table('2020-02-14_infercnv_count_norm_JD1800153SL.txt');
name <- 'JD1800153SL';
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= mycounti,
    annotations_file= paste0('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/2020-02-14_infercnv_type_', name, '.txt'),
    delim="\t",
    gene_order_file="/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt",
    ref_group_names=c(NULL)
    );

out_dir=paste0(Sys.Date(), "_infercnv_ref_stroma");
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=out_dir,
    cluster_by_groups=F,
    plot_steps=F,
    mask_nonDE_genes = T,
    include.spike=T  # used for final scaling to fit range (0,2) centered at 1.
    );
save(infercnv_obj, file = paste0(out_dir, '/', 'infercnv_obj_output.rda'));
#### use the least variable cells from initial analysis as ref, refer https://science.sciencemag.org/content/sci/suppl/2016/04/07/352.6282.189.DC1/Tirosh.SM.pdf
setwd('ws');
geneMeans <- rowMeans(seurat.all@data);
mycells <- seurat.all@data[names(geneMeans[geneMeans>0.1]), ];
mygenes <- read.table('./gencode_hg38_gene_pos_replaced_sorted_noHLA.txt');
mycells <- mycells[rownames(mycells)%in%mygenes$V1, ];
mygenes <- mygenes[mygenes$V1%in%rownames(mycells), ];
mycells <- data.frame(t(scale(t(mycells), center = TRUE, scale = FALSE)));
mycells[mycells>3] <- 3;
mycells[mycells<(-3)] <- -3;
saveRDS(mycells, file = generate.filename('cnv0_cells', 'all', 'rds'));
####
mycnv0 <- data.frame(matrix(NA, nrow = nrow(mycells), ncol = ncol(mycells)));
colnames(mycnv0) <- colnames(mycells);
rownames(mycnv0) <- rownames(mycells);
###
for(i in seq(rownames(mycnv0))){
    if(i < 51){
        print(paste0(i, ':first 50'))
        mycnv0[i, ] <- colMeans(mycells[1:101, ])
    }else if(i>(nrow(mycnv0)-50)){
        print(paste0(i, ':last 50'))
        mycnv0[i, ] <- colMeans(mycells[(nrow(mycnv0)-100):nrow(mycnv0), ])
    }else{
        print(i)
        mycnv0[i, ] <- colMeans(mycells[(i-50):(i+50), ], )
    }
};
saveRDS(mycnv0, file = generate.filename('window_cnv0', 'all', 'rds'));
####
#hc <- hclust(dist(t(mycnv0)));
km <- kmeans(t(mycnv0), 16);
####
setwd('../cnv0/');
mycnv0 <- readRDS('../ws/2020-02-21_window_cnv0_all.rds');
cnv0.value <- colMeans(mycnv0^2);
cnv0.value <- cnv0.value[order(cnv0.value)];
save(cnv0.value, file = generate.filename('cnv0', 'level', 'rds'));
cnv0.s <- cnv0.value[gsub('\\.', '-', names(cnv0.value))%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='stroma', ])];
ref <- names(cnv0.s[1:1000]);
mycount <- exp(seurat.all@data) - 1;

for(i in unique(seurat.all@meta.data$orig.ident)){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ]), ref))];
    write.table(mycounti, generate.filename('infercnv_count_norm_cnv0', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i|rownames(seurat.all@meta.data)%in%ref, 'type_scnorm', drop = FALSE];
    itype[ref,] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_cnv0', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');

};
###
for(i in c('JD1800153SL', 'JD1800156SL')){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ]), ref))];
    write.table(mycounti, generate.filename('infercnv_count_norm_cnv0', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i|rownames(seurat.all@meta.data)%in%ref, 'type_scnorm', drop = FALSE];
    itype[ref,] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_cnv0', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');

};
###
####
setwd('ws');
mycnv0 <- readRDS('./2020-02-21_window_cnv0_all.rds');
colnames(mycnv0) <- gsub('\\.', '-', colnames(mycnv0));
cnv0.value <- colMeans(mycnv0^2);
cnv0.value <- cnv0.value[order(cnv0.value)];
cnv0.s <- cnv0.value[gsub('\\.', '-', names(cnv0.value))%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='stroma', ])];
ref <- gsub('\\.', '-', names(cnv0.s[1:1000]));
name <- 'JD1800159SL'
j <- 1;

mycnv <- data.frame(matrix(NA, nrow = nrow(mycnv0), ncol = ncol(mycnv0)));
colnames(mycnv) <- colnames(mycnv0);
rownames(mycnv) <- rownames(mycnv0);
for(j in seq(nrow(mycnv))){
    if(j%%10==0){
        print(j)
    }
    jcnv0 <- data.frame(t(mycnv0[j, ]));
    jcnv0.ref <- t(mycnv0[j, ref])
    max0 <- max(jcnv0.ref);
    min0 <- min(jcnv0.ref);
    jcnv0$cnv <- 0;
    jcnv0[jcnv0[, 1]>(max0+0.2), ]$cnv <- jcnv0[jcnv0[, 1]>(max0+0.2), 1] - max0;
    jcnv0[jcnv0[, 1]<(min0-0.2), ]$cnv <- jcnv0[jcnv0[, 1]<(min0-0.2), 1] - min0;
    mycnv[j, ] <- jcnv0$cnv
}
assign(paste0('cnv0.', name), mycnv0[, rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==name, ])]);
#####
##### using low cnv non-epithelia cells as reference
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
seurat.all@meta.data$cnv <- cnv0[rownames(seurat.all@meta.data)];
cnv0 <- seurat.all@meta.data[, c('cnv', 'type', 'orig.ident')];
cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv0 <- cnv0[order(cnv0$cnv), ];
cnv.ref <- data.frame(matrix(NA, nrow = 0, ncol = 3));
colnames(cnv.ref) <- colnames(cnv0);
for(i in unique(cnv0$type)){
    icnv <- cnv0[cnv0$type==i, ];
    cnv.ref <- rbind(cnv.ref, icnv[1:(nrow(icnv)/20), ]);
};
ref <- rownames(cnv.ref);

mycount <- exp(seurat.all@data) - 1;

for(i in unique(seurat.all@meta.data$orig.ident)){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ]), ref))];
    write.table(mycounti, generate.filename('infercnv_count_raw_ref1272', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i|rownames(seurat.all@meta.data)%in%ref, 'type_scnorm', drop = FALSE];
    itype[ref,] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_ref1272', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};
###
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
seurat.all@meta.data$cnv <- cnv0[rownames(seurat.all@meta.data)];
cnv0 <- seurat.all@meta.data[, c('cnv', 'type', 'orig.ident')];
cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv0 <- cnv0[order(cnv0$cnv), ];
cnv.ref <- data.frame(matrix(NA, nrow = 0, ncol = 3));
colnames(cnv.ref) <- colnames(cnv0);
for(i in unique(cnv0$type)){
    icnv <- cnv0[cnv0$type==i, ];
    cnv.ref <- rbind(cnv.ref, icnv[1:(nrow(icnv)/10), ]);
};
ref <- rownames(cnv.ref);

mycount <- exp(seurat.all@data) - 1;

for(i in unique(seurat.all@meta.data$orig.ident)){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ]), ref))];
    write.table(mycounti, generate.filename('infercnv_count_raw_ref1272', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i|rownames(seurat.all@meta.data)%in%ref, 'type_scnorm', drop = FALSE];
    itype[ref,] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_ref1272', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};
### use low cnv non-epithela
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
seurat.all@meta.data$cnv <- cnv0[rownames(seurat.all@meta.data)];
cnv0 <- seurat.all@meta.data[, c('cnv', 'type', 'orig.ident')];
#cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv0 <- cnv0[order(cnv0$cnv), ];
mycor <- readRDS('2020-03-03_infercnv_correlation_all.rds');
cnv0$cor <- mycor[rownames(cnv0)];
cnv0 <- cnv0[cnv0$cnv<0.04&cnv0$cor<0.4, ];

mycount <- exp(seurat.all@data) - 1;

for(i in unique(seurat.all@meta.data$orig.ident)[6:13]){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])))];
    write.table(mycounti, generate.filename('infercnv_count_raw_ref0', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, 'type_scnorm', drop = FALSE];
    itype[rownames(itype)%in%rownames(cnv0),] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_ref0', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};
###
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
seurat.all@meta.data$cnv <- cnv0[rownames(seurat.all@meta.data)];
cnv0 <- seurat.all@meta.data[, c('cnv', 'type', 'orig.ident')];
cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv0 <- cnv0[order(cnv0$cnv), ];
cnv.ref <- data.frame(matrix(NA, nrow = 0, ncol = 3));
colnames(cnv.ref) <- colnames(cnv0);
for(i in unique(cnv0$orig.ident)){
    icnv <- cnv0[cnv0$orig.ident==i, ];
    cnv.ref <- rbind(cnv.ref, icnv[1:(nrow(icnv)/10), ]);
};
ref <- rownames(cnv.ref);

mycount <- exp(seurat.all@data) - 1;

for(i in unique(seurat.all@meta.data$orig.ident)){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ]), ref))];
    write.table(mycounti, generate.filename('infercnv_count_raw_ref1272', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i|rownames(seurat.all@meta.data)%in%ref, 'type_scnorm', drop = FALSE];
    itype[ref,] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_ref1272', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};

###
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
seurat.all@meta.data$cnv <- cnv0[rownames(seurat.all@meta.data)];
cnv0 <- seurat.all@meta.data[, c('cnv', 'type', 'orig.ident')];
#cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv0 <- cnv0[order(cnv0$cnv), ];
mycor <- readRDS('2020-03-03_infercnv_correlation_all.rds');
cnv0$cor <- mycor[rownames(cnv0)];
cnv0 <- cnv0[cnv0$cnv<0.05&cnv0$cor<0.5, ];
cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
mycount <- exp(seurat.all@data) - 1;

for(i in unique(seurat.all@meta.data$orig.ident)){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    mycounti <- mycount[, colnames(mycount)%in%unique(c(rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])))];
    write.table(mycounti, generate.filename('infercnv_count_raw_nref0', i, 'txt'), quote = FALSE, col.names = NA, sep = '\t')
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, 'type_scnorm', drop = FALSE];
    itype[rownames(itype)%in%rownames(cnv0),] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_nref0', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};
###
setwd('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/nref0_epi');
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
seurat.all@meta.data$cnv <- cnv0[rownames(seurat.all@meta.data)];
cnv0 <- seurat.all@meta.data[, c('cnv', 'type', 'orig.ident')];
#cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv0 <- cnv0[order(cnv0$cnv), ];
mycor <- readRDS('/cluster/home/sujunc/chensj/scRNA/primary/scran/revision/cnv/raw/v0.8.2/2020-03-03_infercnv_correlation_all.rds');
cnv0$cor <- mycor[rownames(cnv0)];
cnv0 <- cnv0[cnv0$cnv<0.05&cnv0$cor<0.5, ];
#cnv0 <- cnv0[!cnv0$type%in%c('Basal/intermediate', 'Luminal'), ];

for(i in unique(seurat.all@meta.data$orig.ident)){
    #assign(paste0("mycount.", i), mycount[, colnames(mycount)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, ])]);
    itype <- seurat.all@meta.data[seurat.all@meta.data$orig.ident==i, 'type_scnorm', drop = FALSE];
    itype[rownames(itype)%in%rownames(cnv0),] <- 'ref';
    rownames(itype) <- gsub('-', '.', rownames(itype))
    write.table(itype, generate.filename('infercnv_type_nref0_epi', i, 'txt'), quote = FALSE, col.names = FALSE, row.names = TRUE, sep = '\t');
};

###
mycnv <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_infercnv_ref_stroma/infercnv.observations.txt', row.names = 1, as.is = TRUE, header=TRUE);
mycnv.ref <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_infercnv_ref_stroma/infercnv.references.txt', row.names = 1, as.is = TRUE, header=TRUE);
colnames(mycnv) <- gsub('\\.', '-', colnames(mycnv));
colnames(mycnv.ref) <- gsub('\\.', '-', colnames(mycnv.ref));
mycnv.all <- cbind(mycnv, mycnv.ref)
mycnv.all <- mycnv.all[, !colnames(mycnv.all)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type%in%c('Basal/intermediate', 'Luminal'), ])]
rv <- rowVars(as.matrix(mycnv.all));
mycnv.all <- mycnv.all[order(rv, decreasing=TRUE)[seq_len(1000)], ];
hc <- hclust(dist(t(mycnv.all)));
saveRDS(hc, file = generate.filename('hc', 'cnv_all', 'rds'));

png(generate.filename('plot_cnv', 'non-epithela', 'png'), width = 10, height = 8.2, units = 'in', res = 300);
pheatmap(mycnv.all, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE)
    #breaks = mybreaks, color = mycol, annotation_col = ann_col);
dev.off();

###
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/raw/v0.8.2');
name <- 'all';
mycnv <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_infercnv_ref_stroma/infercnv.observations.txt', row.names = 1, as.is = TRUE, header=TRUE);
mycnv.ref <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_infercnv_ref_stroma/infercnv.references.txt', row.names = 1, as.is = TRUE, header=TRUE);
colnames(mycnv) <- gsub('\\.', '-', colnames(mycnv));
colnames(mycnv.ref) <- gsub('\\.', '-', colnames(mycnv.ref));
cnv0 <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-29_infercnv_meansquare_all.rds');
names(cnv0) <- gsub('\\.', '-', names(cnv0));
cnv.all <- data.frame(cnv = cnv0, type = seurat.all@meta.data[names(cnv0), ]$type);
cnv.epi <- cnv.all[cnv.all$type%in%c('Basal/intermediate', 'Luminal'), ];
cnv.epi <- cnv.epi[order(-cnv.epi$cnv), ];
cnv.top <- cnv.epi[1:(nrow(cnv.epi)/20), ];
avg.top <- rowMeans(mycnv[, rownames(cnv.top)]);
mycor <- apply(mycnv, 2, function(x) cor(x, avg.top));
mycor.ref <- apply(mycnv.ref, 2, function(x) cor(x, avg.top));
mycor.all <- c(mycor, mycor.ref);
saveRDS(mycor.all, file = generate.filename('infercnv_correlation', name, 'rds'))
assign(paste0('mycor.', name), mycor.all);
assign(paste0('cnv.', name), cnv0);
to.plot <- data.frame(score = get(paste0('cnv.', name)), cor = get(paste0('mycor.', name)));
to.plot$col <- 'black';
to.plot[to.plot$score>0.04&to.plot$cor>0.4, ]$col <- 'red';
to.plot[to.plot$score<0.04&to.plot$cor<0.4, ]$col <- 'blue';
to.plot$type <- seurat.all@meta.data[match(gsub('\\.', '-', rownames(to.plot)), rownames(seurat.all@meta.data)), ]$type_scnorm;
mypch <- c(3,16,17)
names(mypch) <- c("epithelia", "Unknown", "stroma");
to.plot$pch <- mypch[to.plot$type]
create.scatterplot(
    formula = cor~score,
    data = to.plot,
    col = to.plot$col,
    pch = to.plot$pch,
    xlab.label = '',
    ylab.label = '',
    xaxis.fontface = 'plain', 
    yaxis.fontface = 'plain',
    abline.h = 0.4,
    abline.v = 0.04,
    style = 'Nature',
    filename = generate.filename('cnv_cor_score', name, 'pdf')
    )