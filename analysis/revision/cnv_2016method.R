library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/ws/per');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
name <- 'all';
seurat.all@meta.data$type <- gsub('Macrophage|Myeloid', 'Monolytic', seurat.all@meta.data$type);
seurat.all@meta.data$type <- gsub('Myofibroblast', 'Fibroblast', seurat.all@meta.data$type);
####
mycnv0 <- readRDS('../2020-02-21_window_cnv0_all.rds');
colnames(mycnv0) <- gsub('\\.', '-', colnames(mycnv0));
cnv0.value <- colMeans(mycnv0^2);
cnv0.value <- cnv0.value[order(cnv0.value)];
hc <- hclust(dist(t(mycnv0[, rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='stroma', ])])));
saveRDS(hc, generate.filename('hclust_cnv0', 'stroma', 'rds'));

for(name in unique(seurat.all@meta.data$orig.ident)){
    geneMeans <- rowMeans(seurat.all@data[, seurat.all@meta.data$orig.ident==name]);
    mycells <- seurat.all@data[, seurat.all@meta.data$orig.ident==name][names(geneMeans[geneMeans>0.1]), ];
    mygenes <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt');
    mycells <- mycells[rownames(mycells)%in%mygenes$V1, ];
    mygenes <- mygenes[mygenes$V1%in%rownames(mycells), ];
    mycells <- data.frame(t(scale(t(mycells), center = TRUE, scale = FALSE)));
    mycells <- log2(exp(mycells));
    mycells[mycells>3] <- 3;
    mycells[mycells<(-3)] <- -3;
    saveRDS(mycells, file = generate.filename('cnv0_cells', name, 'rds'));
    ####
    mycnv0 <- data.frame(matrix(NA, nrow = nrow(mycells), ncol = ncol(mycells)));
    colnames(mycnv0) <- colnames(mycells);
    rownames(mycnv0) <- rownames(mycells);
    mycnv <- mycnv0;
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
    saveRDS(mycnv0, file = generate.filename('window_cnv0', name, 'rds'));
    ####
    colnames(mycnv0) <- gsub('\\.', '-', colnames(mycnv0));
    cnv0.value <- colMeans(mycnv0^2);
    cnv0.value <- cnv0.value[order(cnv0.value)];
    cnv0.s <- cnv0.value[gsub('\\.', '-', names(cnv0.value))%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='stroma', ])];
    hc <- hclust(dist(t(mycnv0)));
    saveRDS(hc, generate.filename('hclust', 'cnv0', 'rds'));

    ref <- cutree(hc, 5)[rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='stroma'&seurat.all@meta.data$orig.ident==name, ])];
    ref <- names(ref[ref==5])

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
    };
    saveRDS(mycnv, file = generate.filename('cnv_cells', name, 'rds'));
};
