library(BoutrosLab.plotting.general);
library(Seurat);
library(pheatmap);
library(VennDiagram);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
source('~/chensj/scRNA/script/myfunctions/test_kegg.R');
source('~/chensj/scRNA/script/myfunctions/plot_enrich_manual.R');
seurat.all <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/normalize_data/objects/2019-07-25_seurat_manual_all.rds');
mycnv <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_infercnv_ref_stroma/infercnv.observations.txt');
colnames(mycnv) <- gsub('\\.', '-', colnames(mycnv));
myorder <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/2019-07-24_infercnv_ref_stroma/infercnv.observation_groupings.txt', header = TRUE);
#mygenes <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt');
mygenes <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/cnv/gencode_hg38_gene_pos_replaced_sorted.txt');
to.plot <- mycnv[rownames(mycnv)%in%mygenes[mygenes$V2=='chr6', ]$V1, match(rownames(myorder), colnames(mycnv))];
to.plot <- to.plot[, colnames(to.plot)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='Unknown', ])];
saveRDS(to.plot, generate.filename('plotdata', 'chr6_unknown', 'rds'))
mybreaks <- unique(c(seq(0, 0.5, 0.1), seq(0.5, 1.5, 0.01), seq(1.5, 2, 0.1)));
p <- create.heatmap(
	x = to.plot,
	filename = NULL,
	colour.scheme = c('blue', 'white', 'red'),
	total.colours = length(mybreaks),
	at = mybreaks,
	cluster.dimensions = 'none',
	row.colour = 'black',
	#col.lines = gaps_col,
	col.lwd = 2,
	#force.grid.col = TRUE,
	);
png(generate.filename('plotcnv', 'chr6_unknown', 'png'), width = 10, height = 8.2, units = 'in', res = 300);
print(p);
dev.off();

to.plot <- readRDS('2020-03-05_plotdata_chr6_unknown.rds');
to.plot <- t(to.plot);
to.plot <- to.plot[seq(nrow(to.plot), 1), ];
saveRDS(to.plot, generate.filename('plotdata_format', 'chr6_unknown', 'rds'));

to.plot <- readRDS('2020-03-05_plotdata_format_chr6_unknown.rds');
to.plot <- data.frame(to.plot);
to.plot <- to.plot[rownames(to.plot)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type_scnorm=='Unknown', ]), ];
mean.lum <- colMeans(to.plot[rownames(to.plot)%in%rownames(seurat.all@meta.data[seurat.all@meta.data$type=='Luminal', ]), ]);
mybreaks <- unique(c(seq(0, 0.5, 0.1), seq(0.5, 1.5, 0.01), seq(1.5, 2, 0.1)));
mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
ann_row <- data.frame(seurat.all@meta.data[, c('type_scnorm', 'type')]);
ann_row <- ann_row[match(rownames(to.plot), rownames(ann_row)), ];
ann_col <- data.frame(isHLA = ifelse(grepl('^HLA', colnames(to.plot)), 'HLA', 'Not'),
	group = ifelse(colnames(to.plot)%in%names(mean.lum[mean.lum<0.5]), 'Y', 'N'));
rownames(ann_col) <- colnames(to.plot);
ann_col <-data.frame(chr = rep('chr6', ncol(to.plot)), 
	loss = ifelse(colnames(to.plot)%in%names(mean.lum[mean.lum<0.5]), 'Y', 'N'));
rownames(ann_col) <- colnames(to.plot);
colnames(ann_row)[1] <- 'type_broad';
###
png(generate.filename('plotcnv', 'chr6_unknown', 'png'), width = 3, height = 8.2, units = 'in', res = 300);
pheatmap((to.plot), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
	breaks = mybreaks, color = mycol, annotation_col = ann_col, annotation_row = ann_row);
dev.off();
###

annot.keg <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
ego <- test_enrich(annot.keg, rownames(ann_col[ann_col$loss=='Y', ]))
ego.go <- test_enrich(annot.go, rownames(ann_col[ann_col$loss=='Y', ]))
save(ego, ego.go, file = generate.filename('enrich', 'chr6', 'rda'));
ego.go <- ego.go[ego.go$FDR<0.05, ];
ego.go <- ego.go[order(ego.go$FDR), ];
plot_enrich(list(go = ego.go), 'chr6_lossgenes', mybg = rownames(res), width = 15, height = 10,spot.size.function = function(x) {x}, key.sizes = c(2, 4));
#####
term.genes <- sapply(ego.go[1:15, ]$Intersect, function(x) strsplit(x, ',')[[1]]);
term.genes <- unique(unlist(term.genes));
all.genes <- rownames(ann_col[ann_col$loss=='Y', ])
###
cols=brewer.pal(9, "Set1");
venn.plot <- venn.diagram(list(loss = all.genes, term = term.genes ),
                          col = "transparent",
                          fill=rev(cols[1:2]),
                          alpha=0.5,
                          cex=1.4,
                          cat.fontface=3, 
                          filename=NULL,
                          cat.cex = 1.5,
                          cat.fontfamily = "serif")
 
pdf(generate.filename('loss_gene', 'term', 'pdf'))
grid.draw(venn.plot)
dev.off()
