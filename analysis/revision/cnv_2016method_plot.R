library(BoutrosLab.plotting.general);
library(gtools);
i <- 'JD1800154SL';
mycnv <- readRDS(paste0('2020-02-21_cnv_cells_', i, '.rds'));
mygene <- read.table('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/cnv/gencode_hg38_gene_pos_replaced_sorted_noHLA.txt', as.is = TRUE);
mygene <- mygene[mygene$V1%in%rownames(mycnv), ];
mygene <- mygene[mixedorder(mygene$V2), ];
mychr <- table(mygene$V2);
mychr <- mychr[mixedorder(names(mychr))];
gaps_col <- sapply(seq(length(mychr)), function(x) sum(mychr[1:x]));
mybreaks <- unique(c(seq(-0.5, 0.5, 0.01), seq(0.5, 3, 0.1), seq(-3, -0.5, 0.1)));
mycol <- colorRampPalette(c('blue', 'white', 'red'))(length(mybreaks));
hc <- hclust(dist(t(mycnv)));
myorder <- hc$order;
names(myorder) <- hc$labels;
myorder <- myorder[order(myorder)]
to.plot <- mycnv[match(mygene$V1, rownames(mycnv)), match(names(myorder), colnames(mycnv))];
p <- create.heatmap(
	x = to.plot,
	filename = NULL,
	colour.scheme = c('blue', 'white', 'red'),
	total.colours = length(mybreaks),
	at = mybreaks,
	cluster.dimensions = 'none',
	row.colour = 'black',
	col.lines = gaps_col,
	col.lwd = 2,
	grid.col = TRUE,
	force.grid.col = TRUE,
	);
png(generate.filename('plotcnv', i, 'png'), width = 10, height = 8.2, units = 'in', res = 300);
print(p);
dev.off();
