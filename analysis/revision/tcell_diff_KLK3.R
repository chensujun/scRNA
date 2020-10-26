library(BoutrosLab.plotting.general);
library(Seurat);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/tcell/function');
seurat.t <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/tcell/raw_data/GraphClust.seuset.rds');
###
mydat <- seurat.t@data;
to.plot <- data.frame(all = rowMeans(mydat),
	pos = rowMeans(mydat[, mydat['KLK3', ]>0]),
	neg = rowMeans(mydat[, mydat['KLK3', ]==0]),
	pos2 = rowMeans(mydat[, mydat['KLK3', ]>2]));

nexp <- data.frame(all = rowSums(mydat>0),
	pos = rowSums(mydat[, mydat['KLK3', ]>0]>0),
	neg = rowSums(mydat[, mydat['KLK3', ]==0]>0),
	pos2 = rowSums(mydat[, mydat['KLK3', ]>2]>0));
###
jseurat <- SetIdent(seurat.t, ident.use = ifelse(seurat.t@data['KLK3', ]>0, 'pos', 'neg'));
jm.out <- FindMarkers(jseurat, ident.1 = 'pos', ident.2 = 'neg', logfc.threshold = 0);
save(jm.out, file = generate.filename('diff_klk3vsnon', 'all', 'rds'))

iseurat <- SubsetData(jseurat, cells.use = rownames(jseurat@meta.data[jseurat@meta.data$res.0.8%in%c(2,3,5), ]));
jm.out <- FindMarkers(iseurat, ident.1 = 'pos', ident.2 = 'neg', logfc.threshold = 0);
save(jm.out, file = generate.filename('diff_klk3vsnon', 'c235', 'rds'));
