library(BoutrosLab.plotting.general);
library(Seurat);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/create.heatmap_modi_key.R');
name <- 'vCAF';
mscore <- read.table('../endo/objects/Manual/Manual.ModuleExp.txt', row.names = 1, sep = '\t', header = TRUE);
mgene <- read.table('../endo/objects/Manual/Manual.module.txt', header = TRUE);
rownames(mscore) <- gsub(' ', '_', rownames(mscore));
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');

mscore$maxid <- apply(mscore, 1, function(x) grep(max(x), x));

module.e <- mscore[rowSums(mscore[, c(1, 2)]>0)==2&mscore$maxid%in%c(1, 2), ];
module.e <- module.e[apply(module.e[, -7], 1, function(x) min(x[c(1, 2)])>max(x[-c(1, 2)])), ];

module.v1 <- mscore[mscore$Cluster2>0&mscore$maxid==3, ];
module.v1 <- module.v1[module.v1$Cluster2>0.1&rowMeans(module.v1[, 1:2])<0, ];
module.v1 <- module.v1[rowMeans(module.v1[, c(4:6)])<0, ];

module.v2 <- mscore[rowSums(mscore[, c(4, 5)]>0)==2&mscore$maxid%in%c(4, 5), ];
module.v2 <- module.v2[apply(module.v2[, -7], 1, function(x) min(x[c(4, 5)])>max(x[-c(4, 5)])), ];

module.v3 <- mscore[mscore$Cluster5>0&mscore$maxid==6, ];
module.v3 <- module.v3[module.v3$Cluster5>0.1&rowMeans(module.v3[, 1:5])<0, ];
module.v3 <- module.v1[rowMeans(module.v1[, c(4:6)])<0, ];
gene.v3 <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.v3)), ]$id);
keg.v3 <- test_enrich(annot.keg, gene.v3, mybg = NULL);
go.v3 <- test_enrich(annot.go, gene.v3, mybg = NULL);

to.plot <- mscore[c(rownames(module.e), rownames(module.v1), rownames(module.v2),rownames(module.v3)), -7];
to.plot <- data.frame(t(scale(t(to.plot))));
create.heatmap(
	to.plot,
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	colour.scheme = c('blue', 'white', 'red'),
	print.colour.key = TRUE,
	colourkey.cex = 1,
	yaxis.lab = NA,
	xaxis.lab = gsub('Cluster', '', colnames(to.plot)),
	xaxis.rot = 0,
	yaxis.fontface = 'plain',
	yaxis.cex = 1,
	xaxis.fontface = 'plain',
	xaxis.cex = 1,
	filename = generate.filename('modules_unique', name, 'pdf'),
	width = 4
	);

test_module <- function(mymodule, mgene, annot.keg, annot.go){
	gene.mymodule <- as.vector(mgene[mgene$module%in%gsub('Module_', '', mymodule), ]$id);
	keg.mymodule <- test_enrich(annot.keg, gene.mymodule, mybg = NULL);
	go.mymodule <- test_enrich(annot.go, gene.mymodule, mybg = NULL);
	myresult <- list(keg = keg.mymodule, go = go.mymodule);
	return(myresult);
};

keg.all <- go.all <- list();
for(i in rownames(to.plot)){
	print(i);
	iresult <- test_module(i, mgene, annot.keg, annot.go);
	keg.all[[i]] <- iresult[['keg']];
	go.all[[i]] <- iresult[['go']];
};
save(keg.all, go.all, file = generate.filename('enrich_modules', name, 'rda'));
###
