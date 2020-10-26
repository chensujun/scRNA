library(BoutrosLab.plotting.general);
source('~/svn/singleCell/myfunctions/test_kegg.R');
source('~/svn/singleCell/myfunctions/create.heatmap_modi_key.R');
name <- 'tcell';
mscore <- read.table('GraphClust/GraphClust.ModuleExp.txt', row.names = 1, sep = '\t', header = TRUE);
mgene <- read.table('GraphClust/GraphClust.module.txt', header = TRUE);
rownames(mscore) <- gsub(' ', '_', rownames(mscore));
annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
####
mscore$maxid <- apply(mscore, 1, function(x) grep(max(x), x));

module.cd8e <- mscore[rowSums(mscore[, c(3,4,6)]>0)==3&mscore$maxid%in%c(3,4,6), ];
module.cd8e <- module.cd8e[apply(module.cd8e[, -8], 1, function(x) min(x[c(3,4,6)])>max(x[-c(3,4,6)])), ];
gene.cd8e <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.cd8e)), ]$id);
keg.cd8e <- test_enrich(annot.keg, gene.cd8e, mybg = NULL);
go.cd8e <- test_enrich(annot.go, gene.cd8e, mybg = NULL);
###
module.cd8n <- mscore[rowSums(mscore[, c(2, 5)]>0)==2&mscore$maxid%in%c(2, 5), ];
module.cd8n <- module.cd8n[apply(module.cd8n[, -8], 1, function(x) min(x[c(2,5)])>max(x[-c(2,5)])), ];
gene.cd8n <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.cd8n)), ]$id);
keg.cd8n <- test_enrich(annot.keg, gene.cd8n, mybg = NULL);
go.cd8n <- test_enrich(annot.go, gene.cd8n, mybg = NULL);
###
module.treg <- mscore[mscore$Cluster6>0.1&mscore$maxid==7, ];
module.treg <- module.treg[rowSums(module.treg[, -8]>0)<3&rowSums(module.treg[, -8]<0.1)>5, ];
module.treg <- module.treg[order(-module.treg$Cluster6), ]
gene.treg <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.treg)), ]$id);
keg.treg <- test_enrich(annot.keg, gene.treg, mybg = NULL);
go.treg <- test_enrich(annot.go, gene.treg, mybg = NULL);

module.conv <- mscore[mscore$Cluster0>0&mscore$maxid==1, ];
module.conv <- module.conv[module.conv$Cluster0>0.1&rowMeans(module.conv[, -8])<0, ];
gene.conv <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.conv)), ]$id);
keg.conv <- test_enrich(annot.keg, gene.conv, mybg = NULL);
go.conv <- test_enrich(annot.go, gene.conv, mybg = NULL);

module.klk <- mscore[mscore$Cluster5>0&mscore$maxid==6, ];
module.klk <- module.klk[module.klk$Cluster5>0.1&rowMeans(module.klk[, -8])<0, ];
gene.klk <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.klk)), ]$id);
keg.klk <- test_enrich(annot.keg, gene.klk, mybg = NULL);
go.klk <- test_enrich(annot.go, gene.klk, mybg = NULL);
save(keg.treg, go.treg, keg.cd8n, keg.cd8e, keg.conv, keg.klk, go.cd8n, go.cd8e, go.conv, go.klk, 
	file = generate.filename('enrich_subtype', name, 'rda'));
save(module.treg, module.cd8n, module.cd8e, module.conv, module.klk, file = generate.filename('modules_unique', name, 'rda'))
###
module.klkm <- mscore[mscore$Cluster3>0&mscore$maxid==4, ];
module.klkm <- module.klkm[module.klkm$Cluster3>0.1&rowMeans(module.klkm[, -8])<0, ];
module.klkm <- module.klkm[!rownames(module.klkm)%in%rownames(module.cd8e), ];
gene.klkm <- as.vector(mgene[mgene$module%in%gsub('Module_', '', rownames(module.klkm)), ]$id);
keg.klkm <- test_enrich(annot.keg, gene.klkm, mybg = NULL);
go.klkm <- test_enrich(annot.go, gene.klkm, mybg = NULL);
save(keg.treg, go.treg, keg.cd8n, keg.cd8e, keg.conv, keg.klk, go.cd8n, go.cd8e, go.conv, go.klk, keg.klkm, go.klkm,
	file = generate.filename('enrich_subtype', name, 'rda'));
save(module.treg, module.cd8n, module.cd8e, module.conv, module.klk, module.klkm, file = generate.filename('modules_unique', name, 'rda'))

###
to.plot <- mscore[c(rownames(module.treg), rownames(module.conv), rownames(module.cd8n),rownames(module.cd8e), rownames(module.klk), rownames(module.klkm)), -8];
to.plot <- data.frame(t(scale(t(to.plot))));
to.plot <- to.plot[, c(7, 1, 2, 5, 3, 4, 6)]
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
### test keg and go for each module
test_module <- function(mymodule, mgene, annot.keg, annot.go){
	gene.mymodule <- as.vector(mgene[mgene$module%in%gsub('Module_', '', mymodule), ]$id);
	keg.mymodule <- test_enrich(annot.keg, gene.mymodule, mybg = NULL);
	go.mymodule <- test_enrich(annot.go, gene.mymodule, mybg = NULL);
	myresult <- list(keg = keg.mymodule, go = go.mymodule);
	return(myresult);
};

keg.all <- go.all <- list();
for(i in rownames(to.plot)){
	#for(i in rownames(module.klkm)){
	print(i);
	iresult <- test_module(i, mgene, annot.keg, annot.go);
	keg.all[[i]] <- iresult[['keg']];
	go.all[[i]] <- iresult[['go']];
};
save(keg.all, go.all, file = generate.filename('enrich_modules', name, 'rda'));
###
filter_term <- function(mytest){
	itest <- mytest[mytest$FDR<0.25&mytest$nPath>10&mytest$P.val<0.01, ];
	itest <- itest[order(-itest$Enrichment), ];
	return(itest);
};

for(i in names(go.all)){
	#go.all[[i]] <- filter_term(go.all[[i]]);
	print(i)
	print(go.all[[i]][1:10, ])
};

for(i in names(go.all)){
	igo <- go.all[[i]]
	ifelse(min(igo$FDR<0.05), print(i), print(''))
};

### test module 61 gene
ref.gene <- read.xls('raw_data/12943_2014_285_MOESM3_ESM.xls');
ref.gene$log2FoldChange <- as.numeric(as.vector(ref.gene$log2FoldChange));
gene.reg <- read.table('~/tmp/plot/TRANSFAC_and_JASPAR_PWMs_table.txt', sep = '\t', quote = '', header = TRUE);


 "Module_9"   "Module_75"  "Module_94"  "Module_99"  "Module_102"
aa <- go.all[['Module_102']]
aa <- aa[aa$nPath>10&aa$FDR<0.25&aa$P.val<0.01, ]
