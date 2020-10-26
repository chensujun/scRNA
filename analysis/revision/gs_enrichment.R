library(BoutrosLab.plotting.general);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/GS');
source('~/chensj/scRNA/script/myfunctions/test_kegg.R');
source('~/chensj/scRNA/script/myfunctions/plot_enrich_manual.R');
annot.keg <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_kegg.rds');
annot.go <- readRDS('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/crpc_norm/raw_data/2019-07-18_database_go.rds');
deg <- readRDS('celltype.deg.detail.rds');
myekeg <- myego <- list()
for(i in names(deg)){
	mygene <- rownames(deg[[i]])
	ego <- test_enrich(annot.keg, mygene);
	ego.go <- test_enrich(annot.go, mygene);
	assign(paste0('ego.', i), ego);
	assign(paste0('ego.go.', i), ego.go);
	myekeg[[i]] <- ego;
	myego[[i]] <- ego.go
}
save(myego, myekeg, file = generate.filename('enrich', 'GS_type', 'rda'));
ego.go <- ego.go[ego.go$FDR<0.05, ];
ego.go <- ego.go[order(ego.go$FDR), ];
plot_enrich(myego, 'GS_diff_go', mybg = rownames(res), nterms = 5, width = 20, height = 20,spot.size.function = function(x) {x/2}, key.sizes = c(2, 4));
plot_enrich(myekeg, 'GS_diff_keg', mybg = rownames(res), nterms = 5, width = 15, height = 20,spot.size.function = function(x) {x/2}, key.sizes = c(2, 4));
to.plot.ego <- to.plot.keg <- list();
#for(i in c('DC', 'Fib S2', 'Mast', 'TAM', 'Luminal')){
for(i in c('DC', 'Fib S1', 'Luminal', 'Mast', 'TAM')){
	to.plot.ego[[i]] <- myego[[i]];
	to.plot.keg[[i]] <- myekeg[[i]];
};
plot_enrich(to.plot.ego, 'GS_diff_go5', mybg = rownames(res), nterms = 5, width = 10, height = 10,spot.size.function = function(x) {x/2}, key.sizes = c(2, 4));
plot_enrich(to.plot.keg, 'GS_diff_keg5', mybg = rownames(res), nterms = 5, width = 7.6, height = 8.2,spot.size.function = function(x) {x/2}, key.sizes = c(2, 4));
