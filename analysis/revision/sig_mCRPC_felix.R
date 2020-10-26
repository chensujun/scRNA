setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/sig');
d1 <- read.table('/cluster/projects/hansengroup/Public_Dataset_Hub/changhai/rnaseq/filter_genes.FPKM.xls', header = TRUE, row.names = 1);
d1 <- d1[, grep('^T', colnames(d1))];
c1 <- read.csv('/cluster/projects/hansengroup/Public_Dataset_Hub/changhai/clinical/clinical.csv')

zhengli jieguo 
IF for co-culture 
signature, weight and etc. changhai cohort significance 
###
#use qusage to calculate score enrichment
library(qusage);
run_qusage <- function(mydata, nm = 'samples', gs, my.seed=100){
		labels <- colnames(mydata);
		clust <- list()
		clust.comp <- list()
		  for (i in labels){
		    t <- labels
		    t[!(t %in% i)] <- "REST"
		    clust[i] <- list(i=t)
		    rm(t)
		    clust.comp[i] <- paste0(i,"-REST")
		  }
		  

		for (i in labels){
		    assign(paste0("results.",i),qusage(mydata,unlist(clust[i]),unlist(clust.comp[i]),gs))
		    
		 }

}
setwd('/Users/sujunchen/Document/labs/2020/scRNA/signature');
mydata <- readRDS('data/2019-01-05_tcga_fpkm.rds');
colnames(mydata) <- gsub('\\.', '_', colnames(mydata));
mydata <- as.matrix(mydata)
gs <- list(basal = read.table('/Users/sujunchen/Document/labs/2020/scRNA/signature/data/basal.in.all.cell.markers.posi.genes.txt')$V1,
	cycle = read.table('/Users/sujunchen/Document/labs/2020/scRNA/signature/data/cycle.in.all.cell.markers.posi.genes.txt')$V1);
gs[[1]] <- intersect(rownames(mydata), gs[[1]]);
gs[[2]] <- intersect(rownames(mydata), gs[[2]]);
mydata <- log2(mydata+1)