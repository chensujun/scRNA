####
#### the following script is adapted from Gervaise et al
####https://git.biohpc.swmed.edu/StrandLab/sc-TissueMapper_Pr/blob/master/r.scripts/sc-TissueMapper.R

run_qusage <- function(seurat.all, nm = 'clusters', gs, my.seed=100){
	# Run QuSAGE for each of the clusters for given geneset(s)
	#Inputs:
	#seurat.all = a seurat object from scran
	#nm = name of the run
	#gs = gene set to test, a list

	#Outputs:
	#results.cor = correlation table
	#results.clust.id = correlation results

	ds <- min(table(seurat.all@meta.data$cluster));
	number.clusters <- length(unique(seurat.all@meta.data$cluster))
	cell.sample <- NULL
	for (i in 1:number.clusters){
	cell <- rownames(seurat.all@meta.data[seurat.all@meta.data$cluster==i, ])
	if (length(cell)>ds & ds!=0){
		set.seed(my.seed)
		#set.seed(length(cell))
		rnd <- sample(1:length(cell),ds)
	 	cell <- cell[rnd]
	}
	cell.sample <- c(cell.sample,cell)
	}

	mydata <- seurat.all@data[, colnames(seurat.all@data)%in%cell.sample];
	labels <- paste0('Cluster_', seurat.all@meta.data[colnames(mydata), ]$cluster);
	  
	clust <- list()
	clust.comp <- list()
	  for (i in 1:number.clusters){
	    t <- labels
	    t[!(t %in% paste0("Cluster_",i))] <- "REST"
	    clust[i] <- list(i=t)
	    rm(t)
	    clust.comp[i] <- paste0("Cluster_",i,"-REST")
	  }
	  

	for (i in 1:number.clusters){
	    assign(paste0("results.",i),qusage(mydata,unlist(clust[i]),unlist(clust.comp[i]),gs))
	    
	 }

	results.cor <- NULL
	results.cor <- qsTable(results.1, number = length(gs));
	results.cor$Cluster <- 1
	for (i in 2:number.clusters){
		qs <- qsTable(get(paste0("results.",i)), number = length(gs))
		qs$Cluster <- i
		results.cor <- rbind(results.cor,qs)
	  };

	results.cor <- results.cor[,-3]
	rownames(results.cor) <- NULL
	results.clust.id <- NULL
	if (max(results.cor[results.cor[,4]==1,][,2])>=0){
	results.clust.id <- results.cor[results.cor[,4]==1,][which.max(results.cor[results.cor[,4]==1,][,2]),]
	} else {
	results.clust.id$pathway.name <- "Unknown"
	results.clust.id$log.fold.change <- 0
	results.clust.id$FDR <- 0
	results.clust.id$Cluster <- 1
	results.clust.id <- as.data.frame(results.clust.id)
	}
	for (i in 2:number.clusters){
	if (max(results.cor[results.cor[,4]==i,][,2])>=0){
	  results.clust.id <- rbind(results.clust.id,results.cor[results.cor[,4]==i,][which.max(results.cor[results.cor[,4]==i,][,2]),])
	} else {
	  results.clust.id <- rbind(results.clust.id,data.frame(pathway.name="Unknown",log.fold.change=0,FDR=0,Cluster=i))
	}
	}
	rownames(results.clust.id) <- NULL
  	save(list=c('gs', 'number.clusters', 'results.clust.id', 'results.cor', paste0('results.', seq(number.clusters))), file = paste0(Sys.Date(), '_', nm, '_qusage_results_', length(gs),'_',  my.seed, '.rda'));
	#Determine axes for correlation plots
	max.x.rg <- 0
	min.x.rg <- 0
	max.y.rg <- 0
	for (i in 1:number.clusters){
	qs <- get(paste0("results.",i))
	if (max(qs$path.mean)>max.x.rg){
	  max.x.rg <- max(qs$path.mean)
	}
	if (min(qs$path.mean)<min.x.rg){
	  min.x.rg <- min(qs$path.mean)
	}
	if (max(qs$path.PDF)>max.y.rg){
	  max.y.rg <- max(qs$path.PDF)
	}}
	#Plot correlation plots by geneset
	for (i in 1:length(gs)){
		dir.create(file.path('./', 'analysis'));
		dir.create(file.path('./analysis', 'cor'));
	pdf(paste0("./analysis/cor/QuSAGE_",nm,".",names(gs)[i],".pdf"), width = number.clusters/2)
	for (j in 1:number.clusters){
	  qs <- get(paste0("results.",j))
	  if (j==1){
	    plotDensityCurves(qs,path.index=i,col=viridis(number.clusters)[j],main=names(gs)[i],xlim=c(min.x.rg-0.05,max.x.rg+0.05),ylim=c(0,50*ceiling(max.y.rg/50)),xlab="Gene Set Activation",lwd=5,cex.main=2.5,cex.axis=1.5,cex.lab=2)
	  } else {
	    plotDensityCurves(qs,path.index=i,add=TRUE,col=viridis(number.clusters)[j],lwd=5)
	}}
	leg <- paste0("Cluster ",1:number.clusters)
	legend("topright",legend=leg,lty=1,col=viridis(number.clusters),lwd=5,cex=2, bty="n",pt.cex=2, ncol = ceiling(number.clusters/10))
	box(lwd=5)
	dev.off()
	}
	results <- list(
		results.cor = results.cor, 
		results.clust.id = results.clust.id
		);
	return(results);
};

