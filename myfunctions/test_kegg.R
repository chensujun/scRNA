library(BoutrosLab.plotting.general);
#rownames(seurat.cr@data);
#grep -w gene /.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/Homo_sapiens.GRCh38.91.gtf |cut -f 9 |sed 's/.*gene_name//g'|sed 's/;.*//g'| sed 's/"//g' > /.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/GRCh38_gene_name.txt
### format go and kegg annotation data
#gene.go <- readRDS('~/circRNA/star-circ/circRNA_landscape/rnaseq_landscape/scRNA/2019-01-05_gp_hsapiens.GO.Name.rds');
#ref.go <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/data/ref/gp_ref/hsapiens.GO.NAME.ann', sep = '\t', quote = '');
#gene.go$name <- ref.go[match(gene.go$ont, ref.go$V1), ]$V2;
#saveRDS(gene.go, file = generate.filename('database', 'go', 'rds'));
##annot.keg <- readRDS('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-05_kegg_human_20190613.rds');
#annot.keg$name <- annot.keg$ont;
#annot.keg$ont <- ref.keg[match(annot.keg$ont, ref.keg$PathwayTerm), ]$PathwayID;
#saveRDS(annot.keg, file = generate.filename('database', 'kegg', 'rds'));
###/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/scran/crpc_norm/raw_data/2019-07-18_database_go.rds
test_enrich <- function(annot.keg, jup, mybg = NULL, underep = FALSE){
	if(is.null(mybg)){
		mybg <- read.table('/.mounts/labs/cpcgene/private/projects/rnaseq_landscape/scRNA/data/GRCh38_gene_name.txt');
		mybg <- intersect(mybg$V1, annot.keg$gene);
	};
	mydata <- data.frame(listgene = ifelse(mybg%in%jup, 1, 0));
	rownames(mydata) <- mybg;
	mydata <- mydata[rownames(mydata)%in%annot.keg$gene, drop = FALSE, ];
	if(underep){
		mypaths <- unique(annot.keg[annot.keg$gene%in%mybg, ]$ont);
		}else{
			mypaths <- unique(annot.keg[annot.keg$gene%in%jup, ]$ont)
		}
	mytest <- data.frame(matrix(NA, ncol = 10, nrow = length(mypaths)));
	colnames(mytest) <- c('TermName', 'nBG', 'nPath', 'nDiff', 'nIntersect', 'P.val', 'Enrichment', 'OR', 'FDR','Intersect');
	rownames(mytest) <- mypaths;
	j <- 0
	for(i in mypaths){
		j <- j + 1;
		if(j%%1000==0){
			print(paste0('test ', j, ' out of ', length(mypaths), ' terms...'))
		};
		mydata$pathgene <- ifelse(rownames(mydata)%in%annot.keg[annot.keg$ont==i, ]$gene, 1, 0);
		itest <- fisher.test(mydata$pathgene, mydata$listgene);
		nf <- nrow(mydata[mydata$pathgene==1&mydata$listgene==1, ]);
		n <- nrow(mydata[mydata$listgene==1, ]);
		Nf <- nrow(mydata[mydata$pathgene==1, ]);
		N <- nrow(mydata);
		enrich <- (nf/n)/(Nf/N);
		mytest[i, 1] <- as.vector(unique(annot.keg[annot.keg$ont==i, ]$name));
		mytest[i, 2:8] <- c(N, Nf, n, nf, itest$p.value, enrich, itest$estimate);
		mytest[i, 10] <- paste(rownames(mydata[mydata$pathgene==1&mydata$listgene==1, ]), collapse = ',');
	};
	mytest$FDR <- p.adjust(mytest$P.val, method = 'BH');
	mytest <- mytest[order(mytest$FDR), ];
	return(mytest);
};

run_enrich <- function(annot.keg, annot.go, jup, mybg.keg = NULL, mybg.go = NULL){
	ego.keg <- test_enrich(annot.keg, jup, mybg.keg);
	ego.go <- test_enrich(annot.go, jup, mybg.go);
};