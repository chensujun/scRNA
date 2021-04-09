require("scran")
require("Seurat")

csce <- readRDS('2019-03-21_all_20k_annotated.rds');
# cell type annotation
mycolData <- readRDS('2019-03-27_annotate_type_colData_manual.rds');
mycolData.2 <- readRDS('2019-07-11_seurat_coldata_all.rds');
csce@colData <- mycolData;
table(rownames(csce@colData)==rownames(mycolData.2))
csce@colData$type <- mycolData.2$type
pros13<-as.seurat(csce)
pros13@data=log(2^pros13@data)
# saveRDS(pros13, "pros13.identifyDEGsComparing2allKindsOfCells.20190822.rds")
pros13=readRDS("pros13.identifyDEGsComparing2allKindsOfCells.20190822.rds")
finer=readRDS("2020-06-27_celltype_all_pros13_refinedCellTypeAnno_fromHunhun.rds")

aEC.cells.list=rownames(finer)[finer$type=="aEC"]
basal.cells.list=readRDS("basal.cells.list.txt")
cycle.cells.list=readRDS("cycle.cells.list.txt")
pros13@meta.data$basal=rep("other", length(pros13@cell.names))
pros13@meta.data$basal[pros13@cell.names %in% basal.cells.list]="basal"
table(pros13@meta.data$basal)
pros13@meta.data$cycle=rep("other", length(pros13@cell.names))
pros13@meta.data$cycle[pros13@cell.names %in% cycle.cells.list]="cycle"
table(pros13@meta.data$cycle)
pros13@meta.data$aEC=rep("other", length(pros13@cell.names))
pros13@meta.data$aEC[pros13@cell.names %in% aEC.cells.list]="aEC"
table(pros13@meta.data$aEC)

#### get basal/intermediate cells's markers
pros13=SetAllIdent(object = pros13, id = "basal")
markers.pros13=FindMarkers(pros13, ident.1 = "basal",ident.2="other", only.pos=FALSE, logfc.threshold=log(1.5), min.pct = 0.1)
# saveRDS(markers.pros13, "basal.in.all.cell.markers")
# markers.pros13=readRDS("basal.in.all.cell.markers")
markers.pros13=markers.pros13[markers.pros13$p_val_adj<0.05, ]
markers.pros13=markers.pros13[order(markers.pros13$avg_logFC, decreasing=TRUE), ]
table(markers.pros13$avg_logFC>0)
basal.posi.genes=rownames(markers.pros13)[markers.pros13$avg_logFC>0]
str(basal.posi.genes)

write.table(basal.posi.genes, file="basal.in.all.cell.markers.posi.genes.txt", col.name=F, row.name=F, quote=F, sep='\t')

#### get cycle-enhanced cells's markers
pros13=SetAllIdent(object = pros13, id = "cycle")
markers.pros13=FindMarkers(pros13, ident.1 = "cycle",ident.2="other", only.pos=FALSE, logfc.threshold=log(1.5), min.pct = 0.1)
saveRDS(markers.pros13, "cycle.in.all.cell.markers")
markers.pros13=markers.pros13[markers.pros13$p_val_adj<0.05, ]
markers.pros13=markers.pros13[order(markers.pros13$avg_logFC, decreasing=TRUE), ]
table(markers.pros13$avg_logFC>0)
cycle.posi.genes=rownames(markers.pros13)[markers.pros13$avg_logFC>0]
write.table(cycle.posi.genes, file="cycle.in.all.cell.markers.posi.genes.txt", col.name=F, row.name=F, quote=F, sep='\t')

#### get cycle-enhanced cells's markers
pros13=SetAllIdent(object = pros13, id = "aEC")
markers.pros13=FindMarkers(pros13, ident.1 = "aEC",ident.2="other", only.pos=FALSE, logfc.threshold=log(1.5), min.pct = 0.1)
saveRDS(markers.pros13, "aEC.in.all.cell.markers")
markers.pros13=markers.pros13[markers.pros13$p_val_adj<0.05, ]
markers.pros13=markers.pros13[order(markers.pros13$avg_logFC, decreasing=TRUE), ]
table(markers.pros13$avg_logFC>0)
aEC.posi.genes=rownames(markers.pros13)[markers.pros13$avg_logFC>0]
write.table(aEC.posi.genes, file="aEC.in.all.cell.markers.posi.genes.txt", col.name=F, row.name=F, quote=F, sep='\t')



# Rscript km_sig_tcga.R sb mylist.txt mean ; Rscript km_sig_EBioMedicine2015.R sb mylist.txt mean ; Rscript km_sig_MCTP.Nature.2012.R sb mylist.txt mean ; Rscript km_sig_MSKCC.Cancer.Cell.2010.R sb mylist.txt mean ; Rscript km_sig_cpc.R sb mylist.txt mean ; Rscript km_sig_SwedishWatchfulWaiting.GSE16560.R sb mylist.txt mean ; Rscript km_sig_Glinsky.R sb mylist.txt mean ; Rscript km_sig_EU.R sb mylist.txt mean ;


