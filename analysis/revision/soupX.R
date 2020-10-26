library(BoutrosLab.plotting.general);
library(SoupX);
library(ggplot2);
library(Matrix);
setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision');
name <- 'all';
dataDirs = c("sample13/");
sc = load10X(dataDirs, keepDroplets = TRUE);
sc = estimateSoup(sc);
all_dr <- readRDS('2020-02-11_embedding_cluster_all.rds');
all_dr$KLK3 = sc$toc["KLK3", rownames(all_dr)]
pdf(generate.filename('exp_KLK3', name, 'pdf'));
gg <- ggplot(all_dr, aes(tSNE_1, tSNE_2)) + geom_point(aes(colour = KLK3 > 0));
plot(gg);
dev.off();
###
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
    "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC");

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes));

pdf(generate.filename('cells', 'estimateSoup', 'pdf'));
plotMarkerMap(sc, geneSet = igGenes, DR = all_dr, useToEst = useToEst)
dev.off();

all_dr$useToEst <- useToEst[rownames(all_dr), 1];

pdf(generate.filename('cells', 'estimateSoup', 'pdf'));
gg <- ggplot(all_dr, aes(tSNE_1, tSNE_2)) + geom_point(aes(colour = useToEst));
plot(gg);
dev.off();

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), 
    clusters = setNames(all_dr$cluster, rownames(all_dr)))

sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst);
out <- adjustCounts(sc);
saveRDS(out, file = generate.filename('soupx_corrected', 'ig', 'rds'))
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10);
all_dr$klk3_c <- out['KLK3', rownames(all_dr)];

pdf(generate.filename('corrected_exp_KLK3', name, 'pdf'));
gg <- ggplot(all_dr, aes(tSNE_1, tSNE_2)) + geom_point(aes(colour = klk3_c > 0));
plot(gg);
dev.off();
saveRDS(all_dr, file = generate.filename('embedding', 'corrected_klk3', 'rds'))
###
t_dr <- readRDS('2020-02-11_embedding_cluster_t.rds');
name <- 'tcell';
t_dr$KLK3 = sc$toc["KLK3", rownames(t_dr)]
pdf(generate.filename('exp_KLK3', name, 'pdf'));
gg <- ggplot(t_dr, aes(tSNE_1, tSNE_2)) + geom_point(aes(colour = KLK3 > 0));
plot(gg);
dev.off();
###
pdf(generate.filename('mm_KLK3', name, 'pdf'));
gg <- plotMarkerMap(sc, 'KLK3', t_dr);
plot(gg);
dev.off();
endo_dr <- readRDS('2020-02-11_embedding_cluster_endo.rds')
endo_dr$THY1 = sc$toc["THY1", rownames(endo_dr)]
pdf(generate.filename('mm_THY1', 'endo', 'pdf'));
gg <- plotMarkerMap(sc, 'THY1', endo_dr);
plot(gg);
dev.off();

pdf(generate.filename('mm_KLK3', 'all', 'pdf'));
gg <- plotMarkerMap(sc, 'KLK3', all_dr);
plot(gg);
dev.off();

all_dr$THY1 = sc$toc["THY1", rownames(all_dr)]
pdf(generate.filename('mm_THY1', 'all', 'pdf'));
gg <- plotMarkerMap(sc, 'THY1', all_dr);
plot(gg);
dev.off();

pdf(generate.filename('plotgene', paste0(name, '_KLK3'), 'pdf'), width = 6, height = 6);
FeaturePlot(seurat.all, reduction.use = 'tsne', features.plot = 'KLK3', pt.size = 1, nCol = 1, min.cutoff = 0.1,
                cols.use = c('grey', 'red'), vector.friendly = TRUE, no.axes = TRUE);
dev.off();
