### plot PAM50 subtypes score across clusters
### centroids file download from pam50 website
pam50=read.table( "pam50_centroids.txt", head=TRUE,  stringsAsFactors=FALSE )
# pam50[pam50<0]=0
### which gene names are not in single cell data
rownames(pam50)[!rownames(pam50) %in% rownames(epiTumor@data)]
rownames(pam50)[rownames(pam50)=="CDCA1"]="NUF2"
rownames(pam50)[rownames(pam50)=="KNTC2"]="NDC80"
rownames(pam50)[rownames(pam50)=="ORC6L"]="ORC6"
# epiTumor=ScaleData(object = epiTumor)
exp.m.pam50<-as.matrix(epiTumor@data[ rownames(pam50),])
table(rownames(exp.m.pam50)==rownames(pam50))
# TRUE
#   50
# calculate subtypes score for each cell
# epiTumor@meta.data$lum.a.score=apply(exp.m.pam50, 2, function(x) cor(x, pam50$"LumA"))
epiTumor@meta.data$lum.a.score=apply(exp.m.pam50, 2, function(x) cor(x, pam50$"LumA", method = "spearman"))
epiTumor@meta.data$lum.a.score[is.na(epiTumor@meta.data$lum.a.score)]=0

epiTumor@meta.data$lum.b.score=apply(exp.m.pam50, 2, function(x) cor(x, pam50$"LumB", method = "spearman"))
epiTumor@meta.data$lum.b.score[is.na(epiTumor@meta.data$lum.b.score)]=0

epiTumor@meta.data$basal.score=apply(exp.m.pam50, 2, function(x) cor(x, pam50$"Basal", method = "spearman"))
epiTumor@meta.data$basal.score[is.na(epiTumor@meta.data$basal.score)]=0
