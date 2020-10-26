get_differential <- function(epiTumor,clusters.1,clusters.2,patients,pct = 0.1,which.label){
  ####
  temp.info=c()
  patients.iden=epiTumor@meta.data$fig.patient %in% patients
  temp.cells=epiTumor@cell.names[patients.iden]
  temp.epiTumor=SubsetData(object = epiTumor, cells.use = temp.cells)
  temp.ident=rep(0,length(temp.epiTumor@cell.names))
  clusters.1.iden=temp.epiTumor@meta.data[,which.label] %in% clusters.1
  clusters.2.iden=temp.epiTumor@meta.data[,which.label] %in% clusters.2
  
  temp.ident[clusters.1.iden]=1
  temp.ident[clusters.2.iden]=2
  temp.epiTumor@meta.data$"temp.ident"=temp.ident
  temp.epiTumor=SetAllIdent(object = temp.epiTumor, id = "temp.ident")
  temp.markers.all=FindMarkers(object = temp.epiTumor, ident.1 = 1,ident.2=2, only.pos=FALSE);
  ####get fc of two groups
  temp.markers.all$mean1 <- rowMeans(as.matrix(temp.epiTumor@data[rownames(temp.markers.all), idata@ident==1]));
  temp.markers.all$mean2 <- rowMeans(as.matrix(temp.epiTumor@data[rownames(temp.markers.all), idata@ident==2]));
  temp.markers.all$lfc <- temp.markers.all$mean1 - temp.markers.all$mean1;
  temp.markers.all <- temp.markers.all[order(temp.markers.all$lfc), ];
  ####
  if(nrow(temp.markers.all<=50)){
    plot.genes <- rownames(temp.markers.all)
  }else{
    plot.genes <-  rownames(temp.markers.all)[c(1:25, (nrow(temp.markers.all)-24):nrow(temp.markers.all))]
  };

  filename=paste(paste("cluster.",paste(as.vector(clusters.1),collapse = '.'),sep=""),paste("clusters.",paste(as.vector(clusters.2),collapse = '.'),sep=""),sep="_v.s._")
  filename=paste(filename,"_patient.",paste(patients,collapse="."),sep="")
  print("all.markers::")
  str(all.markers)
  print("WhichCells(temp.epiTumor,c(clusters.1,clusters.2))::")
  print(table(temp.epiTumor@ident %in% c(1,2)))

  temp.epiTumor=SetAllIdent(object = temp.epiTumor, id = which.label)
  temp.order=c(sort(clusters.1),sort(clusters.2))
  pdf(paste("heatmap.",filename,".pdf",sep=""))
  print(DoHeatmap(temp.epiTumor, genes.use = all.markers,
    slim.col.label = TRUE,remove.key = TRUE,cells.use =WhichCells(temp.epiTumor,c(clusters.1,clusters.2)), 
    group.order=temp.order, group.cex=5))
  dev.off();

  return(temp.info)
}