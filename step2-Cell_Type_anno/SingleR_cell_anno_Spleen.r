library(SingleR)
library(celldex)
library(Seurat)
index="P_SingleR_Tissue"
RDS=readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/kangjingmin/project/01.Mt_sc/10.Article/04.Distribution/01.Old/03.Cluster_SignalC/01.P_TSNE_50/P_combined.harmony_GeneActivityLogNormalize.rds")
Ref=MouseRNAseqData()
DefaultAssay(RDS) <- 'RNA'
#Ref=ImmGenData() 
cluster=RDS@meta.data$peaks_snn_res.0.4
singleR=GetAssayData(RDS,slot="data")
pred=SingleR(test=singleR,ref=Ref,labels=Ref$label.main,method="cluster",clusters=cluster)
cellType=data.frame(ClusterID=levels(RDS@meta.data$peaks_snn_res.0.4),celltype=pred$labels)
print (unique(cellType))
RDS@meta.data$CellType=cellType[match(cluster,cellType$clusterID),'celltype']
saveRDS(RDS, file=paste0(index,"_combined.harmony_FindClusters.rds") )


index="P_SingleR_Imm"
RDS=readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/kangjingmin/project/01.Mt_sc/10.Article/04.Distribution/01.Old/03.Cluster_SignalC/01.P_TSNE_50/P_combined.harmony_GeneActivityLogNormalize.rds")
DefaultAssay(RDS) <- 'RNA'
Ref=ImmGenData() 
cluster=RDS@meta.data$peaks_snn_res.0.4
singleR=GetAssayData(RDS,slot="data")
pred=SingleR(test=singleR,ref=Ref,labels=Ref$label.main,method="cluster",clusters=cluster)
cellType=data.frame(ClusterID=levels(RDS@meta.data$peaks_snn_res.0.4),celltype=pred$labels)
print (unique(cellType))
RDS@meta.data$CellType=cellType[match(cluster,cellType$clusterID),'celltype']
saveRDS(RDS, file=paste0(index,"_combined.harmony_FindClusters.rds") )


