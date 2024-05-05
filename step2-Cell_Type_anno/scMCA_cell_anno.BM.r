library(scMCA)
library(Seurat)
scATAC=readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/kangjingmin/project/01.Mt_sc/10.Article/04.Distribution/01.Old/03.Cluster_SignalC/02.GS_TSNE_50/GS_combined.harmony_GeneActivityLogNormalize.rds")

DefaultAssay(scATAC) <- 'RNA'

scATAC_result <- scMCA(scdata = scATAC@assays$RNA@data, numbers_plot = 3)
out=as.data.frame(unlist(scATAC_result$scMCA))
out$`unlist(scATAC_result$scMCA)`=as.character(out$`unlist(scATAC_result$scMCA)`)
scATAC@meta.data$cell_type=out[match(rownames(scATAC@meta.data),rownames(out)),1]
out_meta=scATAC@meta.data
table_list=list()

index="scMCA_GS"
scATAC$seurat_clusters=scATAC$peaks_snn_res.0.4

for(i in 0:(length(unique(scATAC$seurat_clusters))-1)){
    a=i+1
    sub=subset(out_meta,out_meta$seurat_clusters==i)
    tab=as.data.frame(table(sub$cell_type))
    tab_order=tab[order(tab[,2],decreasing = T),]
    tab_order$Cluster=i
    aa=cbind(i,tab_order[1,],tab_order[2,])
    print (aa)
    table_list[[a]]=tab_order[1,]
  }

  sum=do.call(rbind,table_list)
  use=out_meta[,c(22,23)]
  use$seurat_clusters=as.character(use$seurat_clusters)
  use$ID=rownames(out_meta)
  colnames(sum)=c("predicated.cell.type","Freq","seurat_clusters")
  res=merge(use,sum,by="seurat_clusters")
  scATAC@meta.data$predicated.cell.type=res[match(rownames(out_meta),res$ID),4]
  scATAC@meta.data$predicated.cell.type=paste0(scATAC@meta.data$seurat_clusters,":",scATAC@meta.data$predicated.cell.type)

  library(RColorBrewer)
  library(ggplot2)
  getPalette = colorRampPalette(brewer.pal(12, "Set3"))
  a=getPalette(length(unique(scATAC@meta.data$predicated.cell.type)))
  c3=DimPlot(object = scATAC, label = TRUE, group.by = "predicated.cell.type") +scale_color_manual(values = a)
  png(paste0(paste0(index,"_combined.harmony_FindClusters.png")),width = 990,height = 479)
  print(c3)
  dev.off()
    table_cell.type=as.data.frame(table(as.character(scATAC@meta.data$predicated.cell.type)))
  colnames(table_cell.type)=c("predicated.cell.type","number")
  table_cell.type$ratio=table_cell.type$number/colSums(table_cell.type[2])
  order_table_cell.type=table_cell.type[order(table_cell.type$number,decreasing = T),]  
  write.table(order_table_cell.type,paste0(index,"_Table2.celltypecount.csv"),sep = ",",quote = FALSE,row.names = FALSE)

  saveRDS(scATAC, file=paste0(index,"_combined.harmony_FindClusters.rds") )



