#!/usr/bin/Rscript
parser = argparse::ArgumentParser(description="Script to clustering scATAC data")
parser = argparse::ArgumentParser(description=cat("Script to merge scATAC library, dimensionality reduction and clustering using Signac(reference pipeline https://satijalab.org/signac/articles/merging.html)\n\n# Authors:   Maven\n# Contact information:   mawen3@genomics.cn\n# Date: 2022-04-7\n# R package version: ggplot2_3.3.5        cowplot_1.1.1 data.table_1.14.0 future_1.23.0 GenomicRanges_1.44.0 GenomeInfoDb_1.28.4 IRanges_2.26.0       S4Vectors_0.30.2 BiocGenerics_0.38.0 SeuratObject_4.0.4 Seurat_4.0.5 Signac_1.4.0\n\n"))
parser$add_argument('-R','--report', help='input path of report file dir')
parser$add_argument('-s','--samples',default="Samples" ,help='sample name[default = \"%(default)s\"]')
parser$add_argument('-M','--marker', help='marker list')
parser$add_argument('-T','--tissue', help='tissue')
parser$add_argument('-D','--dims', help='dims using for dimensionality reduction')
parser$add_argument('-O','--out', help='out directory')

args = parser$parse_args()
library(dplyr)
library(data.table)
library(Seurat)
library(patchwork)
library(Signac)
library(GenomeInfoDb)
library(BSgenome)
#library(BSgenome.OSativa.NCBI.IRGSPv1.0)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(TxDb.Osativa.NCBI.IRGSPv1.0)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(stringr)
library(IRanges)
library(GenomicRanges)
library(harmony)

tissue <- args$tissue
Marker=args$marker
inputPath <-args$report
outdir <- args$out
sampleFile <- args$samples
dir.create(outdir)
########################################
setwd(outdir)
samples <- read.table(sampleFile,header=T)
print(samples)
snATAC=readRDS(paste0(tissue,"_combined.harmony_GeneActivityLogNormalize.rds"))

DefaultAssay(snATAC) <- 'RNA'

################################ plot marker genes ################################
inGenelist <- Marker
Genelist <- read.table(inGenelist,sep="\t",stringsAsFactors=F,header=F)
print(rownames(snATAC@assays$RNA@counts))
if(ncol(Genelist)==1){
        Genelist <- Genelist[which(Genelist$V1 %in% rownames(snATAC@assays$RNA@counts)),]
        genes <- as.character(Genelist)
        genes_title <- as.character(Genelist)
        }else{
        Genelist <- Genelist[which(Genelist$V1 %in% rownames(snATAC@assays$RNA@counts)),]  
        genes <- as.character(Genelist$Gene)
        genes_title <- as.character(Genelist$V2)
        }
# plot 20220524.combined_markers #
print (Genelist)
genes <- unique(Genelist)
# 气泡图
pdf(paste0(tissue,"_combined.harmony_20220527.combined_markers.DotPlot.pdf"),width=80,height=10)
p <- DotPlot(snATAC, features=genes, group.by='seurat_clusters')+RotatedAxis()
print(p)
dev.off()

# 小提琴图
pdf(paste0(tissue,"_combined.harmony_20220527.combined_markers.VlnPlot.pdf"),width=20,height=550)
p <- VlnPlot(snATAC, features=genes, pt.size=0)
print(p)
dev.off()

########### Find differentially accessible peaks between clusters ###########
DefaultAssay(snATAC) <- 'peaks'
da_peaks <- FindAllMarkers(snATAC,min.pct = 0.2,test.use = 'LR')
head(da_peaks)

write.table(da_peaks,paste0(tissue,".Table1.combined.harmony_by_dataset.all_da_peaks.xls"),sep = "\t",quote = FALSE,row.names = TRUE)

openRegion1 <- rownames(da_peaks)
closest_genes1 <- ClosestFeature(snATAC, regions = openRegion1)

da_peaks$query_region <- rownames(da_peaks)
closest_genes1.join <- full_join(da_peaks,closest_genes1,by="query_region")

write.table(closest_genes1.join,paste0(tissue,".Table2.combined.harmony_by_dataset.da_peaks.closest_genes.xls"),sep = "\t",quote = FALSE,row.names = FALSE)

#
sub_da_peak=subset(da_peaks,da_peaks$p_val_adj<0.01 & da_peaks$avg_log2FC > 1)

openRegion <- rownames(sub_da_peak)
closest_genes <- ClosestFeature(snATAC, regions = openRegion)

head(closest_genes)

write.table(closest_genes,paste0(tissue,".Table3.combined.harmony_by_dataset.top_da_peaks.closest_genes.csv"),sep = ",",quote = FALSE,row.names = FALSE)

saveRDS(snATAC, file =paste0(tissue,"_combined.harmony_by_dataset.closest_genes.rds") )

# 与所有marker比较，找到da_peak在marker TSS附近的
DefaultAssay(snATAC) <- 'peaks'
da_peaks <- read.table(paste0(tissue,".Table1.combined.harmony_by_dataset.all_da_peaks.xls"),sep = "\t")
closest_genes <- ClosestFeature(snATAC, regions = rownames(da_peaks))
da_peaks$query_region <- rownames(da_peaks)
closest_genes.join <- full_join(da_peaks,closest_genes,by="query_region")
Genelist$gene_id <- Genelist
da_peaks <- full_join(closest_genes.join,Genelist,by="gene_id")

write.table(da_peaks,paste0(tissue,".Table2.combined.harmony_by_dataset.da_peaks.20220527.combined_markers.xls"),sep = "\t",quote = FALSE,row.names = FALSE)

##################################################################
################## 2022.06.21 plot 感兴趣的基因 ##################
snATAC <- readRDS(paste0(tissue,"_combined.harmony_GeneActivityLogNormalize.rds"))
snATAC
# An object of class Seurat 
# 99043 features across 23338 samples within 2 assays 
# Active assay: peaks (63377 features, 63377 variable features)
 # 1 other assay present: RNA
 # 3 dimensional reductions calculated: lsi, tsne, harmony
DefaultAssay(snATAC) <- 'RNA'
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.55, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# # 20220621.interestin_genes
# Anther
# "Os01g0831000","Os04g0396500","Os04g0168100"

genes <-c("Clec4g","Fabp4","Pecam1","Egfl7",
"Retnlg","S100a8","S100a9",
"Cd79a","Cd79b",
"Clec4f","Vsig4",
"Alb","Apoe",
"Epcam","Tspan8","Spp1",
"Nkg7",
"Xcl1","Ncr1",
"Cd3g","Cd3e",
"Mki67","Top2a",
"Xcr1","Clec9a",
"Cd209a",
"Ccr9",
"Ccr2")
pdf(file=paste0(tissue,"_20220621.snATAC.interesting_genes.DotPlot.pdf"),width=12,height=12)
p <- DotPlot(snATAC, features=genes)+ coord_flip()
print(p)
dev.off()

#  The following requested variables were not found: Os01g0831000, Os04g0396500

# pdf(file=paste0(tissue,"_20220621.snATAC.interesting_genes.StackedVlnPlot.pdf"),width=8,height=10)
# p <- StackedVlnPlot(snATAC, features=c("Os01g0831000","Os04g0396500","Os04g0168100"), pt.size=0)
# print(p)
# dev.off()s

# 小提琴图-有散点
pdf(paste0(tissue,"_20220621.snATAC.interesting_genes.VlnPlot.pdf"),width=10,height=20)
p <- VlnPlot(snATAC, features=genes, pt.size=0.1)
print(p)
dev.off()
# 小提琴图-无散点
pdf(paste0(tissue,"_20220621.snATAC.interesting_genes.VlnPlot1.pdf"),width=10,height=20)
p <- VlnPlot(snATAC, features=genes, pt.size=0,ncol=1)
print(p)
dev.off()

# TSNE-FeaturePlot
pdf(paste0(tissue,"_20220621.snATAC.interesting_genes.FeaturePlot.pdf"),width=10,height=20)
p <- FeaturePlot(snATAC, order = T, features = genes,reduction = "tsne",pt.size = 0.1,by.col = FALSE,label=T) & theme(legend.position = "right")
print(p)
dev.off()
