#!/usr/bin/Rscript
parser = argparse::ArgumentParser(description="Script to clustering scATAC data")
parser = argparse::ArgumentParser(description=cat("Script to merge scATAC library, dimensionality reduction and clustering using Signac(reference pipeline https://satijalab.org/signac/articles/merging.html)\n\n# Authors:   Maven\n#  R package version: ggplot2_3.3.5        cowplot_1.1.1 data.table_1.14.0 future_1.23.0 GenomicRanges_1.44.0 GenomeInfoDb_1.28.4 IRanges_2.26.0       S4Vectors_0.30.2 BiocGenerics_0.38.0 SeuratObject_4.0.4 Seurat_4.0.5 Signac_1.4.0\n\n"))
parser$add_argument('-R','--report', help='input path of report file dir')
parser$add_argument('-s','--samples',default="Samples" ,help='sample name[default = \"%(default)s\"]')
parser$add_argument('-a','--ann',default=" " ,help='annotation for genome. for scmtATAC ,please use older version!!!')
parser$add_argument('-M','--marker', help='marker list')
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(stringr)
library(IRanges)
library(GenomicRanges)
library(harmony)
anno=args$ann
tissue <- "all_sample"
Marker=args$marker
inputPath <-args$report
outdir <- args$out
sampleFile <- args$samples
dir.create(outdir)
########################################
setwd(outdir)
samples <- read.table(sampleFile,header=T)
print(samples)
 cat(tissue,"has",nrow(samples),"libraries to merge.\n")
 peakSets = list()
 grSets =list()
 for(i in 1:nrow(samples)){
### 1. creating a common peak set ####
  sampleName <- samples$fullName[i]
  LibID <- samples$LibID[i]
  print (i)
  cat("\n\n","****** Creating a common peak set of the No.",i,"sample:",sampleName,"******\n")
  # read in peak sets
  peak_file <- paste0(inputPath,"/",LibID,"/out/",LibID,"_peaks.narrowPeak")
  peak_set <- as.data.frame(fread(peak_file,sep="\t",header = F))
  peak_set <- peak_set[,c(1,2,3)]
  colnames(peak_set) <- c("chr", "start", "end")
  cat("\n","Read in the peak bed file:",peak_file,"\n")
  peakSets[[sampleName]] <- peak_set
  print (sampleName)
  grSets[[sampleName]] <- makeGRangesFromDataFrame(peakSets[[sampleName]])
}
 ##########################################################################
 # Create a unified set of peaks to quantify in each dataset
print (grSets)
merged_grange=grSets[[1]]
for(i in 2:nrow(samples)){
    print (i)
   merged_grange <- c(merged_grange,grSets[[i]])
}
combined.peaks <- reduce(x=merged_grange,with.revmap=TRUE)
 # Filter out bad peaks based on length
 peakwidths <- width(combined.peaks)
 combined.peaks <- combined.peaks[peakwidths < 8000 & peakwidths > 20]
 print (combined.peaks)
 
 write.table(combined.peaks,paste0(tissue,".combined.peak_region.xls"),sep = "\t",quote = FALSE,row.names = FALSE)
 
 metaSets = list()
 fragsSets = list()
 mtxSets = list()
 countSets =list()
 ATACassay = list()
 snATACs = list()
 for(i in 1:nrow(samples)){
  tissue <- samples$tissue[i]
  sampleName <- samples$fullName[i]
  LibID <- samples$LibID[i]
  stage <- samples$stage[i]
  print (samples$LibID[i])
  cat("\n\n","****** Processing the No.",i,"sample:",sampleName,"******\n")
 ### 2. Create Fragment objects ###
  # load metadata
  peak_file <- paste0(inputPath,"/",LibID,"/out/",LibID,"_peaks.narrowPeak")
  cell_file <- paste0(inputPath,"/",LibID,"/d2cfile/","list.",LibID,".txt")
  fragment_file <- paste0(inputPath,"/",LibID,"/d2cfile/",LibID,".fragments.tsv.gz")
  cat("\n","Read in the peak count file:",peak_file,"\n")
  peak_set <- as.data.frame(fread(peak_file,sep="\t",header = F))
  peak_set <- peak_set[,c(1,2,3)]
  colnames(peak_set) <- c("chr", "start", "end")
  grfile <- makeGRangesFromDataFrame(peak_set)
  cellsFile <- read.table(cell_file,header=F)
  frag <- CreateFragmentObject(path=fragment_file,cells=cellsFile$V1)
  my_counts <- FeatureMatrix(fragments = frag,features = grfile, cells=cellsFile$V1,sep=c(":","-"))
  countSets[[sampleName]] = my_counts
  
  m <- paste0(inputPath,"/",samples$LibID[i],"/d2cfile/",LibID,".Metadata.tsv")
  cat("\n","Read in the meta file:",m,"\n")
  raw_meta=read.table(m,header=T,comment.char="")
  rownames(raw_meta)=raw_meta[,1]
  raw_meta=raw_meta[,-1]
  my_meta=raw_meta
  colnames(my_meta)[2]="uniqueNuclearFrags"
  my_counts <- my_counts[grep(paste0("^","chrM"),rownames(my_counts),invert=T),]
  my_meta$FRIP=colSums(my_counts)/my_meta$uniqueNuclearFrags 
  metaSets[[sampleName]]= my_meta
  
  # perform an initial filtering of low count cells
  metaSets[[sampleName]] <- metaSets[[sampleName]][metaSets[[sampleName]]$totalFrags > 500, ]
  
  # create fragment objects
  m <- paste0(inputPath,"/",samples$LibID[i],"/d2cfile/",LibID,".fragments.tsv.gz")
  cat("\n","Read in the fragment file:",m,"\n")
  fragsSets[[sampleName]] <- CreateFragmentObject(
  path = m,
  cells = rownames(metaSets[[sampleName]])
  )
 ### 3. Quantify peaks in each dataset ###
  mtxSets[[sampleName]] <- FeatureMatrix(
  fragments = fragsSets[[sampleName]],
  features = combined.peaks,
  cells = rownames(metaSets[[sampleName]])
  )
 ### 4. Create the objects ###
  ATACassay[[sampleName]] <- CreateChromatinAssay(
    counts=countSets[[sampleName]],
	sep = c(":", "-"),
	#genome=Genome,
	fragments = fragsSets[[sampleName]],
	min.cells = 0,
	min.features = 0
	)
  snATACs[[sampleName]] <- CreateSeuratObject(
   counts = ATACassay[[sampleName]], 
   assay = "peaks", 
   meta.data=metaSets[[sampleName]],
   min.cells = 10
  )
  # add information to identify dataset of origin
  snATACs[[sampleName]]$dataset <- as.character(sampleName)
  snATACs[[sampleName]]$tissue <- tissue
  print (stage)
  snATACs[[sampleName]]$stage <- as.character(stage)
  }

### 5. merge objects ###
#############################################################################
  # merge all datasets, adding a cell ID to make sure cell names are unique
  combined2 <- list()
  for(i in 2:nrow(samples)){
  combined2[[i-1]] <- snATACs[[i]]
}

combined <- merge(
   x = snATACs[[1]],
   y = combined2,
   )
 combined[["peaks"]] 
### add annotation   

anno <- as.data.frame(fread(anno))
txdbNew <- makeGRangesFromDataFrame(anno, keep.extra.columns=TRUE, start.field="Start",end.field="End")
Annotation(combined) <- txdbNew

###################### computation ######################
setwd(outdir)

combined <- NucleosomeSignal(object = combined)
combined <- TSSEnrichment(object = combined, fast = FALSE)
saveRDS(combined, file=paste0(tissue,"_combined.TSSEnrichment.rds") )

pdf(file=paste0(tissue,"_TSS_enrichment_plot.pdf"),width=15,height=8)
combined$high.tss <- ifelse(combined$TSS.enrichment > 1.5, 'High', 'Low')
p <- TSSPlot(combined, group.by = 'high.tss')
print(p)
p <- TSSPlot(combined, group.by = 'dataset')
print(p)
p <- TSSPlot(combined, group.by = 'stage')
print(p)
dev.off()

  
# 计算log10_unique Fragments
combined@meta.data$log10_uniqueFrags=log10(combined@meta.data$uniqueNuclearFrags)

pdf(file=paste0(tissue,"_preQC_VlnPlot.pdf"),width=18,height=8)
p <- VlnPlot(combined, group.by="tissue",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0,ncol = 5)&geom_boxplot(width=.1,col="black",outlier.shape = NA,fill="white")
print(p)
p <- VlnPlot(combined, group.by="dataset",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0,ncol = 5)&geom_boxplot(width=.1,col="black",outlier.shape = NA,fill="white")
print(p)
p <- VlnPlot(combined, group.by="stage",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0,ncol = 5)&geom_boxplot(width=.2,col="black",outlier.shape = NA,fill="white")
print(p)
dev.off()

saveRDS(combined, file=paste0(tissue,"_combined.preQC.rds") )

# 过滤
combined <- subset(
x = combined,
subset = uniqueNuclearFrags > 2000 &
uniqueNuclearFrags < 50000 &
TSS.enrichment > 1.2 &
FRIP > 0.2
)

pdf(file=paste0(tissue,"_afterQC_VlnPlot.pdf"),width=18,height=8)
p <- VlnPlot(combined, group.by="tissue",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0,ncol = 5)&geom_boxplot(width=.1,col="black",outlier.shape = NA,fill="white")
print(p)
p <- VlnPlot(combined, group.by="dataset",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0,ncol = 5)&geom_boxplot(width=.1,col="black",outlier.shape = NA,fill="white")
print(p)
p <- VlnPlot(combined, group.by="stage",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0,ncol = 5)&geom_boxplot(width=.1,col="black",outlier.shape = NA,fill="white")
print(p)
dev.off()

####################### Normalization and linear dimensional reduction #######################
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 50)### min.cutoff = 20) #include features in >20 cells in the set of VariableFeatures
# 奇异值分解
combined <- RunSVD(combined)
saveRDS(combined, file=paste0(tissue,"_combined.SVD.rds"))

# 
pdf(file=paste0(tissue,"_Correlation_between_depth_dimensions.pdf"))
p <- DepthCor(combined)
print(p)
dev.off()

####################### Non-linear dimension reduction and clustering ####################### 
combined <- RunTSNE(object = combined, reduction = 'lsi', dims=2:20)  # 这里保留dim多一些比较好
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims=2:20)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
saveRDS(combined, file=paste0(tissue,"_combined.FindClusters.rds") )

pdf(file=paste0(tissue,"_combined_TSNE.pdf"))
p <- DimPlot(object = combined, label = TRUE)+ggtitle(label = "tsne_Cluster.combined")
print(p)
p <- DimPlot(object = combined, group.by = "tissue", label = TRUE)+ggtitle(label = "tsne_Cluster.combined.by tissue")
print(p)
p <- DimPlot(object = combined, group.by = "stage", label = TRUE)+ggtitle(label = "tsne_Cluster.combined.by stage")
print(p)
p <- DimPlot(object = combined, group.by = "dataset", label = TRUE)+ggtitle(label = "tsne_Cluster.combined.by dataset")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_TSNE_LibSplit.pdf"),width=25,height=7)
p <- DimPlot(combined, reduction="tsne", split.by = "dataset",pt.size = 0.2,label=T)+ggtitle(label = "tsne.LibSplit.combined")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_TSNE_stageSplit.pdf"),width=14,height=7)
p <- DimPlot(combined, reduction="tsne", split.by = "stage",pt.size = 0.2,label=T)+ggtitle(label = "tsne.StageSplit.combined")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_TSNE_FeaturePlot.tissue.pdf"),width=36,height=7)
p <- FeaturePlot(combined, split.by="tissue",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0.1,by.col = FALSE,label=T)+ggtitle(label = "tsne_Cluster.combined.by tissue") & theme(legend.position = "right")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_TSNE_FeaturePlot.stage.pdf"),width=20,height=6)
p <- FeaturePlot(combined, split.by="stage", features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0.1,by.col = FALSE)+ggtitle(label = "tsne_Cluster.combined.by stage") & theme(legend.position = "right")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_TSNE_FeaturePlot.Libs.pdf"),width=20,height=12)
p <- FeaturePlot(combined, split.by="dataset", features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0.1,by.col = FALSE)+ggtitle(label = "tsne_Cluster.combined.by dataset") & theme(legend.position = "right")
print(p)
dev.off()

############ harmony去批次 ############
combined <- readRDS(paste0(tissue,"_combined.SVD.rds"))
combined <- RunTSNE(object = combined, reduction = 'lsi', dims=2:30) # 这里保留dim多一些比较好
harmony <- RunHarmony(
  object=combined,
  group.by.vars='dataset',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
saveRDS(harmony, file=paste0(tissue,"_combined.harmony.rds"))

# re-compute the TSNE using corrected LSI embeddings
harmony <- RunTSNE(harmony, dims = 2:30, reduction = 'harmony')
harmony <- FindNeighbors(object = harmony, reduction = 'harmony', dims=2:20) 
res.used <- seq(0.1,1,by=0.1)
res.used
# Loop over and perform clustering of different resolutions 
pdf(file=paste0(tissue,"_combined.harmony_TSNE.pdf"),width=8,height=6)
for(i in res.used){
   harmony <- FindClusters(object = harmony, verbose = T, algorithm = 3,resolution = i)
   p <- DimPlot(object = harmony, label = TRUE)+ggtitle(label = paste0(i,"tsne_Cluster.harmony"))
   print(p)
   p <- DimPlot(object = harmony, group.by = "tissue", label = TRUE)+ggtitle(label = paste0(i,"tsne_Cluster.combined_harmony.by tissue"))
   print(p)
   p <- DimPlot(object = harmony, group.by = "stage", label = TRUE)+ggtitle(label = paste0(i,"tsne_Cluster.combined_harmony.by stage"))
   print(p)
   p <- DimPlot(object = harmony, group.by = "dataset", label = TRUE)+ggtitle(label = paste0(i,"tsne_Cluster.combined_harmony.by dataset"))
   print(p)
}
dev.off()

# Make plot 
library(clustree)
clus.tree.out <- clustree(harmony) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
pdf(file=paste0(tissue,"_combined.harmony_TSNE_resolution_tree.pdf"),width=20,height=10)
print(clus.tree.out)
dev.off()
saveRDS(harmony, file=paste0(tissue,"_combined.harmony_FindClusters.rds") )

pdf(file=paste0(tissue,"_combined_.harmonyTSNE_LibSplit.pdf"),width=25,height=7)
p <- DimPlot(harmony, reduction="tsne", split.by = "dataset",pt.size = 0.2,label=T)+ggtitle(label = "tsne.Libsplit.harmony")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_.harmonyTSNE_stageSplit.pdf"),width=14,height=7)
p <- DimPlot(harmony, reduction="tsne", split.by = "stage",pt.size = 0.2,label=T)+ggtitle(label = "tsne.StageSplit.harmony")
print(p)
dev.off()

pdf(paste0(tissue,"_combined.harmony_TSNE_split.pdf"), length(levels(harmony))*5, 6)
p <- DimPlot(harmony, reduction="tsne", split.by = "seurat_clusters",pt.size = 0.1,label=F)+ggtitle(label = "tsne.split.harmony")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_harmony_TSNE_FeaturePlot.tissue.pdf"),width=36,height=7)
p <- FeaturePlot(harmony, split.by="tissue",features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0.1,by.col = FALSE,label=T)+ggtitle(label = "tsne_Cluster.combined_harmony.by tissue") & theme(legend.position = "right")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_harmony_TSNE_FeaturePlot.stage.pdf"),width=20,height=6)
p <- FeaturePlot(harmony, split.by="stage", features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0.1,by.col = FALSE)+ggtitle(label = "tsne_Cluster.combined_harmony.by stage") & theme(legend.position = "right")
print(p)
dev.off()

pdf(file=paste0(tissue,"_combined_harmony_TSNE_FeaturePlot.Libs.pdf"),width=20,height=12)
p <- FeaturePlot(harmony, split.by="dataset", features = c('log10_uniqueFrags','TSS.enrichment', 'nucleosome_signal','tssProportion','FRIP'),pt.size = 0.1,by.col = FALSE)+ggtitle(label = "tsne_Cluster.combined_harmony.by dataset") & theme(legend.position = "right")
print(p)
dev.off()

################################ Create a gene activity matrix ################################
snATAC <- readRDS(paste0(tissue,"_combined.harmony_FindClusters.rds") )
# Create a gene activity matrix
gene.activities <- GeneActivity(snATAC,extend.upstream = 2000)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
snATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
snATAC <- NormalizeData(
object = snATAC,
assay = 'RNA',
normalization.method = 'LogNormalize',
scale.factor = median(snATAC$nCount_RNA)
)
saveRDS(snATAC, file =paste0(tissue,"_combined.harmony_GeneActivityLogNormalize.rds"))

DefaultAssay(snATAC) <- 'RNA'

################################ plot marker genes ################################
inGenelist <- Marker
Genelist <- read.table(inGenelist,sep="\t",stringsAsFactors=F,header=F)
if(ncol(Genelist)==1){
        Genelist <- Genelist[which(Genelist$V1 %in% rownames(snATAC@assays$RNA@counts)),]
        genes <- as.character(Genelist)
        genes_title <- as.character(Genelist)
        }else{
        Genelist <- Genelist[which(Genelist$V1 %in% rownames(snATAC@assays$RNA@counts)),]  
        genes <- as.character(Genelist$Gene)
        genes_title <- as.character(Genelist$V2)
        }
# plot markers #
genes <- unique(Genelist)
pdf(paste0(tissue,"_combined.harmony.combined_markers.DotPlot.pdf"),width=80,height=10)
p <- DotPlot(snATAC, features=genes, group.by='seurat_clusters')+RotatedAxis()
print(p)
dev.off()

# 小提琴图
pdf(paste0(tissue,"_combined.harmony.combined_markers.VlnPlot.pdf"),width=20,height=150)
p <- VlnPlot(snATAC, features=genes, pt.size=0)
print(p)
dev.off()
