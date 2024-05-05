parser = argparse::ArgumentParser(description=cat("Script to scATAC_ArchR_analysis  \n\n# Authors:   Maven\n# Contact information: # Date: Aug 1,2022\n\n"))
parser$add_argument('-m', '--mgatk', help='input mgatk dir')
parser$add_argument('-S', '--sample', help='samples for analysis ')
parser$add_argument('-N', '--SNP', help='SNP order file')
parser$add_argument('-I', '--index', help='index for output')
parser$add_argument('-R', '--RDS', help='RDS')
parser$add_argument('-O', '--output', default='.', help='Out directory for results [default = .]')

args = parser$parse_args()
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(Signac)
library(RColorBrewer)
library(Rhdf5lib)
library(parallel)
set.seed(1)
w=70
h=100

sample=args$sample
oldsample_data=read.table(sample,header=FALSE);
oldsample=oldsample_data$V1
outputDirectory <- args$output
if (! is.null(args$SNP ) ){
	SNP=read.table(args$SNP,header=FALSE);
}
sample=args$index
### change cell
crc1=readRDS(args$RDS)
info=crc1@meta.data
write.table(info,file=paste0(outputDirectory,"/",sample,"cell_anno.xls"),sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
crc=crc1[,crc1@meta.data$dataset %in% oldsample]

mito.data <- ReadMGATK(dir = args$mgatk)
mito <- CreateAssayObject(counts = mito.data$counts,cells= row.names(mito.data$depth))

mito <- subset(mito, cells = colnames(crc))
crc[["mito"]] <- mito
crc <- AddMetaData(crc, metadata = mito.data$depth, col.name = "mtDNA_depth")
pdf(paste0(outputDirectory,"/",sample,"depth.vlnplot.pdf"),width = 6, height = 4)
VlnPlot(crc, "mtDNA_depth", pt.size = 0.1,group.by = "cellType") + scale_y_log10()
dev.off()

variable.sites <- IdentifyVariants(crc, assay = "mito",low_coverage_threshold = 10,  refallele = mito.data$refallele)

pdf(paste0(outputDirectory,"/",sample,".filter_SNP_QC.pdf"),width = 6, height = 6)
p=VariantPlot(variants = variable.sites,min.cells = 15)
print(p)
dev.off()

high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 50 
)

high=high.conf[,c(1,2,5)]
plotsnp=rownames(high)
print(high)
write.table(high,file=paste0(outputDirectory,"/",sample,"SNP.filter.xls"),sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)
crc[["alleles"]]
DefaultAssay(crc) <- "alleles"

##use the order of SNP id #######################

if(! is.null(args$SNP ) ){
   alleles.view <-SNP[1:10,1]
}

#################################################

alleles.view <- plotsnp[1:10]## c("16012G>A", "4103G>A","15930C>T","10075G>A")
pdf(paste0(outputDirectory,"/",sample,".UMAP.mt.pdf"),width = 13, height = 8)
FeaturePlot(
  object = crc,
  features = alleles.view,
  order = TRUE,
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  cols = c("grey", "darkred"),
  ncol = 4
) 
dev.off()

pdf(paste0(outputDirectory,"/",sample,".mt.heatmap.pdf"),width = 10, height = 13)
DoHeatmap(crc, features = rownames(crc),group.by = "cellType",slot = "data", disp.max =0.05,angle = 90,size=3) +scale_fill_viridis_c()
dev.off()

DefaultAssay(crc) <- "alleles"
crc <- FindClonotypes(crc)
table(Idents(crc))
pdf(paste0(outputDirectory,"/",sample,".mt.heatmap2.pdf"),width = 6, height = 13)
DoHeatmap(crc, features = VariableFeatures(crc), slot = "data", disp.max = 0.05,angle = 90,size=3) +
  scale_fill_viridis_c()
dev.off()
###produce data
features=VariableFeatures(crc)
features <- rev(x = unique(x = features))
cells <- colnames(x = crc)
data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = crc, slot = "data")[features, cells, drop = FALSE])))
write.table(t(data),file=paste0(outputDirectory,"/",sample,".mt_SNP.txt"),quote=FALSE)

G1=data[which(rowSums(data) > 0),]
G2=data[which(rowSums(data) == 0),]

c11=rownames(G1)
c12=rep(1,length(c11))
c1=cbind(c11,c12)
c21=rownames(G2)
c22=rep(0,length(c21))
c2=cbind(c21,c22)
Group=rbind(c1,c2)
write.table(Group,file=paste0(outputDirectory,"/",sample,".mt_clone.txt"),row.names=FALSE,quote=FALSE)
crc@meta.data$SNP=Group[match(colnames(crc),Group[,1]),2]
saveRDS(crc, file=paste0(outputDirectory,"/",sample,".mt.RDS") )

DefaultAssay(crc) <- "peaks"

# find peaks specific to one clone
markers.fast <- FoldChange(crc, ident.1 = 2)
markers.fast <- markers.fast[order(markers.fast$avg_log2FC, decreasing = TRUE), ] # sort by fold change
write.table(markers.fast,file=paste0(outputDirectory,"/",sample,".mt_clone.diff.txt"),row.names=TRUE,quote=FALSE)
head(markers.fast)
pdf(paste0(outputDirectory,"/",sample,".mt.coverage.pdf"),width = 6, height = 24)
CoveragePlot(
  object = crc,
  region = rownames(markers.fast)[1],
  extend.upstream = 2000,
  extend.downstream = 2000
)
dev.off()

