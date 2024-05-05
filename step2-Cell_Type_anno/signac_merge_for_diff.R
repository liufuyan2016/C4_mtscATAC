#!/usr/bin/Rscript
parser = argparse::ArgumentParser(description="Script to clustering scATAC data")
parser = argparse::ArgumentParser(description=cat("Script to merge scATAC library, dimensionality reduction and clustering using Signac(reference pipeline https://satijalab.org/signac/articles/merging.html)\n\n# Authors:   Maven\n# Contact information:   mawen3@genomics.cn\n# Date: 2022-04-7\n# R package version: ggplot2_3.3.5        cowplot_1.1.1 data.table_1.14.0 future_1.23.0 GenomicRanges_1.44.0 GenomeInfoDb_1.28.4 IRanges_2.26.0       S4Vectors_0.30.2 BiocGenerics_0.38.0 SeuratObject_4.0.4 Seurat_4.0.5 Signac_1.4.0\n\n"))
parser$add_argument('-R','--report', help='input path of report file dir')
parser$add_argument('-s','--samples',default="Samples" ,help='sample name[default = \"%(default)s\"]')
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

tissue <- "all_sample"
Marker=args$marker
inputPath <-args$report
outdir <- args$out
sampleFile <- args$samples
dir.create(outdir)
########################################
setwd(outdir)
samples <- read.table(sampleFile,header=T)
tissue <- samples$tissue[1]




###################### computation ######################

################################ Create a gene activity matrix ################################
snATAC <- readRDS(inputPath)

DefaultAssay(snATAC) <- 'RNA'
Idents(snATAC) <- snATAC[["cellType"]]


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
Genelist$gene_id <- Genelist$Gene
da_peaks <- full_join(closest_genes.join,Genelist,by="gene_id")

write.table(da_peaks,paste0(tissue,".Table2.combined.harmony_by_dataset.da_peaks.20220527.combined_markers.xls"),sep = "\t",quote = FALSE,row.names = FALSE)
