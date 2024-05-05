library(Seurat)
library(ArchR)
library(ggplot2)
library(patchwork)
library(cowplot)
library(Signac)
library(RColorBrewer)
library(Rhdf5lib)
library(parallel)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)

crc=readRDS("mm10.mt.RDS")
snp=read.table("select_snp")###final three SNP for split group Mutation and WT
DefaultAssay(crc) <- "alleles"
features=snp[,1]
cells <- colnames(x = crc)
data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = crc, slot = "data")[features, cells, drop = FALSE])))

G1=data[which(rowSums(data) > 0),]
G2=data[which(rowSums(data) == 0),]

c11=rownames(G1)
c12=rep(1,length(c11))
c1=cbind(c11,c12)
c21=rownames(G2)
c22=rep(0,length(c21))
c2=cbind(c21,c22)
Group=rbind(c1,c2)
crc@meta.data$SNP=Group[match(colnames(crc),Group[,1]),2]
saveRDS(crc, file="mm10.mt.add.SNP.rds" )
###
DefaultAssay(crc) <- "peaks"
pfm <- getMatrixSet(
           x = JASPAR2020,
           opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
crc <- AddMotifs(
                 object = crc,
                 genome = BSgenome.Mmusculus.UCSC.mm10,
                 pfm = pfm
)
crc <- RunChromVAR(
             object = crc,
             genome = BSgenome.Mmusculus.UCSC.mm10
        )
###For DEG analysis

DefaultAssay(crc) <- "peaks"
crc@meta.data$cell <- rownames(crc@meta.data)
group.A <- crc@meta.data %>% filter(SNP==1)
group.B <- crc@meta.data %>% filter(SNP==0)
logfold=log(1.25,2)##set fold change
# find peaks specific to one clone
da_peaks <- FindMarkers(
          object = crc,
          ident.1 =group.A$cell,
          ident.2 = group.B$cell,
          only.pos = TRUE,
          test.use = 'LR',
          min.pct = 0.05,
	  logfc.threshold = logfold,
          latent.vars = 'nCount_peaks'
)

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.01, ])
outfile <- "WT_vs_Mutation_da_peaks.xls"
write.table(top.da.peak, file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

enriched.motifs <- FindMotifs(
            object = crc,
            assay = "peaks",
            features = top.da.peak
          )
outfile <- "WT_vs_Mutation_da_peaks_enrich_motif.xls"
write.table(enriched.motifs, file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
DefaultAssay(crc) <- 'chromvar'
differential.activity <- FindMarkers(
             object = crc,
             ident.1 = group.A$cell,
             ident.2 = group.B$cell,
             only.pos = TRUE,
             min.pct = 0.05,
             logfc.threshold = 0,
             mean.fxn = rowMeans,
             fc.name = "avg_diff"
        )
outfile <- "WT_vs_Mutation_da_peaks_enrich_motif_differential.activity.xls"
write.table(differential.activity, file = outfile, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
pdf("top20_different_motif_plot.pdf")
p=MotifPlot(
           object = crc,
           motifs = head(rownames(differential.activity)),
           assay = 'peaks'
       )
print(p)
dev.off()



