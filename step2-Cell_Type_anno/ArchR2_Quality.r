parser = argparse::ArgumentParser(description=cat("Script to scATAC_ArchR_analysis  \n\n# Authors:   Maven\n#  Date: Aug 1,2022\n\n"))
parser$add_argument('-F', '--fragments', help='input directory of fragment file, .tbi files must including this directory')
parser$add_argument('-D', '--metadata', help='input metadata files')
parser$add_argument('-S', '--sample', help='sample name')
parser$add_argument('-P', '--species', help='species name:hg19 mm10')
parser$add_argument('-m', '--min', help='minimum fragment number')
parser$add_argument('-M', '--marker', help='marker gene')
parser$add_argument('-H', '--Harmony', help='Harmony')
parser$add_argument('-T', '--tsne', help='TSNE')
parser$add_argument('-R', '--res', help='res')
parser$add_argument('-O', '--output', default='.', help='Out directory for results [default = .]')

args = parser$parse_args()

library(ArchR)
library(ggplot2)
library(patchwork)
library(cowplot)
library(Signac)
library(RColorBrewer)
library(Rhdf5lib)
library(parallel)
library(Seurat)

set.seed(1)
addArchRThreads(threads = 10)
addArchRGenome(args$species)
w=70
h=100
sample=args$sample
if ( is.null(args$tsne) )  {
	method="UMAP"
        method2="umap"
}else{
	method="TSNE"
	method2="tsne"
}
outputDirectory <- args$output
if ( is.null(args$marker) )  {
	print ("use top10 marker")
}else{
Marker=read.table(args$marker,header=FALSE);
markerGenes  <- Marker$V1
}

if ( is.null(args$min)){
   minF=4000
}else{
	minF=as.numeric(args$min)
}
TSSEnrichment=5###quality important

### input frag
inputdir=args$fragments
first_name <- list.files(inputdir,"tsv.gz$")
sample_name <- gsub(".fragments.tsv.gz","",first_name)
FragmentFiles <- paste(inputdir,first_name,sep="/")
names(FragmentFiles) <- sample_name
rm(first_name)

## create ArrowFiles
print("begin to create ArrowFiles")
ArrowFiles <- createArrowFiles(
 inputFiles = FragmentFiles,
 sampleNames = sample_name,
 minTSS = TSSEnrichment, #Dont","set","this","too","high","because","you","can","always","increase","later
 minFrags = minF,
 addTileMat = TRUE,
 addGeneScoreMat = TRUE,
 #geneAnnotation = geneAnnotation,
 #genomeAnnotation = genomeAnnotation,
 offsetPlus = 0,
 offsetMinus = 0,
)
print("creating ArrowFiles has finished")

## calculate DoubletScore
print("begin to calculate DoubletScore")
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers","to","how","many","cells","near","a "pseudo-doublet" to","count.
    knnMethod = "UMAP", #Refers","to","the","embedding","to","use","for","nearest","neighbor","search","with","doublet","
    LSIMethod = 1
)
print("calculating DoubletScore has finished")
## create ArchR object
print("begin to create ArchR object")

raw_proj <- ArchRProject(
   ArrowFiles = ArrowFiles,
   outputDirectory = outputDirectory,
   copyArrows = T,
  # geneAnnotation = geneAnnotation,
  # genomeAnnotation = genomeAnnotation
 )

print("creating ArchR object has finished")
print("begin to filter cells which nFrags < 3000 & TSSEnrichment < 4 and doublet")
## filter cells which nFrags < 3000 & TSSEnrichment < 2 
proj <- raw_proj[which(raw_proj$nFrags > minF & raw_proj$TSSEnrichment >TSSEnrichment),]
## filter doublet cells
proj <- filterDoublets(proj, filterRatio = 2)
remain_cell=rownames(proj@cellColData)
write.csv(remain_cell,file=paste0(outputDirectory,"/",sample,"final_used_cells",".csv"))
p1 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p5 <- plotFragmentSizes(ArchRProj = proj)
p6 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p5,p6, name = paste0(outputDirectory,"/",sample,"QC-Sample-FragSizes-TSSProfile.pdf"), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
plotPDF(p1,p2,p3,p4, name =paste0(outputDirectory,"/",sample, "-QC-Sample-Statistics.pdf"), ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

saveArchRProject(ArchRProj = proj, outputDirectory = outputDirectory, load = FALSE)


### input frag
#proj<- loadArchRProject(outputDirectory)

## Dimensionality reduction, clustering and visualization
res=as.double(args$res)
print(res)
print("begin dimensionality reduction, clustering and visualization")
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 5,
    clusterParams = list( #See","Seurat::FindSample
        resolution = c(res),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 35000,
    dimsToUse = 1:50,
    corCutOff = 0.7,
    force = T
)
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = res,
    force = T
)
if(method=="UMAP"){
proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = method,
    nNeighbors = 50,
    minDist = 0.5,
    metric = "cosine", force = T
)
}else{
	proj <- addTSNE(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = method,
    nNeighbors = 50,
    minDist = 0.5,
    metric = "cosine", force = T
)
}
cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
print(cM)
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
ggsave(p,file=paste0(outputDirectory, "/",sample,res,".cluster_distribution.heatmap.pdf"))


proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",force = T
)
proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Harmony_Clusters",
    resolution = 0.8,force = T
)
if(method=="UMAP"){
   proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 50,
    minDist = 0.5,
    metric = "cosine",force = T
   )
   x="UMAPHarmony"
}else{
    proj <- addTSNE(
    ArchRProj = proj,
    reducedDims = "Harmony",
    name = "TSNEHarmony",
    nNeighbors = 50,
    minDist = 0.5,
    metric = "cosine",force = T
   )
   x="TSNEHarmony"
}
write.csv(proj@cellColData,file=paste0(outputDirectory,"/",sample,res,"all.infor_cell_level",".csv"),quote =FALSE)

H=args$Harmony
if(H =="yes"){
   Embedding1="Harmony"
   Embedding2=x
   count1=table(proj@cellColData$Harmony_Clusters)
   count=as.data.frame(count1)
   count_remain=count[which(count$Freq>30),1]
   proj=proj[which(proj@cellColData$Harmony_Clusters %in% count_remain),]
  # proj@cellColData$Clusters=proj@cellColData$Harmony_Clusters
   ClusterName="Harmony_Clusters"

}else{
    Embedding2=method
    Embedding1="IterativeLSI"
    ClusterName="Clusters"
}
print("dimensionality reduction, clustering and visualization have been finished")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = Embedding2)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = ClusterName, embedding = Embedding2)
plotPDF(p1,p2, name = paste0(outputDirectory,"/",sample,"-Plot-",method,"-Sample-Clusters.pdf"), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
saveRDS(proj, paste0(outputDirectory, "/",sample,res,"_filtered.rds"))
####
print("begin to get marker gene score of subtype")
AS_markers <- c()
AS_markerList <- c()
groupby=ClusterName ##could add group keys,such as subtype
AS_markers <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = groupby,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

AS_markerList <- getMarkers(AS_markers, cutOff = "FDR <= 0.01 & Log2FC >= 1",n=50)
write.csv(AS_markerList,file=paste0(outputDirectory,"/",sample,"_GeneScore_gene_marker_",res,".csv"))
save(AS_markers,file=paste0(outputDirectory,"/",sample,"_GeneScore_gene_marker_",res,".Rdata"))
AS_markerList2 <- getMarkers(AS_markers, cutOff = "FDR <= 0.01 & Log2FC >= 1",n=10)
if ( is.null(args$marker) )  {
	listgene=unlist(AS_markerList2)
	markerGenes  <- listgene$name
}else{
        listgene=unlist(AS_markerList)
	markerGenes=intersect(markerGenes,listgene$name)
        print("use give marker ")
        print(markerGenes)
}


heatmapGS <- plotMarkerHeatmap(
  seMarker = AS_markers,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  labelMarkers = markerGenes,
 transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = paste0(outputDirectory,"/",sample,res,"GeneScores-Marker-Heatmap"), width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

heatmapGS <- plotMarkerHeatmap(
  seMarker = AS_markers,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = paste0(outputDirectory,"/",sample,res,"GeneScores-known-Marker-Heatmap"), width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

proj <- addImputeWeights(proj)
p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = Embedding2,
 imputeWeights = getImputeWeights(proj)
)
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 5),p2))
plotPDF(plotList = p,
    name = paste0(outputDirectory,"/",sample,res,"-Plot-Marker-Genes.pdf"),
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)
#p <- plotBrowserTrack(
#    ArchRProj = proj,
#    groupBy = groupby,
#    geneSymbol = markerGenes,
#    upstream = 50000,
#    downstream = 50000
#)

#plotPDF(plotList = p,
#    name = paste0(outputDirectory,"/",sample,res,"-Plot-Tracks-Marker-Genes-groupbyCluster.pdf"),
#    ArchRProj = proj,
#    addDOC = FALSE, width = 5, height = 5)



p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "Sample",
    geneSymbol = markerGenes,
    upstream = 100000,
    downstream = 100000
)
plotPDF(plotList = p,
    name = paste0(outputDirectory,"/",sample,res,"-Plot-Tracks-Marker-Genes-groupbySample.pdf"),
    ArchRProj = proj,addDOC = FALSE, width = 14, height = 7)




#####archr2surat
archr2surat= getMatrixFromProject(proj,useMatrix = "GeneScoreMatrix")
ipc_as_mat <- assay(archr2surat)
ipc_as_mat <- ipc_as_mat[, proj$cellNames]
colnames(ipc_as_mat) <- proj$cellNames
rownames(ipc_as_mat) <- archr2surat@elementMetadata$name
# 2 Extracting (harmony corrected) LSI
ipc_arr <- getReducedDims(proj,
                           reducedDims = Embedding1,##Harmony 
                           corCutOff = 0.75,
                           scaleDims = FALSE)
row.names(ipc_arr) <- proj$cellNames
ipc_metad <- as.data.frame(getCellColData(proj))
ipc_metad <- as.data.frame(getCellColData(proj))
row.names(ipc_metad) <- proj$cellNames

umap <- as.matrix(getEmbedding(proj,embedding=Embedding2))##UMAPHarmony
row.names(umap) <- proj$cellNames
colnames(umap) <- c( paste0(method,1), paste0(method,2))

# tsne <- as.matrix(getEmbedding(mouse_atac,embedding="HarmonyTSNE"))
# row.names(tsne) <- mouse_atac$cellNames
# colnames(tsne) <- c("TSNE1", "TSNE2")

 atac.seurat_gene <- CreateSeuratObject(counts = ipc_as_mat, meta.data = ipc_metad, assay = "GeneScoreScaled")
# atac.seurat_gene[["harmony"]] <- CreateDimReducObject(embeddings = ipc_harmony, key = "harmony_", assay = "GeneScoreScaled")
 atac.seurat_gene[[method2]] <- CreateDimReducObject(embeddings = method2, key = paste0(method2,"_"), assay = "GeneScoreScaled")
# atac.seurat_gene[["tsne"]] <- CreateDimReducObject(embeddings = tsne, key = "tsne_", assay = "GeneScoreScaled")


#pdf(paste0(outputDirectory,"/",sample,res,"_marker-archR2seyrat.pdf"),width = w, height = h)
#FeaturePlot(atac.seurat_gene, features = markerGenes,slot = "data",cols = c("lightgrey", "purple"), ncol = 3, pt.size = 0.1)
#dev.off()
#pdf(paste0(outputDirectory,"/",sample,res,"_marker-archR2seyrat.vlnplot.pdf"),width = w, height = h)
#VlnPlot(atac.seurat_gene, features = markerGenes,pt.size = 0,slot = "data",ncol = 3,group.by =ClusterName )
#dev.off()
#pdf(paste0(outputDirectory,"/",sample,res,"_marker-archR2seyrat.Dotplot.pdf"),width = w, height = 5)
#DotPlot(object=atac.seurat_gene, features = markerGenes,group.by = ClusterName)
#dev.off()
saveRDS(atac.seurat_gene,file=paste0(outputDirectory,"/",sample,res,"_archR2seyrat.rds"))



#metadata <- as.data.frame(fread(args$metadata))
#metadata$cells <- sapply(strsplit(metadata$V1, "#"),"[",2)
#proj@cellColData$cellnames <- sapply(strsplit(rownames(proj@cellColData),"#"),"[",2)
#proj@cellColData$subtype <- metadata[match(proj@cellColData$cellnames,metadata$cells),"Group_subtype"]
print("begin to get marker peaks of subtype")

proj <- addGroupCoverages(
  ArchRProj = proj ,
  groupBy = groupby,
  minCells = 50, 
  maxCells = 1000, ## maximum number of cells per replicate
  maxFragments = 50*1e6, ## maximum number of fragments per cell group
  minReplicates = 2, ## minimum number of replicates per group
  maxReplicates = 5, ## maximum number of replicates per group
  sampleRatio = 0.8, ## f
  useLabels = TRUE,
  force = TRUE,
  thread = 1
)
pathToMacs2 <- "/jdfsbjcas1/ST_BJ/PUB/Pipeline/DNBelabC4scATAC/V2.1/bin/macs2"
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = groupby,
    pathToMacs2 = pathToMacs2,
    #genomeAnnotation = genomeAnnotation,
    #geneAnnotation = geneAnnotation,
    #genomeSize = 2.5e9,
    reproducibility = "2", 
    peaksPerCell = 1000,
    maxPeaks = 200000,
    minCells = 40,
    excludeChr = c("chrY", "chrMT"), 
    method = "q",
    cutOff = 0.01, 
    extendSummits = 250,
    force = T
)

proj2 <- addPeakMatrix(proj)
saveArchRProject(ArchRProj = proj2, outputDirectory = outputDirectory, load = FALSE)

peaks=getPeakSet(proj2)
print(peaks)
saveRDS(proj2, paste0(outputDirectory, "/",sample,res,"_addpeaks.rds"))

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj2,
    useMatrix = "PeakMatrix",
    groupBy =groupby ,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",

)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")



write.csv(markerList,file=paste0(outputDirectory, "/",sample,"_Clusters_marker_peaks",res,".csv"))
saveRDS(markersPeaks, paste0(outputDirectory, "/",sample,"_marker_peaks.rds"))

write.csv(proj2@peakSet,file=paste0(outputDirectory,"/",sample,res,"all.infor_peak_level",".csv"),quote =FALSE)

eatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
plotPDF(heatmapPeaks, name = paste0(outputDirectory, "/",sample,"_Peak-Marker-Heatmap"), width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)



pv <- plotMarkers(seMarker = markersPeaks, name = "Erythroid", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = paste0(outputDirectory, "/",sample,"_Erythroid-Markers-MA-Volcano"), width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)
#p <- plotBrowserTrack(
#    ArchRProj = proj2, 
#    groupBy = "Clusters2", 
#    geneSymbol = c("GATA1"),
#    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
#    upstream = 50000,
#    downstream = 50000
#)
#plotPDF(p, name = paste0(outputDirectory, "/",sample,"_Plot-Tracks-With-Features"), width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)




###############
print ("motif analysis")
#proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", name = "Motif")
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
saveRDS(proj2, paste0(outputDirectory, "/",sample,res,"_addmotifs.rds"))

enrichMotifs  <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)



df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
plotPDF(ggUp, ggDo, name = paste0(outputDirectory, "/",sample,res,"Erythroid-vs-Progenitor-Markers-Motifs-Enriched"), width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = paste0(outputDirectory, "/",sample,res,"_Subtype-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)
