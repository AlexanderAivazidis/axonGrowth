### Exploratory analysis on axon growth dataset:

## Load libraries:

library(readxl)
library(dplyr)
library(Seurat)
library(myUtils)
library(EnsDb.Mmusculus.v75)
library(lubridate)
library(ComplexHeatmap)

## Load and format Data:

directory = '/home/jovyan/axonGrowth/'
setwd(directory)
axonGrowth_counts = read.csv('counts.csv', row.names = 1)
rownames(axonGrowth_counts) = unlist(lapply(1:length(rownames(axonGrowth_counts)), function(x) strsplit(rownames(axonGrowth_counts)[x], split = '\\.')[[1]][1]))
symbols = mapIdsMouse(rownames(axonGrowth_counts), 'ENSEMBL', 'SYMBOL')
mitogenes <- genes(EnsDb.Mmusculus.v75, filter = ~ seq_name == "MT")$gene_id
percent.mt = colSums(axonGrowth_counts[rownames(axonGrowth_counts) %in% mitogenes,])/colSums(axonGrowth_counts)
rownames(axonGrowth_counts) = symbols
sum(unlist(lapply(1:dim(axonGrowth_counts)[1], function(x) substring(rownames(axonGrowth_counts)[x], 1,3))) != 'ENS')/dim(axonGrowth_counts)[1]
metadata = read_xlsx('Summary-neurons_microcoverslips_RNAseq.xlsx', col_types = c('text', 'text', 'date', 'date', 'date', 'date', 'date', rep('text', 7)), skip = 1)
metadata = metadata[10:105,]
axonGrowth = CreateSeuratObject(counts = axonGrowth_counts, project = 'axonGrowth', min.cells = 0, min.features = 0)
axonGrowth$TotalNeuriteLength = as.numeric(metadata$`Total Neurite Length (um)`)
axonGrowth$LongestNeuriteLength = as.numeric(metadata$`Longest Neurite Length (um)`)
axonGrowth$TimeSincePlating = unlist(lapply(1:dim(metadata)[1], function(x) interval(ymd_hms(metadata$`Time of plating (mm/dd/yyyy hh:mm)`[x]), ymd_hms(metadata$`Cell collection time (mm/dd/yyyy hh:mm)`[x]))/dhours(1)))
axonGrowth$TimeSincePlating = unlist(lapply(1:dim(metadata)[1], function(x) interval(ymd_hms(metadata$`Time of plating (mm/dd/yyyy hh:mm)`[x]), ymd_hms(metadata$`Cell collection time (mm/dd/yyyy hh:mm)`[x]))/dhours(1)))
axonGrowth$TimeInAIRMEM = unlist(lapply(1:dim(metadata)[1], function(x) interval(ymd_hms(metadata$`AIR-MEM change (mm/dd/yyyy hh:mm)`[x]), ymd_hms(metadata$`Cell collection time (mm/dd/yyyy hh:mm)`[x]))/dhours(1)))

## QC:
axonGrowth[["percent.mt"]] <- percent.mt
VlnPlot(axonGrowth, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(axonGrowth, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(axonGrowth, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
axonGrowth <- subset(axonGrowth, subset = nFeature_RNA > 1000 & nFeature_RNA < 20000 & percent.mt < 2.5 & nCount_RNA > 10^5 & nCount_RNA < 10^6)
VlnPlot(axonGrowth, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Normalize:
axonGrowth <- NormalizeData(axonGrowth, normalization.method = "LogNormalize", scale.factor = 10000)

## Highly variable features:
axonGrowth <- FindVariableFeatures(axonGrowth, selection.method = "vst", nfeatures = 200)

# Identify the 10 most highly variable genes:
top10 <- head(VariableFeatures(axonGrowth), 10)

# plot variable features with and without labels:
plot1 <- VariableFeaturePlot(axonGrowth)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# Scaling:
all.genes <- rownames(axonGrowth)
axonGrowth <- ScaleData(axonGrowth, features = all.genes)

# PCA:
axonGrowth <- RunPCA(axonGrowth, features = VariableFeatures(object = axonGrowth))

print(axonGrowth[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(axonGrowth, dims = 1:2, reduction = "pca")
DimPlot(axonGrowth, reduction = "pca")
DimHeatmap(axonGrowth, dims = 1, cells = round(dim(axonGrowth)[2])/2, balanced = TRUE)
DimHeatmap(axonGrowth, dims = 1:15, cells = round(dim(axonGrowth)[2])/2, balanced = TRUE)

# Determine Dimensionality of the dataset:

axonGrowth <- JackStraw(axonGrowth, num.replicate = 100)
axonGrowth <- ScoreJackStraw(axonGrowth, dims = 1:20)
JackStrawPlot(axonGrowth, dims = 1:20)
ElbowPlot(axonGrowth, ndims = 40)

# Clustering:
axonGrowth <- FindNeighbors(axonGrowth, dims = 1:5)
axonGrowth <- FindClusters(axonGrowth, resolution = 0.5)
head(Idents(axonGrowth), 5)

# Non-linear dimensional reduction:
axonGrowth <- RunUMAP(axonGrowth, dims = 1:5)
DimPlot(axonGrowth, reduction = "umap")

saveRDS(axonGrowth, file = "axonGrowth_SeuratTutorial.rds")

# Cluster biomarkers:
cluster1.markers <- FindMarkers(axonGrowth, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
axonGrowth.markers <- FindAllMarkers(axonGrowth, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
axonGrowth.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

VlnPlot(axonGrowth, features = c("Actb", "Tuba1a"))
VlnPlot(axonGrowth, features = c("ENSMUSG00000075015", "ENSMUSG00000095186"), slot = "counts", log = TRUE)

FeaturePlot(axonGrowth, features = c("Actb", "Tuba1a","ENSMUSG00000075015", "ENSMUSG00000095186"))


DoHeatmap(axonGrowth, features = top25$gene, cells = 1:dim(axonGrowth)[2] + NoLegend())

as.matrix(top25)

## Find Clusters or PCs that correlate with neurite length:

FeaturePlot(axonGrowth, features = c('TotalNeuriteLength', 'LongestNeuriteLength'))

correlationScores = data.frame(matrix(0,15,7))
colnames(correlationScores) = c('LongestNeuriteLength', 'TotalNeuriteLength', 'TimeSincePlating', 'TimeInAirMem', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')
  
for (i in 1:15){
  correlationScores[i,] = c(cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$LongestNeuriteLength, use = "complete.obs"), cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$TotalNeuriteLength, use = "complete.obs"),
          cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$TimeSincePlating, use = "complete.obs"),cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$TimeInAIRMEM, use = "complete.obs"),
          cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$nCount_RNA, use = "complete.obs"),cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$nFeature_RNA, use = "complete.obs"),
          cor(axonGrowth@reductions$pca@cell.embeddings[,i], axonGrowth$percent.mt, use = "complete.obs"))
}

temp = as.matrix(correlationScores)
rownames(temp) = paste('PCA', 1:15, sep = ' ')

Heatmap(temp, cluster_rows = FALSE)

head(sort(Loadings(axonGrowth, reduction = 'pca')[,c('PC_2', 'PC_3')], decreasing = TRUE), 25)

write.csv(sort(Loadings(axonGrowth, reduction = 'pca')[,c('PC_2', 'PC_3')], decreasing = TRUE), file = 'PC4_loadings.csv')

dev.off()
plot(axonGrowth@reductions$pca@cell.embeddings[,2], axonGrowth$TotalNeuriteLength, pch = '.', cex = 3, xlab = 'PCA2-score', ylab = 'TotalNeuriteLength')
plot(axonGrowth@reductions$pca@cell.embeddings[,2], axonGrowth$TotalNeuriteLength, pch = '.', cex = 3, xlab = 'PCA3-score', ylab = 'TotalNeuriteLength')
plot(axonGrowth@reductions$pca@cell.embeddings[,2], axonGrowth$LongestNeuriteLength, pch = '.', cex = 3, xlab = 'PCA2-score', ylab = 'LongestNeuriteLength')
plot(axonGrowth@reductions$pca@cell.embeddings[,3], axonGrowth$LongestNeuriteLength, pch = '.', cex = 3, xlab = 'PCA3-score', ylab = 'LongestlNeuriteLength')



