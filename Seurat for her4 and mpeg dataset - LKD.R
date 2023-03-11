library(ggplot2)
library(dplyr)
library(Seurat)
library(cowplot)
library(scales)

##for mpeg dataset

mpegn.data<-Read10X(data.dir="11609BC0003TXP-N__MPEG_Na_ve/filtered_feature_bc_matrix/") #replace with path to naive dataset folder
mpegl.data<-Read10X(data.dir = "11609BC0004TXP-N__MPEG_Lesioned/filtered_feature_bc_matrix/") #replace with path to Lesioned dataset folder

mpegn<-CreateSeuratObject(counts = mpegn.data, project = "mpegnaive", min.cells = 3, min.features = 200, assay="RNA")
mpegl<-CreateSeuratObject(counts = mpegl.data, project = "mpegles", min.cells = 3, min.features = 200, assay="RNA")
mpeg.combined <- merge(mpegn, y = mpegl, project = "MPEGCombined")


mito.genes <- grep(pattern = "^mt-", x = rownames(mpeg.combined@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(mpeg.combined@assays[["RNA"]][mito.genes, ])/Matrix::colSums(mpeg.combined@assays[["RNA"]])
mpeg.combined <- AddMetaData(object = mpeg.combined, metadata = percent.mito, col.name = "percent.mito") 
mpeg.comb<- subset(x = mpeg.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 0.2 )


mpeg.comb.list <- SplitObject(mpeg.comb)
mcl<-mpeg.comb.list

mcl <- lapply(X = mcl, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

mcl.anchors <- FindIntegrationAnchors(object.list = mcl, dims = 1:20)
mcl.combined <- IntegrateData(anchorset = mcl.anchors, dims = 1:20)


DefaultAssay(mcl.combined) <- "integrated"
mcl.combined <- ScaleData(mcl.combined, verbose = FALSE)
mcl.combined <- RunPCA(mcl.combined, npcs = 30, verbose = FALSE)
mcl.combined <- RunUMAP(mcl.combined, reduction = "pca", dims = 1:20)
mcl.combined <- FindNeighbors(mcl.combined, reduction = "pca", dims = 1:20)
mcl.combined <- FindClusters(mcl.combined, resolution = 0.3)
mcl.combined <- RunTSNE(mcl.combined, dims = 1:20)

Idents(mcl.combined) <- "seurat_clusters"

DimPlot(mcl.combined, reduction = "tsne", label=TRUE, pt.size=2.5, label.size = 8)


hern.data<-Read10X(data.dir="11609BC0001TXP-N__Her4_Na_ve/outs/filtered_feature_bc_matrix/") #replace with path to Naive dataset folder
herl.data<-Read10X(data.dir = "11609BC0002TXP-N__Her4_Lesioned/outs/filtered_feature_bc_matrix/") #replace with path to Lesioned dataset folder

hern<-CreateSeuratObject(counts = hern.data, project = "hernaive", min.cells = 3, min.features = 200, assay="RNA")
herl<-CreateSeuratObject(counts = herl.data, project = "herles", min.cells = 3, min.features = 200, assay="RNA")

her.combined <- merge(hern, y = herl, add.cell.ids = c("Naive", "Lesioned"), project = "Her4Combined")
mito.genes <- grep(pattern = "^mt-", x = rownames(her.combined@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(her.combined@assays[["RNA"]][mito.genes, ])/Matrix::colSums(her.combined@assays[["RNA"]])
her.combined <- AddMetaData(object = her.combined, metadata = percent.mito, col.name = "percent.mito") 

her.comb4<- subset(x = her.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mito >  -Inf & percent.mito < 0.05 )
her.comb.list4 <- SplitObject(her.comb4)
hcl4<-her.comb.list4

hcl4 <- lapply(X = hcl4, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


hcl4.anchors <- FindIntegrationAnchors(object.list = hcl4, dims = 1:20)
hcl4.combined <- IntegrateData(anchorset = hcl4.anchors, dims = 1:20)

DefaultAssay(hcl4.combined) <- "integrated"
hcl4.combined <- ScaleData(hcl4.combined, verbose = FALSE)
hcl4.combined <- RunPCA(hcl4.combined, npcs = 30, verbose = FALSE)
hcl4.combined <- RunUMAP(hcl4.combined, reduction = "pca", dims = 1:20)
hcl4.combined <- RunTSNE(hcl4.combine3d, dims = 1:20)
hcl4.combined <- FindNeighbors(hcl4.combined, reduction = "pca", dims = 1:20)
hcl4.combined <- FindClusters(hcl4.combined, resolution = 0.5)

DimPlot(hcl4.combined,reduction = "tsne", label=TRUE, pt.size = 2.5)



