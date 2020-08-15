library(reticulate)                                                       
scanorama <- import("scanorama")
library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(simspec)
library(harmony)
library(liger)

library(RANN)
library(lisi)
library(doParallel)

# retrieve data
library(SeuratData)
InstallData("pbmcsca")
data("pbmcsca")
pbmcsca <- NormalizeData(pbmcsca) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

# integration
## CSS
library(simspec)
source("~/Tools/scripts/dimension_reduction.r")
pbmcsca <- cluster_sim_spectrum(pbmcsca, label_tag="Method", cluster_resolution = 0.4, corr_method = "pearson", spectrum_type = "corr_kernel", lambda = 50)
pbmcsca[['csspca']] <- CreateDimReducObject(trunc_PCA(Embeddings(pbmcsca, "css"), npcs = 10), key = "CSSPCA_")
pbmcsca <- cluster_sim_spectrum(pbmcsca, label_tag="Method", cluster_resolution = 0.4, corr_method = "spearman", reduction.key = "CSS2_", reduction.name = "css2")
pbmcsca <- RunUMAP(pbmcsca, reduction="csspca", dims = 1:ncol(Embeddings(pbmcsca, "csspca")), reduction.name="umap_css", reduction.key="UMAPCSS_")
pbmcsca <- RunUMAP(pbmcsca, reduction="css", dims = 1:ncol(Embeddings(pbmcsca, "css")), reduction.name = "umap_cssk", reduction.key="UMAPCSSK_")
pbmcsca <- RunUMAP(pbmcsca, reduction="css2", dims = 1:ncol(Embeddings(pbmcsca, "css2")), reduction.name = "umap_cssz", reduction.key="UMAPCSSZ_")
pbmcsca <- FindNeighbors(pbmcsca, reduction = "csspca", dims = 1:ncol(Embeddings(pbmcsca, "csspca")))
pbmcsca <- FindClusters(pbmcsca)
pbmcsca$cluster_css <- pbmcsca@active.ident

## harmony
library(harmony)
pbmcsca <- RunHarmony(pbmcsca, "Method", dims.use = 1:20, max.iter.harmony = 50)
pbmcsca <- RunUMAP(pbmcsca, reduction="harmony", dims = 1:ncol(Embeddings(pbmcsca, "harmony")), reduction.name="umap_harmony", reduction.key="UMAPHARMONY_")
pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:ncol(Embeddings(pbmcsca, "harmony")))
pbmcsca <- FindClusters(pbmcsca)
pbmcsca$cluster_harmony <- pbmcsca@active.ident

## fastMNN
pbmcsca_mnn <- RunFastMNN(object.list = SplitObject(pbmcsca, split.by = "Method"))
pbmcsca_mnn <- RunUMAP(pbmcsca_mnn, reduction="mnn", dims = 1:30)
pbmcsca_mnn <- FindNeighbors(pbmcsca_mnn, reduction = "mnn", dims = 1:30)
pbmcsca_mnn <- FindClusters(pbmcsca_mnn)
pbmcsca[['mnn']] <- CreateDimReducObject(Embeddings(pbmcsca_mnn,"mnn")[colnames(pbmcsca),], key="MNN_")
pbmcsca[['umap_mnn']] <- CreateDimReducObject(Embeddings(pbmcsca_mnn,"umap")[colnames(pbmcsca),], key="UMAPMNN_")
pbmcsca$cluster_mnn <- pbmcsca_mnn@active.ident[colnames(pbmcsca)]

## scanorama
scanorama <- import('scanorama')
mat <- setNames(lapply(SplitObject(pbmcsca, split.by = "Method"), function(x) as.matrix(t(x@assays$RNA@data[VariableFeatures(pbmcsca),]))), NULL)
gene_list <- lapply(unique(pbmcsca$Method), function(x) VariableFeatures(pbmcsca))
mat_integrated <- scanorama$integrate(mat, gene_list)
dr_scanorama <- do.call(rbind, mat_integrated[[1]])
rownames(dr_scanorama) <- do.call(c, lapply(mat, rownames))
dr_scanorama <- dr_scanorama[colnames(pbmcsca),]
pbmcsca[['scanorama']] <- CreateDimReducObject(dr_scanorama, key="SCANORAMA_")
pbmcsca <- RunUMAP(pbmcsca, reduction="scanorama", dims = 1:ncol(Embeddings(pbmcsca,"scanorama")), reduction.name="umap_scanorama", reduction.key="UMAPSCANORAMA_")
pbmcsca <- FindNeighbors(pbmcsca, reduction = "scanorama", dims = 1:ncol(Embeddings(pbmcsca, "scanorama")))
pbmcsca <- FindClusters(pbmcsca)
pbmcsca$cluster_scanorama <- pbmcsca@active.ident

## seurat
pbmcsca_samples <- lapply(SplitObject(pbmcsca, split.by = "Method"), function(obj) NormalizeData(obj) %>% FindVariableFeatures(nfeatures = 3000))
pbmcsca_anchors <- FindIntegrationAnchors(object.list = pbmcsca_samples, dims = 1:30, anchor.features = 3000)
pbmcsca_integrated <- IntegrateData(anchorset = pbmcsca_anchors, dims = 1:30)
pbmcsca_integrated <- ScaleData(pbmcsca_integrated) %>% RunPCA(npcs = 50, verbose = F) %>% RunUMAP(dims = 1:20)
pbmcsca_integrated <- FindNeighbors(pbmcsca_integrated, reduction = "pca", dims = 1:20)
pbmcsca_integrated <- FindClusters(pbmcsca_integrated)
pbmcsca[['umap_seurat']] <- CreateDimReducObject(Embeddings(pbmcsca_integrated,"umap")[colnames(pbmcsca),], key="UMAPSEURAT_")
pbmcsca[['pca_seurat']] <- CreateDimReducObject(Embeddings(pbmcsca_integrated, "pca")[colnames(pbmcsca),], key="PCASEURAT_")
pbmcsca$cluster_seurat <- pbmcsca_integrated@active.ident[colnames(pbmcsca)]

## liger
library(liger)
pbmcsca_liger <- ScaleData(pbmcsca, split.by = "Method", do.center = FALSE)
pbmcsca_liger <- RunOptimizeALS(pbmcsca_liger, k = 20, lambda = 5, split.by = "Method")
pbmcsca_liger <- RunQuantileAlignSNF(pbmcsca_liger, split.by = "Method")
pbmcsca_liger <- RunUMAP(pbmcsca_liger, dims = 1:ncol(pbmcsca_liger[["iNMF"]]), reduction = "iNMF")
pbmcsca_liger <- FindNeighbors(pbmcsca_liger, reduction = "iNMF", dims = 1:ncol(Embeddings(pbmcsca_liger, "iNMF"))) %>% FindClusters()
pbmcsca[['umap_liger']] <- CreateDimReducObject(Embeddings(pbmcsca_liger,"umap")[colnames(pbmcsca),], key="UMAPLIGER_")
pbmcsca[['liger']] <- CreateDimReducObject(Embeddings(pbmcsca_liger, "iNMF")[colnames(pbmcsca),], key="LIGER_")
pbmcsca$cluster_liger <- pbmcsca_liger@active.ident[colnames(pbmcsca)]



# kNN-based benchmark
k <- 50
knn_dr_pbmcsca <- setNames(lapply(c("pca","csspca","css","css2","mnn","scanorama","harmony","pca_seurat","liger"), function(dred){
  embed <- Embeddings(pbmcsca, dred)
  knn <- RANN::nn2(embed, k = k+1)$nn.idx
  return(knn)
}), c("pca","csspca","css","css2","mnn","scanorama","harmony","pca_seurat","liger"))

metrics_pbmcsca <- lapply(knn_dr_pbmcsca, function(knn){
  same_ct <- matrix(rep(pbmcsca$CellType[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) pbmcsca$CellType[ix])
  same_method <- matrix(rep(pbmcsca$Method[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) pbmcsca$Method[ix])
  res <- data.frame(rowMeans(same_ct & same_method),
                    rowMeans(same_ct & ! same_method),
                    rowMeans(!same_ct),
                    row.names = colnames(pbmcsca))
  return(res)
})

layout(matrix(1:2, nrow=1))
boxplot(sapply(metrics_pbmcsca, "[", 2),
        outline=F, frame=F, col = "#303030",
        names = names(metrics_pbmcsca), las=2, ylab="Proportion")
boxplot(sapply(metrics_pbmcsca, "[", 3),
        outline=F, frame=F, col = "#303030",
        names = names(metrics_pbmcsca), las=2, ylab="Proportion")


# LISI-based benchmark
registerDoParallel(9)
lisi_batch_pbmcsca <- foreach(dr = c("pca","csspca","css","css2","mnn","scanorama","harmony","pca_seurat","liger"), .combine = list, .multicombine = T) %dopar%{
  compute_lisi(Embeddings(pbmcsca, dr), pbmcsca@meta.data, "Method")
}
lisi_ct_pbmcsca <- foreach(dr = c("pca","csspca","css","css2","mnn","scanorama","harmony","pca_seurat","liger"), .combine = list, .multicombine = T) %dopar%{
  compute_lisi(Embeddings(pbmcsca, dr), pbmcsca@meta.data, "CellType")
}
stopImplicitCluster()

layout(matrix(1:2,nrow=1))
boxplot(do.call(cbind, lisi_batch_pbmcsca),
        frame=F, las=2, names = c("pca","csspca","css","css2","scanorama","mnn","harmony","pca_seurat","liger"), outline=F)
boxplot(length(unique(pbmcsca$CellType))+1-do.call(cbind, lisi_ct_pbmcsca),
        frame=F, las=2, names = c("pca","csspca","css","css2","scanorama","mnn","harmony","pca_seurat","liger"), outline=F)

