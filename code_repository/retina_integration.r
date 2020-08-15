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
library(doParallel)

# create Seurat object
counts <- readMM("GSE138002_Final_matrix.mtx.gz")
meta <- read.csv("GSE138002_Final_barcodes.csv.gz", sep=";", row.names=1)
genes <- read.csv("GSE138002_genes.csv.gz", sep=";", row.names=1)
rownames(counts) <- make.unique(genes$gene_short_name)
colnames(counts) <- rownames(meta)
seurat <- CreateSeuratObject(counts = counts, meta.data = meta)

# preprocess without integration
seurat <- subset(seurat, subset = sample != "24_Day") # D24 sample contained too few cells
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)

# integration
## CSS
seurat <- cluster_sim_spectrum(seurat, label_tag="sample")
seurat <- run_PCA(seurat, reduction="css", npcs = 20)
seurat <- RunUMAP(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat,"css")), reduction.name="umap_css", reduction.key="UMACSS_")

## MNN
seurat_mnn <- RunFastMNN(object.list = SplitObject(seurat, split.by = "sample"))
seurat_mnn <- RunUMAP(seurat_mnn, reduction="mnn", dims = 1:30)
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"mnn")[colnames(seurat),], key="MNN_")
seurat[['umap_mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"umap")[colnames(seurat),], key="UMAPMNN_")

## Scanorama
mat <- setNames(lapply(SplitObject(seurat, split.by = "sample"), function(x) as.matrix(t(x@assays$RNA@data[VariableFeatures(seurat),]))), NULL)
gene_list <- lapply(unique(seurat$sample), function(x) VariableFeatures(seurat))
mat_integrated <- scanorama$integrate(mat, gene_list)
dr_scanorama <- do.call(rbind, mat_integrated[[1]])
rownames(dr_scanorama) <- do.call(c, lapply(mat, rownames))
dr_scanorama <- dr_scanorama[colnames(seurat),]
seurat[['scanorama']] <- CreateDimReducObject(dr_scanorama, key="SCANORAMA_")
seurat <- RunUMAP(seurat, reduction="scanorama", dims = 1:ncol(Embeddings(seurat,"scanorama")), reduction.name="umap_scanorama", reduction.key="UMAPSCANORAMA_")

## Harmony
seurat <- RunHarmony(seurat, group.by.vars="sample", max.iter.harmony = 50)
seurat <- RunUMAP(seurat, reduction="harmony", dims = 1:50, reduction.name="umap_harmony", reduction.key="UMAPHARMONY_")

## Seurat v3
seurat_samples <- lapply(SplitObject(seurat, split.by = "sample"), function(obj) NormalizeData(obj) %>% FindVariableFeatures(nfeatures = 3000))
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_samples[c(2:24)], dims = 1:30, anchor.features = 3000)
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)
seurat_integrated <- ScaleData(seurat_integrated) %>% RunPCA(npcs = 50, verbose = F) %>% RunUMAP(dims = 1:20)
seurat[['umap_seurat']] <- CreateDimReducObject(Embeddings(seurat_integrated,"umap")[colnames(seurat),], key="UMAPSEURAT_")
seurat[['pca_seurat']] <- CreateDimReducObject(Embeddings(seurat_integrated, "pca")[colnames(seurat),], key="PCASEURAT_")

## LIGER
seurat_liger <- ScaleData(seurat, split.by = "sample", do.center = FALSE)
seurat_liger <- RunOptimizeALS(seurat_liger, k = 20, lambda = 5, split.by = "sample")
seurat_liger <- RunQuantileAlignSNF(seurat_liger, split.by = "sample")
seurat_liger <- RunUMAP(seurat_liger, dims = 1:ncol(seurat_liger[["iNMF"]]), reduction = "iNMF")
seurat[['umap_liger']] <- CreateDimReducObject(Embeddings(seurat_liger,"umap")[colnames(seurat),], key="UMAPLIGER_")
seurat[['liger']] <- CreateDimReducObject(Embeddings(seurat_liger, "iNMF")[colnames(seurat),], key="LIGER_")



# kNN-based benchmark (cell type vs. time points)
k <- 50
registerDoParallel(8)
knn_dr <- setNames(foreach(dred = c("pca","css","scanorama","mnn","harmony","pca_seurat","liger"), .combine = list, .multicombine=T) %dopar%{
  embed <- Embeddings(seurat, dred)
  knn <- RANN::nn2(embed, k = k+1)$nn.idx
  return(knn)
}, c("pca","css","scanorama","mnn","harmony","pca_seurat","liger"))
stopImplicitCluster()

age_num <- seurat$age
age_num[seurat$sample_type=="Retinal Organoid"] <- as.numeric(gsub("_Day","",age_num[seurat$sample_type=="Retinal Organoid"]))
age_num[grep("Hgw",age_num)] <- as.numeric(gsub("Hgw","",age_num[grep("Hgw",age_num)]))*7
age_num[grep("Hpnd",age_num)] <- as.numeric(gsub("Hpnd","",age_num[grep("Hpnd",age_num)]))+280
age_num[grep("Adult",age_num)] <- 365*20
age_num <- as.numeric(age_num)
age_rank <- setNames(rank(sort(unique(age_num))), sort(unique(age_num)))[as.character(age_num)]

prop_knn_ct_tp <- lapply(knn_dr, function(knn){
  same_ct <- matrix(rep(seurat$umap2_CellType[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) seurat$umap2_CellType[ix])
  same_tp <- matrix(rep(seurat$age[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) seurat$age[ix])
  prev_tp <- matrix(rep(age_rank[knn[,1]]-1, ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) age_rank[ix]) 
  next_tp <- matrix(rep(age_rank[knn[,1]]+1, ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) age_rank[ix])
  
  prop_same_ct <- rowMeans(same_ct)
  prop_nearby_tp <- rowMeans((prev_tp | next_tp) & !same_tp)
  prop_faraway_tp <- rowMeans(!(same_tp | prev_tp | next_tp))
  res <- data.frame(prop_same_ct,
                    prop_nearby_tp,
                    prop_faraway_tp)
  return(res)
})

layout(matrix(c(1:3), nrow=1))
boxplot(sapply(prop_knn_ct_tp, "[", 2),
        frame = F, names = names(knn_dr), las = 2, outline = F)
boxplot(sapply(prop_knn_ct_tp, function(x) (x[,2])/(x[,2]+x[,3])),
        frame = F, names = names(knn_dr), las = 2, outline = F)
boxplot(1-sapply(prop_knn_ct_tp, "[", 1),
        frame = F, names = names(knn_dr), las = 2, outline = F)


# kNN-based benchmark (intermediate)
registerDoParallel(8)
prop_knn_ac_hc_pre <- foreach(knn = knn_dr, .combine = list, .multicombine = T) %dopar%{
  idx <- which(seurat$umap2_CellType == "AC/HC_Precurs")
  return(t(apply(knn[idx,-1], 1, function(ix) table(seurat$umap2_CellType[ix])/(ncol(knn)-1))))
}
prop_knn_bip_pr_pre <- foreach(knn = knn_dr, .combine = list, .multicombine = T) %dopar%{
  idx <- which(seurat$umap2_CellType == "BC/Photo_Precurs")
  return(t(apply(knn[idx,-1], 1, function(ix) table(seurat$umap2_CellType[ix])/(ncol(knn)-1))))
}
stopImplicitCluster()

barplot(sapply(prop_knn_ac_hc_pre[-3], function(x) colMeans(x))[c(8,1,2,6,3:5,7,9:11),], border=NA, col=rep(c("#000000","#909090","#bdbdbd","#efefef"), c(1,1,2,7)))
barplot(sapply(prop_knn_bip_pr_pre[-3], function(x) colMeans(x))[c(8,3,4,5,10,1:2,6:7,9,11),], border=NA, col=rep(c("#000000","#909090","#bdbdbd","#efefef"), c(1,1,3,6)))
