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
library(presto)

# import data and preprocess
counts <- readMM("cerebral_organoids_Kanton_2019_timecourse/matrix.mtx.gz")
meta <- read.table("cerebral_organoids_Kanton_2019_timecourse/cells.tsv.gz", sep="\t")
features <- read.table("cerebral_organoids_Kanton_2019_timecourse/features.tsv.gz", sep="\t", stringsAsFactors = F)
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- rownames(meta)
seurat <- CreateSeuratObject(counts = counts, meta.data = meta) %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)


# integration
## CSS
seurat <- cluster_sim_spectrum(seurat, label_tag="sample", cluster_resolution = 0.6)
seurat <- RunUMAP(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat,"css")), reduction.name="umap_css", reduction.key="UMACSS_")

## MNN
library(dplyr)
seurat_mnn <- CreateSeuratObject(seurat@assays$RNA@counts, meta.data = seurat@meta.data) %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 5000)
seurat_mnn <- RunFastMNN(object.list = SplitObject(seurat_mnn, split.by="sample"))
seurat_mnn <- RunUMAP(seurat_mnn, reduction="mnn", dims = 1:30)
seurat_mnn <- FindNeighbors(seurat_mnn, reduction = "mnn", dims = 1:30)
seurat_mnn <- FindClusters(seurat_mnn, resolution=0.3)
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"mnn")[colnames(seurat),], key="MNN_")
seurat[['umap_mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"umap")[colnames(seurat),], key="UMAPMNN_")

## scanorama
mat <- setNames(lapply(SplitObject(seurat, split.by = "sample"), function(x) as.matrix(t(x@assays$RNA@data[VariableFeatures(seurat),]))), NULL)
gene_list <- lapply(unique(seurat$sample), function(x) VariableFeatures(seurat))
mat_integrated <- scanorama$integrate(mat, gene_list)
dr_scanorama <- do.call(rbind, mat_integrated[[1]])
rownames(dr_scanorama) <- do.call(c, lapply(mat, rownames))
dr_scanorama <- dr_scanorama[colnames(seurat),]
seurat[['scanorama']] <- CreateDimReducObject(dr_scanorama, key="SCANORAMA_")
seurat <- RunUMAP(seurat, reduction="scanorama", dims = 1:ncol(Embeddings(seurat,"scanorama")), reduction.name="umap_scanorama", reduction.key="UMAPSCANORAMA_")

## Harmony
seurat <- RunHarmony(seurat, "sample")
seurat <- RunUMAP(seurat, reduction="harmony", reduction.name="umap_harmony", reduction.key="UMAPHARMONY_", dims = 1:50)

## Seurat
seurat_samples <- SplitObject(seurat, split.by="sample")
seurat_samples <- lapply(seurat_samples, function(obj) NormalizeData(obj) %>% FindVariableFeatures(nfeatures = 5000))
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_samples, dims = 1:30, anchor.features = 5000)
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)
seurat_integrated <- ScaleData(seurat_integrated) %>% RunPCA(npcs = 50, verbose = F) %>% RunUMAP(dims = 1:20)
seurat[['umap_Seurat']] <- CreateDimReducObject(Embeddings(seurat_integrated,"umap")[colnames(seurat),], key="UMAPSEURAT_")
seurat[['pca_Seurat']] <- CreateDimReducObject(Embeddings(seurat_integrated, "pca")[colnames(seurat),], key="PCASEURAT_")

## LIGER
data_samples <- lapply(seurat_samples, function(x) x@assays$RNA@counts)
liger_samples <- createLiger(data_samples)
liger_samples <- normalize(liger_samples)
liger_samples <- selectGenes(liger_samples, var.thresh=0.1, do.plot = T)
liger_samples <- scaleNotCenter(liger_samples)
liger_samples <- optimizeALS(liger_samples, k=20, thresh = 5e-5, nrep = 3)
liger_samples <- quantileAlignSNF(liger_samples, resolution = 0.4, small.clust.thresh = 20)
liger_samples <- runUMAP(liger_samples)
coord <- liger_samples@tsne.coords
coord <- coord[colnames(seurat),]
seurat[['umap_LIGER']] <- CreateDimReducObject(coord, key="UMAPLIGER_")
seurat[['LIGER']] <- CreateDimReducObject(liger_samples@H.norm[colnames(seurat),], key="LIGER_")



# kNN-based benchmark
library(RANN)
k <- 50
knn_dr <- setNames(lapply(c("pca","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"), function(dred){
  embed <- Embeddings(seurat, dred)
  knn <- RANN::nn2(embed, k = k+1)$nn.idx
  return(knn)
}), c("pca","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"))

## how well cells of different lines mixed (in relative to time point mixing)
stages <- as.character(seurat$stage_group)
stages[stages %in% c("iPSCs","EB")] <- "PSC"
line_mixing <- sapply(knn_dr, function(knn){
  same_line <- apply(knn[,-1], 2, function(x) seurat$line[x] == seurat$line[knn[,1]])
  same_stage <- apply(knn[,-1], 2, function(x) stages[x] == stages[knn[,1]])
  ratio_line <- (rowSums(!same_line & same_stage)+1)/(rowSums(same_line & same_stage)+1)
  log(ratio_line)
})
line_mixing_relative <- sapply(knn_dr, function(knn){
  same_line <- apply(knn[,-1], 2, function(x) seurat$line[x] == seurat$line[knn[,1]])
  same_stage <- apply(knn[,-1], 2, function(x) stages[x] == stages[knn[,1]])
  ratio_line <- (rowSums(!same_line & same_stage)+1)/(rowSums(same_line & same_stage)+1)
  ratio_fp <- (rowSums(same_line & !same_stage)+1)/(rowSums(same_line & same_stage)+1)
  log(ratio_line) - log(ratio_fp)
})

# how well cells of nearby time point mixed (in relative to other time points)
stages <- as.character(seurat$stage_group)
stages[stages %in% c("iPSCs","EB")] <- "PSC"
stages_ordered <- c("PSC","Neuroectoderm","Neuroepithelium","Organoid-1M","Organoid-2M","Organoid-4M")
registerDoParallel(20)
connectivity_stages <- sapply(knn_dr, function(knn){
  odds <- foreach(i = 1:nrow(knn), .combine = c) %dopar%{
    next_stage <- setNames(c(stages_ordered[-1],NA), stages_ordered)[stages[knn[i,1]]]
    prev_stage <- setNames(c(NA,stages_ordered[-length(stages_ordered)]), stages_ordered)[stages[knn[i,1]]]
    num_same_stage <- sum(stages[knn[i,-1]] == stages[knn[i,1]])+1
    num_nearby_stage <- sum(stages[knn[i,-1]] %in% c(next_stage, prev_stage))+1
    logodds <- -log((num_same_stage * sum(stages %in% c(next_stage, prev_stage))) / (num_nearby_stage * sum(stages == stages[knn[i,1]])))
    return(logodds)
  }
})
connectivity_stages_relative <- sapply(knn_dr, function(knn){
  odds <- foreach(i = 1:nrow(knn), .combine = c) %dopar%{
    next_stage <- setNames(c(stages_ordered[-1],NA), stages_ordered)[stages[knn[i,1]]]
    prev_stage <- setNames(c(NA,stages_ordered[-length(stages_ordered)]), stages_ordered)[stages[knn[i,1]]]
    num_same_stage <- sum(stages[knn[i,-1]] == stages[knn[i,1]])+1
    num_nearby_stage <- sum(stages[knn[i,-1]] %in% c(next_stage, prev_stage))+1
    num_others <- ncol(knn) - 1 - (num_same_stage-1) - (num_nearby_stage-1) + 1
    logodds <- -log((num_same_stage * sum(stages %in% c(next_stage, prev_stage))) / (num_nearby_stage * sum(stages == stages[knn[i,1]])))
    logodds_fp <- -log((num_same_stage * sum(!stages %in% c(next_stage, prev_stage, stages[knn[i,1]]))) / (num_others * sum(stages == stages[knn[i,1]])))
    return(logodds - logodds_fp)
  }
})
stopImplicitCluster()

layout(matrix(1:4,nrow=1))
boxplot(line_mixing, las=2, frame=F, outline=F)
boxplot(line_mixing_relative, las=2, frame=F, outline=F)
abline(h = 0, lty=2)
boxplot(connectivity_stages, las = 2, frame=F, outline=F)
boxplot(connectivity_stages_relative, las = 2, frame=F, outline=F)



# benchmark: intermediate cell states
## PSC and Necto
num_stage_nn_PSC_necto <- lapply(knn_dr, function(knn){
  idx <- which(seurat$stage_group %in% c("Neuroectoderm","iPSCs","EB"))
  t(apply(knn[idx, -1], 1, function(idx){
    stage_nn <- seurat$stage_group[idx]
    num_PSC <- sum(stage_nn %in% c("iPSCs","EB"))
    num_Necto <- sum(stage_nn %in% c("Neuroectoderm"))
    return(c(PSC = num_PSC, Necto = num_Necto))
  }))
})
idx_intermediate_PSC_necto <- lapply(num_stage_nn_PSC_necto, function(x)
  which(seurat$stage_group %in% c("Neuroectoderm","iPSCs","EB"))[which(x[,1]>15 & x[,2]>15)])
idx_intermediate_PSC_necto_union <- as.numeric(names(which(table(unlist(idx_intermediate_PSC_necto))>1)))
idx_nn_intermediate_PSC_necto <- do.call(c, lapply(knn_dr, function(knn) as.numeric(knn[idx_intermediate_PSC_necto_union,-1])))
idx_nn_intermediate_PSC_necto <- as.numeric(names(which(table(idx_nn_intermediate_PSC_necto)>10)))
idx_nn_intermediate_PSC_necto <- setdiff(idx_nn_intermediate_PSC_necto, idx_intermediate_PSC_necto_union)
idx_nn_intermediate_PSC <- intersect(idx_psc, idx_nn_intermediate_PSC_necto)
idx_nn_intermediate_necto <- intersect(idx_nect, idx_nn_intermediate_PSC_necto)

DE_PSC_necto <- wilcoxauc(X=seurat@assays$RNA@data[,c(idx_nn_intermediate_PSC, idx_nn_intermediate_necto)],
                          y=rep(c("PSC","necto"), c(length(idx_nn_intermediate_PSC), length(idx_nn_intermediate_necto))))
psc_high <- DE_PSC_necto$feature[DE_PSC_necto$group == "PSC" & DE_PSC_necto$auc > 0.7 & DE_PSC_necto$pct_in - DE_PSC_necto$pct_out > 30 & DE_PSC_necto$logFC>log(1.2) & DE_PSC_necto$padj<0.01]
necto_high <- DE_PSC_necto$feature[DE_PSC_necto$group == "necto" & DE_PSC_necto$auc > 0.7 & DE_PSC_necto$pct_in - DE_PSC_necto$pct_out > 30 & DE_PSC_necto$logFC>log(1.2) & DE_PSC_necto$padj<0.01]
scores_psc_necto <- colMeans(seurat@assays$RNA@data[psc_high,]) - colMeans(seurat@assays$RNA@data[necto_high,])

boxplot(c(list(scores_psc_necto[idx_nn_intermediate_PSC]),
          lapply(idx_intermediate_PSC_necto, function(x) scores_psc_necto[x]),
          list(scores_psc_necto[idx_nn_intermediate_necto])),
        frame = F, outline=F, names=c("NEpith",names(idx_intermediate_PSC_necto),"Org1m"), las = 2, ylab = "NEpith-Org1m score")
points(x = rep(1, length(idx_nn_intermediate_PSC))+runif(length(idx_nn_intermediate_PSC),-0.25,0.25), scores_psc_necto[idx_nn_intermediate_PSC], pch=21, cex=0.8, col="#303030", bg="#31a354", lwd=0.3)
for(i in 1:length(idx_intermediate_PSC_necto))
  points(x = rep(i+1, length(idx_intermediate_PSC_necto[[i]]))+runif(length(idx_intermediate_PSC_necto[[i]]),-0.25,0.25), scores_psc_necto[idx_intermediate_PSC_necto[[i]]], pch=21, cex=0.8, col="#303030", bg="#98BB2A", lwd=0.3)
points(x = rep(length(idx_intermediate_PSC_necto)+2, length(idx_nn_intermediate_necto))+runif(length(idx_nn_intermediate_necto),-0.25,0.25), scores_psc_necto[idx_nn_intermediate_necto], pch=21, cex=0.8, col="#303030", bg="#FFD400", lwd=0.3)

## Nepith and Org-1M 
num_stage_nn_Nepith_org1m <- lapply(knn_dr, function(knn){
  idx <- which(seurat$stage_group %in% c("Neuroepithelium","Organoid-1M"))
  t(apply(knn[idx, -1], 1, function(idx){
    stage_nn <- seurat$stage_group[idx]
    num_Nepith <- sum(stage_nn %in% c("Neuroepithelium"))
    num_Org1m <- sum(stage_nn %in% c("Organoid-1M"))
    return(c(Nepith = num_Nepith, Org1m = num_Org1m))
  }))
})
idx_intermediate <- lapply(num_stage_nn_Nepith_org1m, function(x)
  which(seurat$stage_group %in% c("Neuroepithelium","Organoid-1M"))[which(x[,1]>15 & x[,2]>15)])
idx_intermediate_union <- as.numeric(names(which(table(unlist(idx_intermediate))>1)))
idx_nn_intermediate <- do.call(c, lapply(knn_dr, function(knn) as.numeric(knn[idx_intermediate_union,-1])))
idx_nn_intermediate <- as.numeric(names(which(table(idx_nn_intermediate)>10)))
idx_nn_intermediate <- setdiff(idx_nn_intermediate, idx_intermediate_union)
idx_nn_intermediate_nepi <- intersect(idx_nepi, idx_nn_intermediate)
idx_nn_intermediate_org1m <- intersect(idx_org1m, idx_nn_intermediate)

DE_nepi_org1m <- wilcoxauc(X=seurat@assays$RNA@data[,c(idx_nn_intermediate_nepi, idx_nn_intermediate_org1m)],
                           y=rep(c("nepi","org1m"), c(length(idx_nn_intermediate_nepi), length(idx_nn_intermediate_org1m))))
nepi_high <- DE_nepi_org1m$feature[DE_nepi_org1m$group == "nepi" & DE_nepi_org1m$auc > 0.7 & DE_nepi_org1m$pct_in - DE_nepi_org1m$pct_out > 30 & DE_nepi_org1m$logFC>log(1.2) & DE_nepi_org1m$padj<0.01]
org1m_high <- DE_nepi_org1m$feature[DE_nepi_org1m$group == "org1m" & DE_nepi_org1m$auc > 0.7 & DE_nepi_org1m$pct_in - DE_nepi_org1m$pct_out > 30 & DE_nepi_org1m$logFC>log(1.2) & DE_nepi_org1m$padj<0.01]
scores <- colMeans(seurat@assays$RNA@data[nepi_high,]) - colMeans(seurat@assays$RNA@data[org1m_high,])

boxplot(c(list(scores[idx_nn_intermediate_nepi]),
          lapply(idx_intermediate, function(x) scores[x]),
          list(scores[idx_nn_intermediate_org1m])),
        frame = F, outline=F, names=c("NEpith",names(idx_intermediate),"Org1m"), las = 2, ylab = "NEpith-Org1m score", ylim=c(-1.1,1))
points(x = rep(1, length(idx_nn_intermediate_nepi))+runif(length(idx_nn_intermediate_nepi),-0.25,0.25), scores[idx_nn_intermediate_nepi], pch=21, cex=0.8, col="#303030", bg="#31a354", lwd=0.3)
for(i in 1:length(idx_intermediate))
  points(x = rep(i+1, length(idx_intermediate[[i]]))+runif(length(idx_intermediate[[i]]),-0.25,0.25), scores[idx_intermediate[[i]]], pch=21, cex=0.8, col="#303030", bg="#98BB2A", lwd=0.3)
points(x = rep(length(idx_intermediate)+2, length(idx_nn_intermediate_org1m))+runif(length(idx_nn_intermediate_org1m),-0.25,0.25), scores[idx_nn_intermediate_org1m], pch=21, cex=0.8, col="#303030", bg="#FFD400", lwd=0.3)



# DE between PSC-Necto intermediate cells (based on CSS) and PSC/neuroectoderm cells
num_stage_nn_PSC_NEcto <- t(apply(knn_dr$CSS[which(seurat$stage_group %in% c("iPSCs","EB","Neuroectoderm")), -1], 1, function(idx){
  stage_nn <- seurat$stage_group[idx]
  num_PSC <- sum(stage_nn %in% c("iPSCs","EB"))
  num_NEcto <- sum(stage_nn %in% c("Neuroectoderm"))
  return(c(PSC = num_PSC, NEcto = num_NEcto))
}))
idx <- which(seurat$stage_group %in% c("iPSCs","EB","Neuroectoderm"))[num_stage_nn_PSC_NEcto[,2] > 15 & num_stage_nn_PSC_NEcto[,1] > 15]
idx_PSC <- which(seurat$stage_group %in% c("iPSCs","EB","Neuroectoderm"))[num_stage_nn_PSC_NEcto[,1] > 45]
idx_NE <- which(seurat$stage_group %in% c("iPSCs","EB","Neuroectoderm"))[num_stage_nn_PSC_NEcto[,2] > 45]

DE_intermediate_PSC <- wilcoxauc(seurat@assays$RNA@data[,c(idx,idx_PSC)], rep(c("intermediate","PSC"), c(length(idx),length(idx_PSC))))
DE_intermediate_NE <- wilcoxauc(seurat@assays$RNA@data[,c(idx,idx_NE)], rep(c("intermediate","NEcto"), c(length(idx),length(idx_NE))))
