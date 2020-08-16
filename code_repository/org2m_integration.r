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

# retrieve data and preprocess

## Important: the actually analysis was done on the preprocessed data 
## described in [Kanton et al. 2019], which was relied on Seurat v2.
## Therefore, using the following preprocessing would work but return 
## slightly different results (esp. for the cluster labels used to build
## cell-level regional identity predictor)

#counts <- readMM("cerebral_organoids_Kanton_2019_org2m/matrix.mtx.gz")
#meta <- read.table("cerebral_organoids_Kanton_2019_org2m/cells.tsv.gz", sep="\t")
#features <- read.table("cerebral_organoids_Kanton_2019_org2m/features.tsv.gz", sep="\t", stringsAsFactors = F)
#rownames(counts) <- make.unique(features[,2])
#colnames(counts) <- rownames(meta)
#seurat <- CreateSeuratObject(counts = counts, meta.data = meta) %>%
#  NormalizeData() %>%
#  FindVariableFeatures(nfeatures = 5000) %>%
#  ScaleData()

seurat <- readRDS("cerebral_org2m_processed.rds")
seurat <- RunPCA(seurat, npcs = 20, verbose=F) %>%
  RunUMAP(seurat, reduction = "pca", dims = 1:20)


# integration
## RSS
ref_brainspan <- readRDS("ext/brainspan_fetal.rds")
seurat <- ref_sim_spectrum(seurat, ref_brainspan, reduction.name = "RSS", reduction.key = "RSS_")
seurat <- RunUMAP(seurat, reduction="RSS", reduction.name="umap_RSS", reduction.key="UMAPRSS_", dims = 1:ncol(Embeddings(seurat, "RSS")))

## CSS
seurat <- cluster_sim_spectrum(seurat, label_tag="organoid", cluster_resolution = 0.6)
seurat <- RunUMAP(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat,"css")), reduction.name="umap_css", reduction.key="UMACSS_")

## MNN
seurat_mnn <- RunFastMNN(object.list = SplitObject(seurat, split.by = "organoid"))
seurat_mnn <- RunUMAP(seurat_mnn, reduction="mnn", dims = 1:30)
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"mnn")[colnames(seurat),], key="MNN_")
seurat[['umap_mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"umap")[colnames(seurat),], key="UMAPMNN_")

## scanorama
mat <- setNames(lapply(SplitObject(seurat, split.by = "organoid"), function(x) as.matrix(t(x@assays$RNA@data[VariableFeatures(seurat),]))), NULL)
gene_list <- lapply(unique(seurat$organoid), function(x) VariableFeatures(seurat))
mat_integrated <- scanorama$integrate(mat, gene_list)
dr_scanorama <- do.call(rbind, mat_integrated[[1]])
rownames(dr_scanorama) <- do.call(c, lapply(mat, rownames))
dr_scanorama <- dr_scanorama[colnames(seurat),]
seurat[['scanorama']] <- CreateDimReducObject(dr_scanorama, key="SCANORAMA_")
seurat <- RunUMAP(seurat, reduction="scanorama", dims = 1:ncol(Embeddings(seurat,"scanorama")), reduction.name="umap_scanorama", reduction.key="UMAPSCANORAMA_")

## Harmony
seurat <- RunHarmony(seurat, "organoid")
seurat <- RunUMAP(seurat, reduction="harmony", reduction.name="umap_harmony", reduction.key="UMAPHARMONY_", dims = 1:50)

## Seurat integration
seurat_samples <- SplitObject(seurat, split.by="organoid")
seurat_samples <- lapply(seurat_samples, function(obj) FindVariableFeatures(obj, nfeatures = 5000))
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_samples, dims = 1:30, anchor.features = 5000)
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)
seurat_integrated <- ScaleData(seurat_integrated) %>% RunPCA(npcs = 50, verbose = F) %>% RunUMAP(dims = 1:20)
seurat[['umap_Seurat']] <- CreateDimReducObject(Embeddings(seurat_integrated, "umap_integrated")[colnames(seurat),], key="UMAPSEURAT_")
seurat[['pca_Seurat']] <- CreateDimReducObject(Embeddings(seurat_integrated, "pca")[colnames(seurat),], key = "PCASEURAT_")

## LIGER
seurat_samples <- SplitObject(seurat, split.by="organoid")
data_samples <- lapply(seurat_samples, function(x) x@assays$RNA@data)
data_samples <- lapply(data_samples, function(x) exp(x)-1)
liger_samples <- createLiger(data_samples)
liger_samples <- normalize(liger_samples)
liger_samples <- selectGenes(liger_samples, var.thresh=0.3, do.plot = T)
liger_samples <- scaleNotCenter(liger_samples)
liger_samples <- optimizeALS(liger_samples, k=20, thresh = 5e-5, nrep = 3)
liger_samples <- runUMAP(liger_samples, use.raw=T)
liger_samples <- quantileAlignSNF(liger_samples, resolution = 0.4, small.clust.thresh = 20)
liger_samples <- runUMAP(liger_samples)
coord <- liger_samples@tsne.coords
coord <- coord[colnames(seurat),]
seurat[['umap_LIGER']] <- CreateDimReducObject(coord, key="UMAPLIGER_")
seurat[['LIGER']] <- CreateDimReducObject(liger_samples@H.norm[colnames(seurat),], key="LIGER_")



# cell-level annotation
### dorsal-ventral-nontelen
library(glmnet)
library(doMC)
registerDoMC(10)
seurat <- FindNeighbors(seurat, reduction="pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6)

cl_dorsal <- c(2,18,0,12,3,11)
cl_ventral <- c(4,16,7,13,8)
cl_nontelen <- c(19,10,22,20,15)
idx <- c(sample(which(seurat@active.ident %in% cl_dorsal), 1000),
         sample(which(seurat@active.ident %in% cl_ventral), 1000),
         sample(which(seurat@active.ident %in% cl_nontelen), 1000))

x <- as.matrix(seurat@assays$RNA@data[VariableFeatures(seurat),idx])
y <- factor(setNames(rep(c("dorsal","ventral","nontelen"), c(6,5,5)), c(cl_dorsal,cl_ventral,cl_nontelen))[as.character(seurat@active.ident[idx])])
m <- cv.glmnet(x = t(x), y = y, family = "multinomial", parallel = T)
pred_fates <- predict(m, as.matrix(t(seurat@assays$RNA@data[VariableFeatures(seurat),])), type = "response")[,,1]
lab_pred_fates <- factor(colnames(pred_fates)[apply(pred_fates, 1, which.max)])

### NPC-neuron
DE_neuron_NPC <- read.table("ext/DE_NPC_neurons.tsv", sep="\t", stringsAsFactors=F)
high_neuron <- intersect(rownames(seurat), DE_neuron_NPC$feature[which(DE_neuron_NPC$group == "neuron" & DE_neuron_NPC$padj < 0.01 & DE_neuron_NPC$logFC > log(1.2) & DE_neuron_NPC$pct_in > 50 & DE_neuron_NPC$pct_out < 20 & DE_neuron_NPC$auc > 0.6)])
high_NPC <- intersect(rownames(seurat), DE_neuron_NPC$feature[which(DE_neuron_NPC$group == "NPC" & DE_neuron_NPC$padj < 0.01 & DE_neuron_NPC$logFC > log(1.2) & DE_neuron_NPC$pct_in > 50 & DE_neuron_NPC$pct_out < 20 & DE_neuron_NPC$auc > 0.6)])
seurat$neuron_vs_NPC <- colMeans(seurat@assays$RNA@data[high_neuron,]) - colMeans(seurat@assays$RNA@data[high_NPC,])
lab_NPC_neuron <- factor(ifelse(seurat$neuron_vs_NPC > 0, "neuron", "NPC"))

lab_annot <- paste0(lab_pred_fates, "_", lab_NPC_neuron)
seurat$semi_branch <- lab_pred_fates
seurat$semi_celltype <- lab_annot



# kNN-based benchmark
k <- 50
knn_dr_org2m <- setNames(lapply(c("pca","RSS","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"), function(dred){
  embed <- Embeddings(seurat, dred)
  knn <- RANN::nn2(embed, k = k+1)$nn.idx
  return(knn)
}), c("pca","RSS","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"))

metrics_org2m <- lapply(knn_dr_pbmcsca, function(knn){
  same_ct <- matrix(rep(seurat$celltype_reannot[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) seurat$celltype_reannot[ix])
  same_fates <- matrix(rep(seurat$lineage_RSS[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) seurat$lineage_RSS[ix])
  same_org <- matrix(rep(seurat$organoid[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) seurat$organoid[ix])
  same_batch <- matrix(rep(seurat$batch[knn[,1]], ncol(knn)-1), ncol = ncol(knn)-1) == apply(knn[,-1], 2, function(ix) seurat$batch[ix])
  
  res <- data.frame(rowMeans(same_ct & ! same_batch),
                    rowMeans(!same_ct & same_fates),
                    rowMeans(!same_fates),
                    row.names = colnames(seurat))
  return(res)
})

layout(matrix(1:3, nrow=1))
boxplot(sapply(metrics_org2m, "[", 1),
        outline=F, frame=F, col = "#303030", names = names(metrics_org2m), las=2, ylab="Proportion", cex=0.5, pch=16)
boxplot(sapply(metrics_org2m, "[", 2),
        outline=F, frame=F, col = "#bdbdbd", names = names(metrics_org2m), las=2, ylab="Proportion", cex=0.5, pch=16)
boxplot(sapply(metrics_org2m, "[", 3),
        outline=F, frame=F, col = "#efefef", names = names(metrics_org2m), las=2, ylab="Proportion", cex=0.5, pch=16)


# LISI-based benchmark
registerDoParallel(8)
lisi_batch_org2m <- foreach(dr = c("pca","RSS","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"), .combine = list, .multicombine = T) %dopar%{
  compute_lisi(Embeddings(seurat, dr), seurat@meta.data, "organoid")
}
lisi_ct_org2m <- foreach(dr = c("pca","RSS","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"), .combine = list, .multicombine = T) %dopar%{
  compute_lisi(Embeddings(seurat, dr), seurat@meta.data, "celltype_reannot")
}
stopImplicitCluster()

layout(matrix(1:2, nrow=1))
boxplot(do.call(cbind, lisi_batch_org2m), frame=F, las=2, names = c("pca","RSS","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"), outline=F)
boxplot(length(unique(seurat$celltype_reannot))+1-do.call(cbind, lisi_ct_org2m), frame=F, las=2, names = c("pca","RSS","CSS","mnn","scanorama","harmony","pca_Seurat","LIGER"), outline=F)
