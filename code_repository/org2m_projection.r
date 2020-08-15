library(Seurat)
library(simspec)
library(dplyr)
source("get_umap_model.r")

library(RANN)

# load data and split ref and query
seurat <- readRDS("cerebral_org2m_processed.rds")

query_samples <- unique(seurat$organoid[seurat$line == "SC102A1"])
ref_samples <- setdiff(unique(seurat$organoid), query_samples)
seurat_ref <- subset(seurat, subset = organoid %in% ref_samples)
seurat_query <- subset(seurat, subset = organoid %in% query_samples)

# CSS model for ref
seurat_ref <- FindVariableFeatures(seurat_ref, nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = F) %>%
  RunUMAP(dims = 1:20)
css_model_ref <- cluster_sim_spectrum(seurat_ref,
                                      label_tag = "organoid",
                                      cluster_resolution = 0.6,
                                      return_seuratObj = F)
seurat_ref[['CSS']] <- CreateDimReducObject(css_model_ref$sim2profiles, key="CSS_")

umap_model_ref <- get_umap_model(seurat_ref, reduction = "CSS", dims = 1:ncol(Embeddings(seurat_ref, "CSS")))
rownames(umap_model_ref$embedding) <- colnames(seurat_ref)
seurat_ref[['umap_CSS']] <- CreateDimReducObject(umap_model_ref$embedding, key="UMAPCSS_")

# projection for query
## CSS representation
seurat_query <- css_project(seurat_query, css_model_ref, reduction.name = "css_proj")

## UMAP
proj_umap_query <- umap_transform(Embeddings(seurat_query, "css_proj"), umap_model_ref)
rownames(proj_umap_query) <- colnames(seurat_query)
seurat_query[['umap_proj']] <- CreateDimReducObject(proj_umap_query, key="UMAPPROJ_")

## cell type
k <- 50
nn_query_ref <- RANN::nn2(Embeddings(seurat_ref, "CSS"), Embeddings(seurat_query, "css_proj"), k = k)
proj_ct_query <- apply(nn_query_ref$nn.idx, 1, function(x)
  names(which.max(table(seurat_ref$celltype_reannot[x]))))
seurat_query$proj_cl <- proj_cl_query

