# simspec: Similarity Spectrum
An R package to calculate representation of cells in single-cell genomic data, by their similarities to external references (RSS) or cell clusters in the data (CSS). More details of the method are available in the paper **[CSS: cluster similarity spectrum integration of single-cell genomics data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02147-4)**. The manscript is also available in [biorxiv](https://doi.org/10.1101/2020.02.27.968560).

Recent update
------
(240226)
1. Add support to Assay5 in Seurat v5
2. Remove the qlcMatrix dependency

(221101)
1. Implement `estimate_projection_failure` function to estimate failure likelihood of data projection to the given reference for each query cell
2. Add verbose messages to the `transfer_labels` function
3. Support providing cluster labels instead of doing clustering per sample from scratch in `cluster_sim_spectrum`
4. Update verbose message

(220622)
1. Add `min_cluster_num` parameter to `cluster_sim_spectrum` function to exclude samples with too few clusters from the ref profiles
2. Support `ref_sim_spectrum` function to output as a new assay in the Seurat object
3. Update verbose message

(211124)
1. Sparse matrix ranking for Spearman correlation coefficient to speed up calculation and avoid conversion to dense matrix
2. Faster kNN-based label projection

Installation
------
```
install.packages("devtools")
devtools::install_github("quadbiolab/simspec")
```

Usage
------
The more detailed vignette can be seen in https://github.com/quadbiolab/simspec/blob/master/vignette/vignette.md.

The codes to generate resulted reported in the paper are deposited in https://github.com/quadbiolab/simspec/blob/master/code_repository/. Data can be retrieved from Mendeley Data (http://doi.org/10.17632/3kthhpw2pd).

### Reference Similarity Spectrum (RSS)
To calculate RSS, two inputs are required
* Expression matrix of the data (expr)
* Expression matrix of the reference (ref)
```
RSS <- ref_sim_spectrum(expr, ref)
```
A Seurat object can also be the input. In that case, an updated Seurat object with additional dimension reduction ('rss' by default) is returned
```
seurat <- ref_sim_spectrum(seurat, ref)
seurat <- RunUMAP(seurat, reduction = "rss", dims = 1:ncol(Embeddings(seurat, "rss")))
seurat <- FindNeighbors(seurat, reduction = "rss", dims = 1:ncol(Embeddings(seurat, "rss")))
seurat <- FindClusters(seurat)
UMAPPlot(seurat)
```

### Cluster Similarity Spectrum (CSS)
To calculate CSS, two inputs are required
* Expression matrix of the data (expr)
* Labels indicating samples (labels)
```
CSS <- cluster_sim_spectrum(expr, labels = labels)
```
Similarly, a Seurat object can be the input. When a Seurat object is used, the name of a column in the meta.data, which shows labels of samples, should be provided.
*Note: the Seurat object is expected to have variable features defined and PCA run*
```
seurat <- cluster_sim_spectrum(seurat, label_tag = "sample")
seurat <- RunUMAP(seurat, reduction = "css", dims = 1:ncol(Embeddings(seurat, "css")))
seurat <- FindNeighbors(seurat, reduction = "css", dims = 1:ncol(Embeddings(seurat, "css")))
seurat <- FindClusters(seurat)
UMAPPlot(seurat)
```
### CSS for query data projection
CSS representation allows simple and straightforward projection of query data to a reference atlas. To do that, the CSS representation model of the reference data needs to be returned.
```
model <- cluster_sim_spectrum(expr_ref, labels = labels_ref, return_css_only = F)
model <- cluster_sim_spectrum(seurat_ref, label_tag = "sample", return_seuratObj = F)
```
The model is then used to project query data to the same CSS space
```
css_query <- css_project(expr_query, model)
seurat_query <- css_project(seurat_query, model)
```

