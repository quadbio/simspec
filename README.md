# simspec: Similarity Spectrum
An R package to calculate representation of cells in single-cell genomic data, by their similarities to external references (RSS) or cell clusters in the data (CSS)

Installation
------
```
install.packages("devtools")
devtools::install_github("quadbiolab/simspec")
```

Usage
------
The more detailed vignette can be seen in https://github.com/quadbiolab/simspec/blob/master/vignette/vignette.md.

The codes to generate resulted reported in the manuscript (https://doi.org/10.1101/2020.02.27.968560) are deposited in https://github.com/quadbiolab/simspec/blob/master/code_repository/. Data can be retrieved from Mendeley Data (http://doi.org/10.17632/3kthhpw2pd).

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

