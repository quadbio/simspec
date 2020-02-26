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
### Reference Similarity Spectrum (RSS)
To calculate RSS, two inputs are required
..* Expression matrix of the data (expr)
..* Expression matrix of the reference (ref)
```
RSS <- ref_sim_spectrum(expr, ref)
```
A Seurat object can also be the input. In that case, an updated Seurat object with additional dimension reduction ('rss' by default) is returned
```
seurat <- ref_sim_spectrum(seurat, ref)
seurat <- RunUMAP(seurat, reduction = "rss", dims = 1:ncol(Embeddings(seurat, "rss")))
```

### Cluster Similarity Spectrum (CSS)
To calculate CSS, two inputs are required
..* Expression matrix of the data (expr)
..* Labels indicating samples (labels)
```
CSS <- cluster_sim_spectrum(expr, labels = labels)
```
Similarly, a Seurat object can be the input. When a Seurat object is used, the name of a column in the meta.data, which shows labels of samples, should be provided
```
seurat <- cluster_sim_spectrum(seurat, label_tag = "sample")
seurat <- RunUMAP(seurat, reduction = "css", dims = 1:ncol(Embeddings(seurat, "css")))
```
