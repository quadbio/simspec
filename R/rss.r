#'@param ref Expression matrix of the reference data set
#'@param method Correlation method
#'@param scale If TRUE, z-transform is applied to the calculated similarities across reference samples
#'@rdname ref_sim_spectrum
#'@export
#'@method ref_sim_spectrum default
ref_sim_spectrum.default <- function(object, ref, method = c("pearson","spearman"), scale = TRUE){
  method <- match.arg(method)
  
  candidates <- intersect(rownames(object), rownames(ref))
  corr <- cor(as.matrix(object[candidates,]), ref[candidates,], method = method)
  if (scale)
    corr <- t(scale(t(corr)))
  return(corr)
}

#'@param reduction.name Reduction name of thea RSS representation in the returned Seurat object
#'@param reduction.key Reduction key of the RSS representation in the returned Seurat object
#'@rdname ref_sim_spectrum
#'@export
#'@method ref_sim_spectrum Seurat
ref_sim_spectrum.Seurat <- function(object, ref, ..., reduction.name = "rss", reduction.key = "RSS_"){
  input <- object@assays[[DefaultAssay(object)]]@data
  rss <- ref_sim_spectrum.default(input, ref, ...)
  object[[reduction.name]] <- CreateDimReducObject(rss, key = reduction.key, assay = DefaultAssay(object))
  return(object)
}


