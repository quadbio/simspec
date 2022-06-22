#'@param ref Expression matrix of the reference data set
#'@param method Correlation method
#'@param use_fast_rank When the presto package is available, use its rank_matrix function to rank sparse matrix
#'@param scale If TRUE, z-transform is applied to the calculated similarities across reference samples
#'@rdname ref_sim_spectrum
#'@method ref_sim_spectrum default
ref_sim_spectrum.default <- function(object, ref, method = c("pearson","spearman"), use_fast_rank = TRUE, scale = TRUE){
  method <- match.arg(method)
  
  candidates <- intersect(rownames(object), rownames(ref))
  if (method == "pearson"){
    corr <- qlcMatrix::corSparse(object[candidates,], ref[candidates,])
  } else if (method == "spearman"){
    if (require(presto, quietly = T) & use_fast_rank){
      ranked_data <- presto::rank_matrix(object[candidates,])$X_ranked
    } else{
      ranked_data <- rank_input_matrix(object[candidates,])
    }
    corr <- qlcMatrix::corSparse(ranked_data, ref[candidates,])
  }
  if (scale)
    corr <- t(scale(t(corr)))
  rownames(corr) <- colnames(object)
  colnames(corr) <- colnames(ref)
  return(corr)
}



#'@param as_assay When it is TRUE, the output is returned as an Assay object in the Seurat object
#'@param assay.name Assay name of the RSS representation in the returned Seurat object
#'@param reduction.name Reduction name of the RSS representation in the returned Seurat object
#'@param reduction.key Reduction key of the RSS representation in the returned Seurat object
#'@rdname ref_sim_spectrum
#'@method ref_sim_spectrum Seurat
ref_sim_spectrum.Seurat <- function(object, ref, as_assay = FALSE, assay.name = "rss", reduction.name = "rss", reduction.key = "RSS_", ...){
  input <- object@assays[[DefaultAssay(object)]]@data
  rss <- ref_sim_spectrum.default(input, ref, ...)
  if (as_assay){
    object[[assay.name]] <- CreateAssayObject(data = t(rss))
  } else{
    colnames(rss) <- NULL
    object[[reduction.name]] <- CreateDimReducObject(rss, key = reduction.key, assay = DefaultAssay(object))
  }
  return(object)
}


