#'@param vars Data frame of the variables to regress out
#'
#'@rdname regress_out_from_embeddings
#'@export
#'@method regress_out_from_embeddings default
regress_out_from_embeddings.default <- function(object,
                                                vars)
{
  if (requireNamespace("pbapply", quietly = T)){
    apply <- pbapply::pbapply
  }
  
  emb <- apply(object, 2, function(x){
    dat <- data.frame(x_ = x, vars)
    residuals(lm(x_ ~ ., data = dat))
  })
  colnames(emb) <- colnames(object)
  return(emb)
}

#'@param reduction Name of the reduction object to use
#'@param vars_to_regress Variables in the meta.data slot for regressing out
#'@param reduction.name Name of the new reduction object
#'@param reduction,key Key of the returned reduction
#'
#'@rdname regress_out_from_embeddings
#'@export
#'@method regress_out_from_embeddings Seurat
regress_out_from_embeddings.Seurat <- function(object,
                                               reduction,
                                               vars_to_regress,
                                               reduction.name = reduction,
                                               reduction.key = NULL)
{
  emb <- Embeddings(object, reduction)
  vars <- setNames(data.frame(object@meta.data[,vars_to_regress],
                              row.names = colnames(object)),
                   vars_to_regress)
  emb_new <- regress_out_from_embeddings.default(object = emb,
                                                 vars = vars)
  colnames(emb_new) <- NULL
  object[[reduction.name]] <- CreateDimReducObject(emb_new, key = reduction.key, assay = object[[reduction]]@assay.used)
  return(object)
}