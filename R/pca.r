#'@param npcs Number of PCs to compute
#'@param return_rotation Whether to return the rotation matrix 
#'
#'@rdname run_PCA
#'@method run_PCA default
run_PCA.default <- function(object,
                            npcs = 50,
                            return_rotation = F)
{
  require(irlba)
  t_pca <- irlba::irlba(object, nv = npcs)
  rotation <- t_pca$v
  rownames(rotation) <- colnames(object)
  dr <- t_pca$u %*% diag(t_pca$d)
  dr <- dr[,1:npcs]
  rownames(dr) <- rownames(object)
  
  if (return_rotation)
    return(list(rotation = rotation, x = dr))
  return(dr)
}

#'@param reduction The reduction to run PCA on
#'@param npcs Number of PCs to compute
#'@param reduction.name Name of the returned reduction
#'@param reduction.key Key of the returned reduction
#'@rdname run_PCA
#'@method run_PCA Seurat
run_PCA.Seurat <- function(object,
                           reduction,
                           npcs = 10,
                           reduction.name = NULL,
                           reduction.key = NULL)
{
  reduc <- Seurat::Embeddings(object, reduction)
  dr <- run_PCA.default(reduc, npcs = npcs)
  
  if (is.null(reduction.name)) reduction.name <- paste0(reduction, "_pca")
  if (is.null(reduction.key)) reduction.key <- paste0(toupper(reduction), "PCA_")
  object[[reduction.name]] <- Seurat::CreateDimReducObject(dr, key = reduction.key, assay = Seurat::DefaultAssay(object))
  return(object)
}
