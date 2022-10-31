#' Transfer labels from the reference to the query data set
#'
#' A kNN classifier is used to predict labels of the query data, given the reference
#' data with the same data representation
#'
#'@param data_ref Reference data matrix
#'@param data_query Query data matrix
#'@param knn_ref_query KNN matrix, where each row represents one sample in the query data, and each column represents one nearest neighbors in the reference data. If NULL, it is calculated using data_ref and data_query
#'@param label_ref Class labels of each sample in the reference
#'@param k Number of neighbors for classification
#'@param thres_prop_match If no class gets higher supported proportion than this value, the prediction fails for the query sample and NA is returned
#'@param return_all Whether to return all results including the transferred labels and the kNN proportions. Only the transferred labels are returned when it is FALSE
#'@param verbose Whether to output the verbose message
#'@return A vector of the predicted/transferred labels of the query data
#'@export
transfer_labels <- function(data_ref = NULL,
                            data_query = NULL,
                            knn_ref_query = NULL,
                            label_ref,
                            k = 50,
                            thres_prop_match = 0.5,
                            return_all = FALSE,
                            verbose = TRUE){
  if (is.null(knn_ref_query)){
    if (verbose)
      message("Calculating kNN of the query data to the reference data...")
    knn_ref_query <- RANN::nn2(data_ref, data_query, k = k)$nn.idx
  }
  if (! is.factor(label_ref))
    label_ref <- as.factor(label_ref)
  
  if (verbose)
    message("Converting kNN and label information to identity matrices...")
  knn_mat <- sparseMatrix(i = rep(1:nrow(knn_ref_query), ncol(knn_ref_query)),
                          j = as.numeric(knn_ref_query),
                          x = 1,
                          dims = c(nrow(knn_ref_query), length(label_ref)))
  label_ref_mat <- sparseMatrix(i = 1:length(label_ref),
                                j = as.numeric(label_ref),
                                x = 1,
                                dims = c(length(label_ref),length(levels(label_ref))), dimnames = list(names(label_ref),levels(label_ref)))
  
  if (verbose)
    message("Projecting labels by kNN voting...")
  knn_label_mat <- knn_mat %*% label_ref_mat
  df_knn_label <- Matrix::summary(knn_label_mat)
  label_query_proj <- colnames(knn_label_mat)[sapply(split(df_knn_label, df_knn_label$i), function(x) ifelse(max(x$x) < sum(x$x) * thres_prop_match, NA, x$j[which.max(x$x)]) )]
  
  if (verbose)
    message("Done. Returning results...")
  if (return_all){
    label_max_prop <- colnames(knn_label_mat)[sapply(split(df_knn_label, df_knn_label$i), function(x) x$j[which.max(x$x)] )]
    max_prop <- sapply(split(df_knn_label, df_knn_label$i), function(x) max(x$x)/sum(x$x))
    mat_props <- as.matrix(rowNorm(knn_label_mat))
    res <- data.frame(label_query_proj, label_max_prop, max_prop, mat_props) %>%
      setNames(c("transferred_label","max_knn_label","max_knn_prop",paste0("prop_",colnames(knn_label_mat))))
    return(res)
  } else{
    return(label_query_proj)
  }
}

#' Estimate likelihood of projection failure
#' 
#' Calculate the normalized distances of each query sample to the projected reference subset,
#' and use a Gaussian mixture model to estimate the likelihood of projection failure, which is
#' featured by large distance to the reference
#' 
#' @param data_ref Reference data matrix
#' @param data_query Query data matrix
#' @param k Number of neighbors
#' @param seed The seed value for the Gaussian mixture model estimate
#' @param do_plot Whether to plot the fitted mixture model (as density plot)
#' @param verbose Whether to output the verbose message
#' @return A data frame with normalized projection distance and estimated projection failure likelihood
#' @export
estimate_projection_failure <- function(data_ref,
                                        data_query,
                                        k = 50,
                                        seed = NULL,
                                        do_plot = FALSE,
                                        verbose = TRUE){
  if (verbose)
    message("Calculating kNN of the reference data...")
  knn2ref <- RANN::nn2(data = data_ref, query = data_query, k = k)
  if (verbose)
    message("Calculating kNN of the query data to the reference data...")
  knnref <- RANN::nn2(data = data_ref, k = k)
  
  if (verbose)
    message("Calculating normalized distances of query data to the projected reference subset...")
  avg_dist2self_ref <- rowMeans(knnref$nn.dists)
  avg_dist2ref <- rowMeans(knn2ref$nn.dists)
  norm_dist2ref <- avg_dist2ref / rowMeans(apply(knn2ref$nn.idx, 2, function(idx) avg_dist2self_ref[idx]))
  if (verbose)
    message("Estimating the Gaussian mixture model...")
  set.seed(seed)
  m_norm_dist <- mixtools::normalmixEM(norm_dist2ref, lambda = c(0.5,0.5), mu = c(1, max(norm_dist2ref)))
  
  if (do_plot){
    if (verbose)
      message("Plotting the fitted Gaussian mixture model as density plot...")
    mixtools::plot.mixEM(m_norm_dist, 2)
  }
  
  if (verbose)
    message("Done. Returning results...")
  return(data.frame(norm_dist = norm_dist2ref, lh_fail = m_norm_dist$posterior[,2], row.names = rownames(data_query)))
}