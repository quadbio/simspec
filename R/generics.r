#' Calculate Cluster Similarity Spectrum
#'
#' Calculate Cluster Similarity Spectrum (CSS), given expression
#' of the data and cell labels used to distinguish samples. Clustering
#' is applied to cells of each sample separately, similarities of one cell
#' to those clusters are calculated and normalized.
#'
#'@param object An object
#'@rdname cluster_sim_spectrum
#'@export cluster_sim_spectrum
cluster_sim_spectrum <- function(object, ...) {
  UseMethod(generic = 'cluster_sim_spectrum', object = object)
}

#' Prject data to the CSS representation using the given calculation model
#'
#' Calculate CSS representation using a given model. It is used to project
#' new data to the existing reference atlas which is already represented by CSS
#'
#'@param object An object
#'@rdname css_project
#'@export css_project
css_project <- function(object, ...) {
  UseMethod(generic = 'css_project', object = object)
}

#' Calculate Reference Similarity Spectrum
#'
#' Calculate Reference Similarity Spectrum (RSS), given expression of the
#' data and expression of the external reference data set.
#'
#'@param object An object
#'@rdname ref_sim_spectrum
#'@export ref_sim_spectrum
ref_sim_spectrum <- function(object, ...){
  UseMethod(generic = 'ref_sim_spectrum', object = object)
}

#' Principal component analysis on dimension reduced data
#' 
#' Do PCA analysis on dimension reduced data to further reduce
#' data dimensionality
#' 
#' @param object An object
#' @rdname run_PCA
#' @export run_PCA
run_PCA <- function(object, ...){
  UseMethod(generic = 'run_PCA', object = object)
}

#' Regress out variables from cell embeddings
#' 
#' Regress out variables from cell embeddings with a linear model
#' 
#' @param object An object
#' @rdname regress_out_from_embeddings
#' @export regress_out_from_embeddings
regress_out_from_embeddings <- function(object, ...){
  UseMethod(generic = 'regress_out_from_embeddings', object = object)
}
