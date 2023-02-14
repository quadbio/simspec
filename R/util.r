#' Build kNN network
#'
#' Calculate the kNN-network given a data matrix
#'
#'@import Matrix
#'@param data The data matrix
#'@param k The number of neighbors
#'@param dist_type The type of distance. Should be one of "euclidean", "pearson" and "raw"
#'@param use_seurat_snn If TRUE, the FindNeighbors function in Seurat will be called to calculate
#'@param mutual Only consider mutual nearest neighbors. Ignore if use_seurat_snn is TRUE
#'@param jaccard_weighted If TRUE, edges of the resulted graph are weighted by Jaccard index. Ignore if use_seurat_snn is TRUE
#'@param jaccard_prune Prune an edge if its weight is smaller than this value. Ignore if use_seurat_snn is FALSE and jaccard_weighted is FALSE
#'@param return_igraph if TRUE, an igraph object instead of the adjacent matrix will be returned
#'@param verbose If TRUE, progress message is provided
#'@return A Matrix of kNN network adjacant matrix. If return_igraph is TRUE, the corresponding igraph object is returned
#'@export
build_knn_graph <- function(data,
                            k = 20,
                            dist_type = c("euclidean", "pearson", "raw"),
                            use_seurat_snn = TRUE,
                            mutual = TRUE,
                            jaccard_weighted = TRUE,
                            jaccard_prune = 1/15,
                            return_igraph = FALSE,
                            verbose = TRUE){
  dist_type <- match.arg(dist_type)
  
  nn <- matrix(0, nrow = ncol(data), ncol = k)
  if (dist_type == "raw"){
    diag(data) <- NA
    nn <- t(apply(data, 2, function(x) order(x)[1:k]))
  } else if (dist_type %in% c("euclidean", "pearson")){
    if (dist_type == "pearson")
      data <- scale(data)
    nn <- RANN::nn2(t(data), t(data), k = k+1)$nn.idx[,-1]
  }
  if (verbose)
    message("build_knn_graph: found nearest neighbors.")
  
  i <- rep(1:nrow(nn), each = ncol(nn))
  j <- as.integer(t(nn))
  adj <- Matrix::sparseMatrix(i, j, dims = c(ncol(data), ncol(data)))
  
  if(use_seurat_snn & require(Seurat)){
    if (verbose)
      message("build_knn_graph: revoke Seurat to compute SNN.")
    
    adj <- .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn, jaccard_prune)
    if (verbose)
      message("build_knn_graph: done.")
    
  } else{
    if (mutual){
      adj <- adj * Matrix::t(adj)
    } else{
      adj <- adj + Matrix::t(adj)
      adj <- Matrix::sparseMatrix(i = adj@i, p = adj@p, dims = adj@Dim, index1 = FALSE) + 0
    }
    if (verbose)
      message("build_knn_graph: got adjacent matrix.")
    
    if (jaccard_prune > 0 & jaccard_weighted){
      if (verbose)
        message("build_knn_graph: start calculating Jaccard indices...")
      nonzero_idx <- summary(adj)[,1:2]
      
      adj_intersect <- summary(adj * (adj %*% t(adj)))
      adj_union <- colSums(adj[,adj_intersect[,1]] + adj[,adj_intersect[,2]] > 0)
      jaccard <- adj_intersect$x / adj_union
      
      remained <- which(jaccard >= jaccard_prune)
      adj <- Matrix::sparseMatrix(i = nonzero_idx[remained,1], j = nonzero_idx[remained,2], x = jaccard[remained], dims = c(nrow(nn), nrow(nn))) + 0
      if (verbose)
        message("build_knn_graph: done pruning.")
    }
  }
  adj@Dimnames <- list(colnames(data), colnames(data))
  
  if (verbose)
    message("build_knn_graph: done and returning result...")
  if (return_igraph) return(igraph::graph_from_adjacency_matrix(adj, mode = "undirected"))
  return(adj)
}

#' Construct pseudocells by averaging transcriptome of similar cells
#'
#' This function implements the procedure to average the expression profiles
#' of similar cells in order to construct pseudocells
#'
#'@import Matrix
#'@param knn_adj The adjacent matrix of the kNN network
#'@param ratio The downsample rate
#'@param init_dist Maximal distance on the kNN network to group a cell to the selected capital or pseudocell center
#'@param max_dist A cell which is not grouped to any pseudocell territory is grouped to the closest one if its distance to the capital is no more than max_dist
#'@param min_pooled_cells Only pseudocells covering at least this number of cells are valid
#'@param min_seed_num The minimum number of seed capitals to start the procedure
#'@param seed The base seed for randomly assigning a cell to a capital candidate
#'@return A data.frame indicating indices of cells and pseudocells
#'@export
construct_pseudocells <- function(knn_adj,
                                  ratio,
                                  init_dist = 1,
                                  max_dist = 2,
                                  min_pooled_cells = 2,
                                  min_seed_num = 1,
                                  seed = 12345){
  rownames(knn_adj) <- 1:nrow(knn_adj)
  colnames(knn_adj) <- 1:ncol(knn_adj)
  set.seed(seed)
  capitals <- which(runif(nrow(knn_adj)) <= ratio)
  if (length(capitals) == 0){
    set.seed(seed)
    capitals <- sample(1:nrow(knn_adj), min_seed_num)
  }
  
  graph <- igraph::graph_from_adjacency_matrix(knn_adj > 0, mode = "undirected")
  dist_to_capitals <- igraph::distances(graph, v = as.character(1:ncol(knn_adj)), to = as.character(capitals))
  
  selected <- dist_to_capitals <= init_dist
  selected <- t(sapply(1:nrow(selected), function(i){
    sel <- selected[i,]
    if (! is.null(seed))
      set.seed(seed + i)
    if (sum(sel) == 1) return(sel)
    if (sum(sel) > 1) return(1:length(sel) == sample(which(sel),1))
    if (sum(sel) == 0 & min(dist_to_capitals[i,]) <= max_dist) return(1:length(sel) == sample(which(dist_to_capitals[i,] <= max_dist),1))
    return(sel)
  }))
  if (ncol(dist_to_capitals) == 1) selected <- t(selected)
  
  sel_df <- data.frame(idx_raw = 1:nrow(knn_adj), idx_pseudo = -1, stringsAsFactors=F)
  sel_df$idx_pseudo[apply(selected, 1, sum) == 1] <- apply(data.frame(selected[apply(selected, 1, sum) == 1,]), 1, which)
  sel_df$idx_pseudo[sel_df$idx_pseudo %in% unique(sel_df$idx_pseudo)[sapply(unique(sel_df$idx_pseudo), function(x) sum(sel_df$idx_pseudo == x)) < min_pooled_cells]] <- -1
  return(sel_df)
}

#' Generate ranked sparse matrix
#' 
#' This function convert a non-negative matrix into a column-ranked matrix
#' 
#' @import Matrix
#' @param mat The input matrix, can be dense or sparse
#' @return A matrix, whose dense/sparse is the same as the input, with each column ranked
#' @export
rank_input_matrix <- function(mat){
  if (is.matrix(mat) | is.data.frame(mat)){
    ranked_mat <- apply(mat, 2, rank)
  } else{
    df_mat <- Matrix::summary(mat)
    dfs_mat <- split(df_mat, df_mat$j)
    df_mat_ranked <- do.call(rbind, lapply(dfs_mat, function(df){
      num_zeros <- nrow(mat) - nrow(df)
      ranks_nonzero <- rank(df$x)
      df$x <- ranks_nonzero + num_zeros - (1 + num_zeros) / 2
      return(df)
    }))
    ranked_mat <- sparseMatrix(i = df_mat_ranked$i, j = df_mat_ranked$j, x = df_mat_ranked$x, dims = dim(mat), dimnames = dimnames(mat))
  }
  return(ranked_mat)
}

#' Normalize each row of the given matrix to be sum of 1
#'
#' This function normalizes the given matrix in a row-wise manner so that
#' the sum of every row is one
#'
#' @import Matrix
#' @param mat The input matrix
#' @return The normalized matrix
#' @export
rowNorm <- function(mat){
  diag_mat <- Diagonal(x = 1/rowSums(mat))
  res <- diag_mat %*% mat
  dimnames(res) <- dimnames(mat)
  return(res)
}

