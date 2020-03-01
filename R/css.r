#' Build kNN network
#'
#' Calculate the kNN-network given a data matrix
#'
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
    cat("found nearest neighbors.\n")
  
  i <- rep(1:nrow(nn), each = ncol(nn))
  j <- as.integer(t(nn))
  adj <- Matrix::sparseMatrix(i, j, dims = c(ncol(data), ncol(data)))
  
  if(use_seurat_snn & require(Seurat)){
    if (verbose)
      cat("revoke Seurat to compute SNN.\n")
    
    adj <- .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn, jaccard_prune)
    if (verbose)
      cat("done.\n")
    
  } else{
    if (mutual){
      adj <- adj * t(adj)
    } else{
      adj <- adj + t(adj)
      adj <- Matrix::sparseMatrix(i = adj@i, p = adj@p, dims = adj@Dim, index1 = FALSE) + 0
    }
    if (verbose)
      cat("got adjacent matrix.\n")
    
    if (jaccard_prune > 0 & jaccard_weighted){
      if (verbose)
        cat("start calculating Jaccard indices...\n")
      nonzero_idx <- summary(adj)[,1:2]
      
      adj_intersect <- summary(adj * (adj %*% t(adj)))
      adj_union <- colSums(adj[,adj_intersect[,1]] + adj[,adj_intersect[,2]] > 0)
      jaccard <- adj_intersect$x / adj_union
      
      remained <- which(jaccard >= jaccard_prune)
      adj <- Matrix::sparseMatrix(i = nonzero_idx[remained,1], j = nonzero_idx[remained,2], x = jaccard[remained], dims = c(nrow(nn), nrow(nn))) + 0
      if (verbose)
        cat("done pruning.\n")
    }
  }
  adj@Dimnames <- list(colnames(data), colnames(data))
  
  if (verbose)
    cat("done. returning result...\n")
  if (return_igraph) return(igraph::graph_from_adjacency_matrix(adj, mode = "undirected"))
  return(adj)
}

#' Construct pseudocells by averaging transcriptome of similar cells
#'
#' This function implements the procedure to average the expression profiles
#' of similar cells in order to construct pseudocells
#'
#'@param knn_adj The adjacent matrix of the kNN network
#'@param ratio The downsample rate
#'@param init_dist Maximal distance on the kNN netowkr to group a cell to the selected capital or pseudocell center
#'@param max_dist A cell which is not grouped to any pseudocell territory is grouped to the closest one if its distance to the capital is no more than max_dist
#'@param min_pooled_cells Only pseudocells covering at least this number of cells are valid
#'@param min_seed_num The minimum number of seed capitals to start the procedure
#'@return A data.frame indicating indices of cells and pseudocells
#'@export
construct_pseudocells <- function(knn_adj,
								  ratio,
								  init_dist = 1,
								  max_dist = 2,
								  min_pooled_cells = 2,
								  min_seed_num = 1){
  rownames(knn_adj) <- 1:nrow(knn_adj)
  colnames(knn_adj) <- 1:ncol(knn_adj)
  capitals <- which(runif(nrow(knn_adj)) <= ratio)
  if (length(capitals) == 0) capitals <- sample(1:nrow(knn_adj), min_seed_num)
  
  graph <- igraph::graph_from_adjacency_matrix(knn_adj > 0, mode = "undirected")
  dist_to_capitals <- igraph::distances(graph, v = as.character(1:ncol(knn_adj)), to = as.character(capitals))
  
  selected <- dist_to_capitals <= init_dist
  selected <- t(sapply(1:nrow(selected), function(i){
    sel <- selected[i,]
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

#'@param dr Dimension reduction matrix used for clustering. When it is NULL, truncated PCA is run on the expression matrix for dimension reduction
#'@param dr_input Alternative expression matrix used for dimension reduction. Ignore if dr is specified
#'@param num_pcs_compute Number of PCs to calculate. Ignore if dr is specified
#'@param num_pcs_use Number of PCs used for clustering
#'@param labels Labels specifying different samples
#'@param redo_pca If TRUE, PCA is rerun for each sample separately for clustering
#'@param k Number of nearest neighbors of the kNN network used for clustering
#'@param ... Other parameters to build_knn_graph
#'@param cluster_method Method used to apply clustering to the kNN network. By default it calls FindClusters in Seurat using Louvain method. Alternative method is the walktrap community identification algorithm in igraph
#'@param cluster_resolution Resolution of clustering. Ignore if cluster_method is not Seurat
#'@param spectrum_type Method to normalize similarities. "corr_ztransform" uses z-transform; "corr_kernel" introduces correlation kernel to convert similarities to likelihood; "corr_raw" uses no normalization; "nnet" and "lasso" build probabilistic prediction model on the data and estimate likelihoods
#'@param corr_method Type of correlation. Ignore if spectrum_type is "nnet" or "lasso"
#'@param lambda Lambda in the correlation kernel
#'@param threads Number of threads to use. Only useful if spectrum_type is "lasso"
#'@param train_on Type of data used to train the likelihood model. Only useful if spectrum_type is "nnet" or "lasso"
#'@param downsample_ratio Downsample rate. Only useful if train_on is "pseudo" or "rand"
#'@param k_pseudo Number of nearest neighbors used to construct pseudocells. Only useful if train_on is "pseudo"
#'@param logscale_likelihood If TRUE, estimated likelihoods are log-transformed. Ignore if spectrum_type is "corr_ztransform" or "corr_raw"
#'@param merge_spectrum If TRUE, similar similarity spectrums are averaged
#'@param merge_height_prop The height of dendrogram to cut. Ignore if merge_spectrum is FALSE
#'@param spectrum_dist_type Type of distance to construct the dendrogram of spectrums. Ignore if merge_spectrum is FALSE
#'@param spectrum_cl_method Method of hierarchical clustering to construct the dendrogram of spectrums. Ignore if merge_spectrum is FALSE
#'@param return_css_only If FALSE, not only the calculated CSS matrix, but also other information to recalculate the spectrum is returned
#'@param verbose If TRUE, progress message is provided
#'
#'@rdname cluster_sim_spectrum
#'@export
#'@method cluster_sim_spectrum default
cluster_sim_spectrum.default <- function(object, # expression matrix
                                         dr = NULL,
										 dr_input = NULL,
										 num_pcs_compute = 50,
										 num_pcs_use = 20, # for dimension reduction
                                         labels,
										 redo_pca = FALSE,
										 k = 20,
										 ..., # how to separate samples and whether or not to do DR separately
                                         cluster_method = c("Seurat","walktrap"),
										 cluster_resolution = 0.6,
										 spectrum_type = c("corr_ztransform","corr_kernel","corr_raw","nnet","lasso"), # clustering and types of spectrum
                                         corr_method = c("spearman","pearson"),
										 lambda = 1,
										 threads = 1,  # spectrum related parameters
                                         train_on = c("raw","pseudo","rand"),
										 downsample_ratio = 1/10,
										 k_pseudo = 10,
										 logscale_likelihood = F, # parameters of likelihood spectrum
                                         merge_spectrums = FALSE,
										 merge_height_prop = 1/10,
										 spectrum_dist_type = c("pearson", "euclidean"),
										 spectrum_cl_method = "complete", # parameters of spectrum merging
                                         return_css_only = T,
										 verbose = T){
  spectrum_type <- match.arg(spectrum_type)
  cluster_method <- match.arg(cluster_method)
  corr_method <- match.arg(corr_method)
  train_on <- match.arg(train_on)
  spectrum_dist_type <- match.arg(spectrum_dist_type)
  data <- object
  if (is.null(dr_input))
    dr_input <- data
  
  if (is.null(dr)){
    if (verbose)
      cat("No dimension reduction is provided. Start to do truncated PCA...\n")
    
    t_pca <- irlba::irlba(t(dr_input), nv = num_pcs_compute)
    dr <- t_pca$u %*% diag(t_pca$d)
    dr <- dr[,1:num_pcs_use]
    rownames(dr) <- colnames(data)
    
    if (verbose)
      cat("PCA finished.\n")
  }
  
  if (verbose)
    cat("Start to do clustering for each sample...\n")
  labels <- as.factor(labels)
  cl <- lapply(levels(labels), function(x){
    idx <- which(labels == x)
    dr_x <- dr[idx,]
    if (redo_pca){
      if (verbose)
        cat(paste0(">> Redoing truncated PCA on sample ", x, "...\n"))
      t_pca <- irlba::irlba(t(dr_input[,idx]), nv = num_pcs_compute)
      dr_x <- t_pca$u %*% diag(t_pca$d)
      dr_x <- dr_x[,1:num_pcs_use]
      rownames(dr_x) <- colnames(data)[idx]
    }
    knn <- build_knn_graph(t(dr_x), k = k, ..., verbose = verbose)
    rownames(knn) <- colnames(data)[idx]
    colnames(knn) <- colnames(data)[idx]
    
    if (cluster_method == "Seurat" & require(Seurat)){
      cl <- Seurat::FindClusters(Seurat::as.Graph(knn), resolution = cluster_resolution, verbose = verbose)[,1]
    } else if (require(igraph)){
      graph <- igraph::graph_from_adjacency_matrix(knn, mode = "undirected", weighted = T)
      cl <- igraph::walktrap.community(graph)
      cl <- as.factor(setNames(cl$membership, cl$names)[rownames(knn)])
    } else{
      stop("At least one of Seurat and igraph should be installed.")
    }
    
    if (verbose)
      cat(paste0(">> Done clustering of sample ", x, ".\n"))
    return(cl)
  })
  names(cl) <- levels(labels)
  if (verbose)
    cat("Finished clustering.\n")
  
  if (spectrum_type %in% c("corr_ztransform","corr_kernel","corr_raw")){
    cl_profiles <- lapply(levels(labels), function(x){
      idx <- which(labels == x)
      profiles <- sapply(levels(as.factor(cl[[x]])), function(cl_x)
        apply(data[,idx[as.factor(cl[[x]])==cl_x]], 1, mean))
      return(profiles)
    })
    names(cl_profiles) <- levels(labels)
    if (verbose)
      cat("Obtained average profiles of clusters.\n")
    
    if (verbose)
      cat("Start to calculate standardized similarities to clusters...\n")
    sim2profiles <- lapply(cl_profiles, function(profiles){
      cor(as.matrix(data), profiles, method=corr_method)
    })
    
    if (spectrum_type == "nnet"){
      if (! require(nnet)){
        warning("cannot find package nnet, switch spectrum type to corr_ztransform")
        spectrum_type <- "corr_ztransform"
      }
    } else if (spectrum_type == "lasso"){
      if (!require(glmnet)){
        warning("cannot find package glmnet, switch spectrum type to corr_ztransform")
        spectrum_type <- "corr_ztransform"
      }
    }
    
    if (spectrum_type == "corr_ztransform"){
      if (verbose)
        cat("Doing z-transformation...\n")
      sim2profiles <- lapply(sim2profiles, function(sim)
        t(scale(t(sim))))
    } else if (spectrum_type == "corr_kernel"){
      if (verbose)
        cat("Doing kernel transformation...\n")
      sim2profiles <- lapply(sim2profiles, function(sim)
        t(apply(exp(sim * lambda) * exp(-lambda), 1, function(x) x/sum(x))))
    }
  } else if (spectrum_type %in% c("nnet","lasso")){
    if (verbose)
      cat("Start to build multinomial logistic regression models...\n")
    
    models <- lapply(levels(labels), function(x){
      idx_x <- which(labels == x)
      cl_x <- cl[[x]]
      train_x <- t(data[,idx_x])
      train_y <- cl_x
      if (train_on == "rand"){
        sel_idx <- lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          sel_idx <- cl_idx[which(runif(length(cl_idx)) <= downsample_ratio)]
          if (length(sel_idx) == 0) sel_idx <- sample(cl_idx, 1)
          return(sel_idx)
        })
        train_x <- as.matrix(do.call(rbind, lapply(sel_idx, function(idx) t(data[,idx]))))
        train_y <- as.factor(rep(levels(cl_x), sapply(sel_idx, length)))
        if (verbose)
          cat(paste0(">> Randomly selected cells for model training: sample ", x, ".\n"))
        
      } else if (train_on == "pseudo"){
        pseudo_idx <- lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          knn_cl <- build_knn_graph(t(dr[cl_idx,]), k = k_pseudo,
                                    use_seurat_snn = F, mutual = F, jaccard_weighted = F, jaccard_prune = 0,
                                    ..., verbose = F)
          pseudo_idx <- construct_pseudocells(knn_cl, ratio = downsample_ratio, min_seed_num = 2)
          return(pseudo_idx)
        })
        names(pseudo_idx) <- levels(cl_x)
        train_x <- do.call(rbind, lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          idx <- pseudo_idx[[cl_this]]
          data_pseudo <- sapply(sort(unique(idx[idx[,2]!=-1,2])), function(ix) apply(data[,cl_idx[idx[,2] == ix]], 1, mean) )
          return(t(data_pseudo))
        }))
        train_y <- rep(levels(cl_x), sapply(pseudo_idx, function(idx) length(unique(idx[idx[,2]!=-1,2]))))
        if (verbose)
          cat(paste0(">> Constructed pseudocells for training: sample ", x, ".\n"))
      }
      idx <- which(train_y %in% unique(train_y)[sapply(unique(train_y), function(x) sum(train_y == x)) > 2])
      if (verbose & length(idx) != length(train_y))
        cat(paste0(">> Dropped clusters ",unique(train_y[-idx]), " in sample ", x, " due to too few training data.\n"))
      train_x <- train_x[idx,]
      train_y <- train_y[idx]
      
      # train model: nnet::multinom
      m <- NULL
      if (spectrum_type == "nnet"){
        train_dat <- data.frame(y = train_y, train_x)
        m <- nnet::multinom(y ~ ., data = train_dat, MaxNWts = 1E8, trace = verbose)
      } else if (spectrum_type == "lasso"){
        if (threads > 1){
          library(doMC)
          registerDoMC(threads)
        }
        m <- glmnet::cv.glmnet(x = train_x, y = train_y, family = "multinomial", alpha = 1)
      }
      if (verbose)
        cat(paste0(">> Model trained for sample ", x, ".\n"))
      return(m)
    })
    names(models) <- levels(labels)
    if (verbose)
      cat("Models trained, start to produce likelihood spectrum...\n")
    
    sim2profiles <- lapply(models, function(m){
      if (spectrum_type == "lasso"){
        pred <- predict(m, t(data), type = "response")[,,1]
      } else if (spectrum_type == "nnet"){
        pred <- predict(m, data.frame(t(data)), "probs")
      }
      return(pred)
    })
    
    if (logscale_likelihood)
      sim2profiles <- lapply(sim2profiles, function(pred){
        pred <- t(scale(t(log(pred))))
      })
    
    if (verbose)
      cat("Done likelihood estimation.\n")
  }
  
  sim2profiles <- do.call(cbind, sim2profiles)
  sim2profiles_raw <- sim2profiles
  if (merge_spectrums){
    if (verbose)
      cat("Start to merge similar spectrums...\n")
    
    dist_css <- NULL
    if (spectrum_dist_type == "pearson"){
      dist_css <- as.dist(1 - cor(sim2profiles))
    } else if (spectrum_dist_type == "euclidean"){
      dist_css <- dist(t(sim2profiles))
    }
    cl <- hclust(dist_css, method = spectrum_cl_method)
    cl_spectrums <- cutree(cl, h = max(cl$height) * merge_height_prop)
    sim2profiles <- sapply(sort(unique(cl_spectrums)), function(x){
      if(sum(cl_spectrums==x)==1)
        return(sim2profiles[,cl_spectrums==x])
      return(apply(sim2profiles[,cl_spectrums==x], 1, mean))
    })
    #sim2profiles <- t(scale(t(sim2profiles)))
  }
  rownames(sim2profiles) <- colnames(data)
  colnames(sim2profiles) <- paste0("CSS_", 1:ncol(sim2profiles))
  
  if (verbose)
    cat("Done.\n")
  
  if (return_css_only)
    return(sim2profiles)
  
  model <- list(spectrum_type = spectrum_type)
  if (spectrum_type %in% c("corr_ztransform","corr_kernel","corr_raw")){
    model$profiles <- cl_profiles
    model$args <- c(corr_method = corr_method, lambda = lambda)
  } else if (spectrum_type %in% c("nnet","lasso")){
    model$models <- models
    model$args <- c(logscale_likelihood = logscale_likelihood)
  }
  model$merged_spectrum <- seq(1, ncol(sim2profiles_raw))
  if (merge_spectrums)
    model$merged_spectrum <- cl_spectrums
  res <- list(model = model, sim2profiles = sim2profiles)
  return(res)
}

#'@param var_genes Genes used for similarity calculation. If NULL, predefined variable features are used
#'@param use_scale If TRUE, scale.data rather than data slot is used for similarity calculation
#'@param use_dr Name of reduction used for clustering
#'@param dims_use Dimensions in the reduction used for clustering
#'@param label_tag Column in the meta.data slot showing sample labels
#'@param redo_pca_with_data If TRUE, data slot is used to redo PCA for each sample. Ignore if redo_pca is FALSE
#'@param reduction.name Reduction name of the CSS representation in the returned Seurat object
#'@param reduction.key Reduction key of the CSS representation in the returned Seurat object
#'@param return_seuratObj If TRUE, a Seurat object with CSS added as one dimension reduction representation is returned. Otherwise, a list with CSS matrix and the calculation model is returned
#'@rdname cluster_sim_spectrum
#'@export
#'@method cluster_sim_spectrum Seurat
cluster_sim_spectrum.Seurat <- function(object, var_genes = NULL, use_scale = F, use_dr = "pca", dims_use = 1:20,
                                        label_tag, redo_pca = FALSE, redo_pca_with_data = FALSE, k = 20, ...,
                                        cluster_resolution = 0.6, spectrum_type = c("corr_ztransform","corr_kernel","corr_raw","nnet","lasso"),
                                        corr_method = c("spearman","pearson"), lambda = 1, threads = 1,
                                        train_on = c("raw","pseudo","rand"), downsample_ratio = 1/10, k_pseudo = 10, logscale_likelihood = F,
                                        merge_spectrums = FALSE, merge_height_prop = 1/10, spectrum_dist_type = c("pearson", "euclidean"), spectrum_cl_method = "complete",
                                        reduction.name = "css", reduction.key = "CSS_",
                                        return_seuratObj = T, verbose = T){
  spectrum_type <- match.arg(spectrum_type)
  corr_method <- match.arg(corr_method)
  train_on <- match.arg(train_on)
  spectrum_dist_type <- match.arg(spectrum_dist_type)
  if (is.null(var_genes))
    var_genes <- Seurat::VariableFeatures(object)
  
  data <- object[[object@active.assay]]@data[var_genes,]
  if (use_scale)
    data <- object[[object@active.assay]]@scale.data[var_genes,]
  
  dr_input <- NULL
  if (redo_pca){
    if (redo_pca_with_data){
      dr_input <- object[[object@active.assay]]@data[var_genes,]
    } else{
      dr_input <- object[[object@active.assay]]@scale.data[var_genes,]
    }
  }
  
  dr <- object@reductions[[use_dr]]@cell.embeddings[,dims_use]
  labels <- object@meta.data[,label_tag]
  css <- cluster_sim_spectrum.default(object = data, dr = dr, dr_input = dr_input, labels = labels, redo_pca = redo_pca, k = k, ..., cluster_resolution = cluster_resolution,
                                      spectrum_type = spectrum_type, corr_method = corr_method, lambda = lambda, threads = threads,
                                      train_on = train_on, downsample_ratio = downsample_ratio, k_pseudo = k_pseudo, logscale_likelihood = logscale_likelihood,
                                      merge_spectrums = merge_spectrums, merge_height_prop = merge_height_prop, spectrum_dist_type = spectrum_dist_type, spectrum_cl_method = spectrum_cl_method,
                                      return_css_only = return_seuratObj, verbose = verbose)
  if (return_seuratObj){
    object[[reduction.name]] <- CreateDimReducObject(embeddings = css, key = reduction.key, assay = DefaultAssay(object))
    return(object)
  } else{
    return(css)
  }
}

#'@param model Calculation model of the reference CSS representation
#'@rdname css_project
#'@export
#'@method css_project default
css_project.default <- function(object,
                                model){
  model <- model$model
  object <- as.matrix(object)
  
  if (model$spectrum_type %in% c("corr_ztransform","corr_raw")){
    res <- do.call(cbind, lapply(model$profiles, function(ref)
      ref_sim_spectrum.default(object, ref, method = model$args['corr_method'], scale = model$spectrum_type == "corr_ztransform")))
    
  } else if (model$spectrum_type == "corr_kernel"){
    res <- do.call(cbind, lapply(model$profiles, function(ref){
      sim <- ref_sim_spectrum.default(object, ref, method = model$args['corr_method'], scale = FALSE)
      return(t(apply(exp(sim * model$args['lambda']) * exp(-model$args['lambda']), 1, function(x) x/sum(x))))
    }))
    
  } else if (model$spectrum_type %in% c("nnet","lasso")){
    res <- do.call(cbind, lapply(model$models, function(m){
      if (model$spectrum_type == "lasso"){
        require(glmnet)
        pred <- predict(m, t(object), type = "response")[,,1]
      } else if (model$spectrum_type == "nnet"){
        require(nnet)
        pred <- predict(m, data.frame(t(object)), "probs")
      }
      return(pred)
    }))
  }
  
  return(res)
}

#'@param reduction.name Reduction name of the projected CSS representation in the returned Seurat object
#'@param reduction.key Reduction key of the projected CSS representation in the returned Seurat object
#'@rdname css_project
#'@export
#'@method css_project Seurat
css_project.Seurat <- function(object,
                               model,
                               reduction.name = "css_proj",
                               reduction.key = "CSSPROJ_"){
  dat <- object@assays[[DefaultAssay(object)]]@data
  css_proj <- css_project(dat, model)
  rownames(css_proj) <- colnames(object)
  object[[reduction.name]] <- CreateDimReducObject(css_proj, key = reduction.key, assay = DefaultArray(object))
  return(object)
}

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
#'@return A vector of the predicted/transferred labels of the query data
#'@export
transfer_labels <- function(data_ref = NULL, data_query = NULL, knn_ref_query = NULL, label_ref, k = 50, thres_prop_match = 0.5){
  require(RANN)
  if (is.null(knn_ref_query))
    knn_ref_query <- RANN::nn2(data_ref, data_query, k = k)$nn.idx
  if (! is.factor(label_ref))
    label_ref <- as.factor(label_ref)
  
  proj_cl_num_query <- sapply(levels(label_ref), function(cl)
    apply(knn_ref_query, 1, function(idx) sum(label_ref[idx] == cl)))
  colnames(proj_cl_num_query) <- levels(label_ref)
  
  label_query_proj <- colnames(proj_cl_num_query)[apply(proj_cl_num_query, 1, which.max)]
  label_query_proj[which(rowMax(proj_cl_num_query) < ncol(knn_ref_query) * thres_prop_match)] <- NA
  return(label_query_proj)
}

