#'@import Matrix
#'@param dr Dimension reduction matrix used for clustering. When it is NULL, truncated PCA is run on the expression matrix for dimension reduction
#'@param dr_input Alternative expression matrix used for dimension reduction. Ignore if dr is specified
#'@param num_pcs_compute Number of PCs to calculate. Ignore if dr is specified
#'@param num_pcs_use Number of PCs used for clustering
#'@param labels Labels specifying different samples
#'@param redo_pca If TRUE, PCA is rerun for each sample separately for clustering
#'@param k Number of nearest neighbors of the kNN network used for clustering
#'@param min_batch_size The minimal cell number of a batch to be clustered to generate references
#'@param ... Other parameters to build_knn_graph
#'@param cluster_method Method used to apply clustering to the kNN network. By default it calls FindClusters in Seurat using Louvain method. Alternative method is the walktrap community identification algorithm in igraph
#'@param cluster_resolution Resolution of clustering. Ignore if cluster_method is not Seurat
#'@param min_cluster_num The minimal number of clusters to include a sample in the ref profile (default=3)
#'@param spectrum_type Method to normalize similarities. "corr_ztransform" uses z-transform; "corr_kernel" introduces correlation kernel to convert similarities to likelihood; "corr_raw" uses no normalization; "nnet" and "lasso" build probabilistic prediction model on the data and estimate likelihoods
#'@param corr_method Type of correlation. Ignore if spectrum_type is "nnet" or "lasso"
#'@param use_fast_rank When the presto package is available, use its rank_matrix function to rank sparse matrix
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
                                         min_batch_size = k * 2,
                                         ..., # how to separate samples and whether or not to do DR separately
                                         cluster_method = c("Seurat","walktrap"),
                                         cluster_resolution = 0.6,
                                         min_cluster_num = 3,
                                         spectrum_type = c("corr_ztransform","corr_kernel","corr_raw","nnet","lasso"), # clustering and types of spectrum
                                         corr_method = c("spearman","pearson"),
                                         use_fast_rank = TRUE,
                                         lambda = 50,
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
      message("No dimension reduction is provided. Start to do truncated PCA...")
    
    t_pca <- irlba::irlba(t(dr_input), nv = num_pcs_compute)
    dr <- t_pca$u %*% diag(t_pca$d)
    dr <- dr[,1:num_pcs_use]
    rownames(dr) <- colnames(data)
    
    if (verbose)
      cat("PCA finished.\n")
  }
  
  if (verbose)
    message("Start to do clustering for each sample...")
  labels <- as.factor(labels)
  batches <- names(which(table(labels) >= min_batch_size))
  
  cl <- lapply(batches, function(x){
    idx <- which(labels == x)
    dr_x <- dr[idx,]
    if (redo_pca){
      if (verbose)
        message(paste0("  Redoing truncated PCA on sample ", x, "..."))
      t_pca <- irlba::irlba(t(dr_input[,idx]), nv = num_pcs_compute)
      dr_x <- t_pca$u %*% diag(t_pca$d)
      dr_x <- dr_x[,1:num_pcs_use]
      rownames(dr_x) <- colnames(data)[idx]
    }
    knn <- build_knn_graph(t(dr_x), k = k, ..., verbose = verbose > 1)
    rownames(knn) <- colnames(data)[idx]
    colnames(knn) <- colnames(data)[idx]
    
    if (cluster_method == "Seurat" && requireNamespace("Seurat", quietly=T)){
      cl <- Seurat::FindClusters(Seurat::as.Graph(knn), resolution = cluster_resolution, verbose = verbose > 1)[,1]
    } else if (requireNamespace("igraph", quietly=T)){
      graph <- igraph::graph_from_adjacency_matrix(knn, mode = "undirected", weighted = T)
      cl <- igraph::walktrap.community(graph)
      cl <- as.factor(setNames(cl$membership, cl$names)[rownames(knn)])
    } else{
      stop("At least one of Seurat and igraph should be installed.")
    }
    
    if (verbose)
      message(paste0("  Done clustering of sample ", x, "."))
    return(cl)
  })
  names(cl) <- batches
  if (verbose)
    message("Finished clustering.")
  
  idx_toofew_cl_batches <- which(sapply(cl, function(x) length(levels(x)) < min_cluster_num))
  if (length(idx_toofew_cl_batches) > 0){
    if (verbose)
      message(paste0("The following batch(es) have < ",min_cluster_num," clusters and are removed from the ref: ", paste(batches[idx_toofew_cl_batches], collapse=", ")))
    batches <- batches[-idx_toofew_cl_batches]
    cl <- cl[-idx_toofew_cl_batches]
  }
  
  if (spectrum_type == "nnet"){
    if (! requireNamespace("nnet", quietly=T)){
      warning("cannot find package nnet, switch spectrum type to corr_ztransform")
      spectrum_type <- "corr_ztransform"
    }
  } else if (spectrum_type == "lasso"){
    if (! requireNamespace("glmnet", quietly=T)){
      warning("cannot find package glmnet, switch spectrum type to corr_ztransform")
      spectrum_type <- "corr_ztransform"
    }
  }
  
  if (spectrum_type %in% c("corr_ztransform","corr_kernel","corr_raw")){
    cl_profiles <- lapply(batches, function(x){
      idx <- which(labels == x)
      profiles <- sapply(levels(as.factor(cl[[x]])), function(cl_x)
        apply(data[,idx[as.factor(cl[[x]])==cl_x]], 1, mean))
      return(profiles)
    })
    names(cl_profiles) <- batches
    if (verbose)
      message("Obtained average profiles of clusters.")
    
    if (verbose)
      message("Start to calculate standardized similarities to clusters...")
    sim2profiles <- lapply(cl_profiles, function(profiles){
      if (corr_method == "pearson"){
        sims <- qlcMatrix::corSparse(data, Matrix(profiles))
      } else if (corr_method == "spearman"){
        if (use_fast_rank && requireNamespace("presto", quietly=T)){
          ranked_data <- presto::rank_matrix(data)$X_ranked
        } else{
          ranked_data <- rank_input_matrix(data)
        }
        sims <- qlcMatrix::corSparse(ranked_data, profiles)
      }
      return(sims)
    })
    
    if (spectrum_type == "corr_ztransform"){
      if (verbose)
        message("Doing z-transformation...")
      sim2profiles <- lapply(sim2profiles, function(sim)
        t(scale(t(sim))))
    } else if (spectrum_type == "corr_kernel"){
      if (verbose)
        message("Doing kernel transformation...")
      sim2profiles <- lapply(sim2profiles, function(sim)
        t(apply(exp(sim * lambda) * exp(-lambda), 1, function(x) x/sum(x))))
    }
  } else if (spectrum_type %in% c("nnet","lasso")){
    if (verbose)
      message("Start to build multinomial logistic regression models...")
    
    models <- lapply(batches, function(x){
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
          message(paste0("  Randomly selected cells for model training: sample ", x, "."))
        
      } else if (train_on == "pseudo"){
        pseudo_idx <- lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          knn_cl <- build_knn_graph(t(dr[cl_idx,]), k = k_pseudo,
                                    use_seurat_snn = F, mutual = F, jaccard_weighted = F, jaccard_prune = 0,
                                    ..., verbose = verbose > 1)
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
          message(paste0("  Constructed pseudocells for training: sample ", x, "."))
      }
      idx <- which(train_y %in% unique(train_y)[sapply(unique(train_y), function(x) sum(train_y == x)) > 2])
      if (verbose & length(idx) != length(train_y))
        message(paste0("  Dropped clusters ",unique(train_y[-idx]), " in sample ", x, " due to too few training data."))
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
        message(paste0("  Model trained for sample ", x, "."))
      return(m)
    })
    names(models) <- batches
    if (verbose)
      message("Models trained, start to produce likelihood spectrum...")
    
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
      message("Done likelihood estimation.")
  }
  
  sim2profiles <- do.call(cbind, sim2profiles)
  sim2profiles_raw <- sim2profiles
  if (merge_spectrums){
    if (verbose)
      message("Start to merge similar spectrums...")
    
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
    message("Done.")
  
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
cluster_sim_spectrum.Seurat <- function(object,
                                        var_genes = NULL,
                                        use_scale = F,
                                        use_dr = "pca",
                                        dims_use = 1:20,
                                        label_tag,
                                        redo_pca = FALSE,
                                        redo_pca_with_data = FALSE,
                                        k = 20, min_batch_size = k*2,
                                        ...,
                                        cluster_resolution = 0.6,
                                        spectrum_type = c("corr_ztransform","corr_kernel","corr_raw","nnet","lasso"),
                                        corr_method = c("spearman","pearson"),
                                        lambda = 50,
                                        threads = 1,
                                        train_on = c("raw","pseudo","rand"),
                                        downsample_ratio = 1/10,
                                        k_pseudo = 10,
                                        logscale_likelihood = F,
                                        merge_spectrums = FALSE,
                                        merge_height_prop = 1/10,
                                        spectrum_dist_type = c("pearson", "euclidean"),
                                        spectrum_cl_method = "complete",
                                        reduction.name = "css",
                                        reduction.key = "CSS_",
                                        return_seuratObj = T,
                                        verbose = T){
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
  css <- cluster_sim_spectrum.default(object = data, dr = dr, dr_input = dr_input, labels = labels, redo_pca = redo_pca, k = k, min_batch_size = min_batch_size, ..., cluster_resolution = cluster_resolution,
                                      spectrum_type = spectrum_type, corr_method = corr_method, lambda = lambda, threads = threads,
                                      train_on = train_on, downsample_ratio = downsample_ratio, k_pseudo = k_pseudo, logscale_likelihood = logscale_likelihood,
                                      merge_spectrums = merge_spectrums, merge_height_prop = merge_height_prop, spectrum_dist_type = spectrum_dist_type, spectrum_cl_method = spectrum_cl_method,
                                      return_css_only = FALSE, verbose = verbose)
  if (return_seuratObj){
    object[[reduction.name]] <- CreateDimReducObject(embeddings = css$sim2profiles, key = reduction.key, assay = DefaultAssay(object), misc = list(model = css))
    return(object)
  } else{
    return(css)
  }
}

#'@param model Calculation model of the reference CSS representation
#'@param use_fast_rank When the presto package is available, use its rank_matrix function to rank sparse matrix
#'@rdname css_project
#'@export
#'@method css_project default
css_project.default <- function(object,
                                model,
                                use_fast_rank = TRUE){
  model <- model$model
  
  if (model$spectrum_type %in% c("corr_ztransform","corr_raw")){
    res <- do.call(cbind, lapply(model$profiles, function(ref)
      ref_sim_spectrum.default(object, ref, method = model$args['corr_method'], use_fast_rank = use_fast_rank, scale = model$spectrum_type == "corr_ztransform")))
    
  } else if (model$spectrum_type == "corr_kernel"){
    res <- do.call(cbind, lapply(model$profiles, function(ref){
      sim <- ref_sim_spectrum.default(object, ref, method = model$args['corr_method'], use_fast_rank = use_fast_rank, scale = FALSE)
      return(t(apply(exp(sim * model$args['lambda']) * exp(-model$args['lambda']), 1, function(x) x/sum(x))))
    }))
    
  } else if (model$spectrum_type %in% c("nnet","lasso")){
    object <- as.matrix(object)
    res <- do.call(cbind, lapply(model$models, function(m){
      if (model$spectrum_type == "lasso"){
        require(glmnet, quietly=T)
        pred <- predict(m, t(object), type = "response")[,,1]
      } else if (model$spectrum_type == "nnet"){
        require(nnet, quietly=T)
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
                               use_fast_rank = TRUE,
                               reduction.name = "css_proj",
                               reduction.key = "CSSPROJ_"){
  dat <- object@assays[[DefaultAssay(object)]]@data
  css_proj <- css_project(dat, model, use_fast_rank = use_fast_rank)
  rownames(css_proj) <- colnames(object)
  object[[reduction.name]] <- CreateDimReducObject(css_proj, key = reduction.key, assay = DefaultAssay(object))
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
  if (is.null(knn_ref_query))
    knn_ref_query <- RANN::nn2(data_ref, data_query, k = k)$nn.idx
  if (! is.factor(label_ref))
    label_ref <- as.factor(label_ref)
  
  knn_mat <- sparseMatrix(i = rep(1:nrow(knn_ref_query), ncol(knn_ref_query)),
                          j = as.numeric(knn_ref_query),
                          x = 1,
                          dims = c(nrow(knn_ref_query), length(label_ref)))
  label_ref_mat <- sparseMatrix(i = 1:length(label_ref),
                                j = as.numeric(label_ref),
                                x = 1,
                                dims = c(length(label_ref),length(levels(label_ref))), dimnames = list(names(label_ref),levels(label_ref)))
  
  knn_label_mat <- knn_mat %*% label_ref_mat
  df_knn_label <- Matrix::summary(knn_label_mat)
  label_query_proj <- colnames(knn_label_mat)[sapply(split(df_knn_label, df_knn_label$i), function(x) ifelse(max(x$x) < sum(x$x) * thres_prop_match, NA, x$j[which.max(x$x)]) )]
  
  return(label_query_proj)
}

