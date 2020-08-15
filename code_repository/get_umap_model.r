get_umap_model.default <- function(
  object,
  n.neighbors = 30L,
  n.components = 2L,
  metric = 'cosine',
  n.epochs = NULL,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  n_threads = 1,
  verbose = TRUE
) {
  require(uwot)
  umap(
    X = object,
    n_threads = n_threads,
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    fast_sgd = uwot.sgd,
    ret_model = TRUE,
    verbose = verbose
  )
}

get_umap_model.Seurat <- function(
  object,
  dims = 1:20,
  reduction = 'pca', ...
) {
  require(uwot)
  data.use <- Embeddings(object[[reduction]])[, dims]
  get_umap_model.default(data.use, ...)
}

get_umap_model <- function(object, ...) {
  UseMethod(generic = 'get_umap_model', object = object)
}

