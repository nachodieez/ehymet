#' Perform hierarchical clustering for a given combination of indexes, method and distance
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars character vector of indexes names for clustering (e. g. vars1)
#' @param method The agglomerative method to be used in hierarchical clustering
#' @param dist The distance method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
#'
clustInd_hierarch_aux <- function(ind_data, vars, method ="single", dist ="euclidean",
                                  n_cluster = 2, true_labels=NULL){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(get(vars) %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  t0 <- Sys.time()

  # Perform hierarchical clustering
  d <- stats::dist(ind_data[,get(vars)], method = dist)
  met <- stats::hclust(d, method = method)
  clus <- stats::cutree(met, k = n_cluster)

  t1 <- Sys.time()

  # Calculate execution time
  et <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(et))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(et))
  }

  return(res)
}

#' Perform hierarchical clustering for a different combinations of indexes, method and distance
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars_list List of characters representing variable sets for clustering
#' @param method_list List of clustering methods
#' @param dist_list List of distance metrics
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing hierarchical clustering results for each configuration
#' @export
clustInd_hierarch <- function(ind_data, vars_list,
                              method_list = c("single","complete","average","centroid","ward.D2"),
                              dist_list = c("euclidean", "manhattan"), n_cluster=2, true_labels = NULL) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check all elements in vars_list are characters
  if(!all(is.character(vars_list))){
    stop("Elemens in vars_list should be characters")
  }

  # Check if indices, methods and distances lists are provided
  if (!is.character(vars_list) || !is.character(method_list) || !is.character(dist_list) ||
      length(vars_list) == 0 || length(method_list) == 0 || length(dist_list) == 0) {
    stop("Invalid 'method_list' or 'dist_list'. Both must be non-empty character vectors.")
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = vars_list, method = method_list, distance = dist_list)

  cl_list <- list()
  tl_null <- is.null(true_labels)
  nc <- ifelse(tl_null, 1, 4)
  metrics_df <- data.frame(matrix(ncol = nc, nrow = 0))

  for(comb in 1:nrow(parameter_combinations)){

    # Apply hierarchical clustering to each combination
    result <- clustInd_hierarch_aux(ind_data, as.character(parameter_combinations$vars[comb]),
                                    parameter_combinations$method[comb],
                                    parameter_combinations$distance[comb],
                                    n_cluster, true_labels)
    result_name <- paste("hierarch", parameter_combinations$method[comb],
                         parameter_combinations$distance[comb],
                         paste(get(as.character(parameter_combinations$vars[comb])), collapse = ""),
                         sep = "_")
    cl_list[[result_name]] <- result$cluster

    metrics_row <- result$time
    if (!tl_null) {
      metrics_row <- c(result$valid, metrics_row)
    }
    metrics_df <- rbind(metrics_df, metrics_row)
    rownames(metrics_df)[nrow(metrics_df)] <- result_name
  }

  if(tl_null) colnames(metrics_df) <- "Time"
  else colnames(metrics_df) <- c("Purity","Fmeasure","RI","Time")

  result_list <- list("cluster" = cl_list, "metrics" = metrics_df)

  return(result_list)
}


kmeans_mahal <- function(ind_data, n_cluster){

  # Check if input is numeric matrix or array
  if (!(is.numeric(ind_data) || is.matrix(ind_data) || is.array(ind_data) ||
        is.data.frame(ind_data))) {
    stop("Input must be a numeric matrix, array or data frame.")
  }

  # Convert data to matrix
  if (! is.matrix(ind_data)) {
    data_matrix <- as.matrix(ind_data)
  } else {
    data_matrix <- ind_data
  }

  # Cholesky transformation
  c <- chol(stats::var(data_matrix))
  data_transform <- data_matrix %*% solve(c)

  #vector containing the clustering partition
  km <- stats::kmeans(data_transform, centers=n_cluster, iter.max=1000, nstart=100)$cluster
  return(km)
}

#' Perform kmeans clustering for a given combination of indexes and distance
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars character vector of indexes names for clustering (e. g. vars1)
#' @param dist The distance method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_kmeans_aux <- function(ind_data, vars, dist = "euclidean",
                                n_cluster = 2, true_labels = NULL){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(get(vars) %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  # Check if the distance given can be used
  if (!dist %in% c("euclidean", "mahalanobis")) {
    stop("Invalid distance.")
  }

  t0 <- Sys.time()

  if(dist=="euclidean"){
    clus <- stats::kmeans(ind_data[,get(vars)], centers = n_cluster, iter.max = 1000, nstart = 100)$cluster
  }else{
    clus <- kmeans_mahal(ind_data[,get(vars)], n_cluster)
  }
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform hierarchical clustering for a different combinations of indexes, method and distance
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars_list List of characters representing variable sets for clustering
#' @param dist_list List of distance metrics
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing kmeans clustering results for each configuration
#' @export
clustInd_kmeans <- function(ind_data, vars_list,
                            dist_list =  c("euclidean", "mahalanobis"),
                            n_cluster = 2, true_labels = NULL) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check all elements in vars_list are characters
  if(!all(is.character(vars_list))){
    stop("Elemens in vars_list should be characters")
  }

  # Check if indices and distances lists are provided
  if (!is.character(vars_list) || !is.character(dist_list) ||
      length(vars_list) == 0 || length(dist_list) == 0) {
    stop("Invalid 'method_list' or 'dist_list'. Both must be non-empty character vectors.")
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = vars_list, distance = dist_list)

  cl_list <- list()
  tl_null <- is.null(true_labels)
  nc <- ifelse(tl_null, 1, 4)
  metrics_df <- data.frame(matrix(ncol = nc, nrow = 0))

  for(comb in 1:nrow(parameter_combinations)){

    # Apply kmeans to each combination
    result <- clustInd_kmeans_aux(ind_data, as.character(parameter_combinations$vars[comb]),
                                  parameter_combinations$distance[comb],
                                  n_cluster, true_labels)
    result_name <- paste("kmeans", parameter_combinations$distance[comb],
                         paste(get(as.character(parameter_combinations$vars[comb])), collapse = ""),
                         sep = "_")
    cl_list[[result_name]] <- result$cluster

    metrics_row <- result$time
    if (!tl_null) {
      metrics_row <- c(result$valid, metrics_row)
    }
    metrics_df <- rbind(metrics_df, metrics_row)
    rownames(metrics_df)[nrow(metrics_df)] <- result_name
  }

  if(tl_null) colnames(metrics_df) <- "Time"
  else colnames(metrics_df) <- c("Purity","Fmeasure","RI","Time")

  result_list <- list("cluster" = cl_list, "metrics" = metrics_df)

  return(result_list)
}
#' Perform kernel kmeans clustering for a given combination of indexes and distance
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars character vector of indexes names for clustering (e. g. vars1)
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_kkmeans_aux <- function(ind_data, vars, kernel = "rbfdot", n_cluster = 2,
                                 true_labels = NULL, ...){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(get(vars) %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  t0 <- Sys.time()
  kkmeans_out <- kernlab::kkmeans(as.matrix(ind_data[,get(vars)]), centers= n_cluster, kernel = kernel, ...)
  clus <- kkmeans_out@.Data
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform kernel kmeans clustering for a different combinations of indexes and kernel
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars_list List of characters representing variable sets for clustering
#' @param kernel_list List of kernels
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#' @param ... Additional arguments (unused)
#'
#' @return A list containing kkmeans clustering results for each configuration
#' @export
clustInd_kkmeans <- function(ind_data, vars_list, kernel_list = c("rbfdot", "polydot"),
                             n_cluster=2, true_labels = NULL, ...) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check all elements in vars_list are characters
  if(!all(is.character(vars_list))){
    stop("Elemens in vars_list should be characters")
  }

  # Check if indices and distances lists are provided
  if (!is.character(vars_list) || !is.character(kernel_list) ||
      length(vars_list) == 0 || length(kernel_list) == 0) {
    stop("Invalid 'method_list' or 'kernel_list'. Both must be non-empty character vectors.")
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = vars_list, kernel = kernel_list)

  cl_list <- list()
  tl_null <- is.null(true_labels)
  nc <- ifelse(tl_null, 1, 4)
  metrics_df <- data.frame(matrix(ncol = nc, nrow = 0))

  for(comb in 1:nrow(parameter_combinations)){

    # Apply kmeans to each combination
    result <- clustInd_kkmeans_aux(ind_data, as.character(parameter_combinations$vars[comb]),
                                  as.character(parameter_combinations$kernel[comb]),
                                  n_cluster, true_labels, ...)
    result_name <- paste("kkmeans", parameter_combinations$kernel[comb],
                         paste(get(as.character(parameter_combinations$vars[comb])), collapse = ""),
                         sep = "_")
    cl_list[[result_name]] <- result$cluster

    metrics_row <- result$time
    if (!tl_null) {
      metrics_row <- c(result$valid, metrics_row)
    }
    metrics_df <- rbind(metrics_df, metrics_row)
    rownames(metrics_df)[nrow(metrics_df)] <- result_name
  }

  if(tl_null) colnames(metrics_df) <- "Time"
  else colnames(metrics_df) <- c("Purity","Fmeasure","RI","Time")

  result_list <- list("cluster" = cl_list, "metrics" = metrics_df)

  return(result_list)
}

#' Perform support vector clustering for a given combination of indexes and kernel
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars_list List of characters representing variable sets for clustering
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_svc_aux <- function(ind_data, vars, method = "kmeans", n_cluster = 2,
                             true_labels = NULL, ...){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(get(vars) %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  t0 <- Sys.time()

  # create target variable
  nrow_ind_data <- nrow(ind_data)
  r1 <- round(nrow_ind_data/n_cluster)
  r2 <- nrow_ind_data - (n_cluster-1) * r1
  y <-  c(rep(1:(n_cluster-1), each =r1), rep(n_cluster, r2))
  clusterSVM_output <- SwarmSVM::clusterSVM(ind_data[,get(vars)], y, n_cluster, cluster.method=method, ...)
  clus <- clusterSVM_output$label

  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform support vector clustering for a different combinations of indexes and distance
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars_list List of characters representing variable sets for clustering
#' @param method_list List of methods
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#' @param ... Additional arguments (unused)
#'
#' @return A list containing kkmeans clustering results for each configuration
#' @export
clustInd_svc <- function(ind_data, vars_list, method_list = c("kmeans", "kernkmeans"),
                         n_cluster = 2, true_labels = NULL, ...) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check all elements in vars_list are characters
  if(!all(is.character(vars_list))){
    stop("Elemens in vars_list should be characters")
  }

  # Check if indices and distances lists are provided
  if (!is.character(vars_list) || !is.character(method_list) ||
      length(vars_list) == 0 || length(method_list) == 0) {
    stop("Invalid 'method_list' or 'method_list'. Both must be non-empty character vectors.")
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = vars_list, method = method_list)

  cl_list <- list()
  tl_null <- is.null(true_labels)
  nc <- ifelse(tl_null, 1, 4)
  metrics_df <- data.frame(matrix(ncol = nc, nrow = 0))

  for(comb in 1:nrow(parameter_combinations)){

    # Apply svc to each combination
    result <- clustInd_svc_aux(ind_data, as.character(parameter_combinations$vars[comb]),
                                   as.character(parameter_combinations$method[comb]),
                                   n_cluster, true_labels, ...)
    result_name <- paste("svc", parameter_combinations$method[comb],
                         paste(get(as.character(parameter_combinations$vars[comb])), collapse = ""),
                         sep = "_")
    cl_list[[result_name]] <- result$cluster

    metrics_row <- result$time
    if (!tl_null) {
      metrics_row <- c(result$valid, metrics_row)
    }
    metrics_df <- rbind(metrics_df, metrics_row)
    rownames(metrics_df)[nrow(metrics_df)] <- result_name
  }

  if(tl_null) colnames(metrics_df) <- "Time"
  else colnames(metrics_df) <- c("Purity","Fmeasure","RI","Time")

  result_list <- list("cluster" = cl_list, "metrics" = metrics_df)

  return(result_list)
}

#' Perform spectral clustering for a given combination of indexes and kernel
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars Character vector of indexes names for clustering (e. g. vars1 <- c("dtaMEI", "dtaMHI"))
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_spc_aux <- function(ind_data, vars, kernel = "rbfdot", n_cluster = 2,
                             true_labels=NULL, ...){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(get(vars) %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  t0 <- Sys.time()

  # create target variable
  nrow_ind_data <- nrow(ind_data)
  specc_output <- kernlab::specc(as.matrix(ind_data[,get(vars)]), centers = n_cluster,  kernel = kernel, ...)
  clus <- specc_output@.Data

  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform spectral clustering for a different combinations of indexes and kernels
#'
#' @param ind_data Dataframe containing indexes applied to the main data and derivatives
#' @param vars_list List of characters representing variable sets for clustering
#' @param kernel_list List of kernels
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation (if it is not known true_labels is set to NULL)
#' @param ... Additional arguments (unused)
#'
#' @return A list containing kkmeans clustering results for each configuration
#' @export
clustInd_spc <- function(ind_data, vars_list, kernel_list = c("rbfdot", "polydot"),
                         n_cluster = 2, true_labels = NULL, ...) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check all elements in vars_list are characters
  if(!all(is.character(vars_list))){
    stop("Elemens in vars_list should be characters")
  }

  # Check if indices and distances lists are provided
  if (!is.character(vars_list) || !is.character(kernel_list) ||
      length(vars_list) == 0 || length(kernel_list) == 0) {
    stop("Invalid 'method_list' or 'method_list'. Both must be non-empty character vectors.")
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = vars_list, kernel = kernel_list)

  cl_list <- list()
  tl_null <- is.null(true_labels)
  nc <- ifelse(tl_null, 1, 4)
  metrics_df <- data.frame(matrix(ncol = nc, nrow = 0))

  for(comb in 1:nrow(parameter_combinations)){

    # Apply spc to each combination
    result <- clustInd_spc_aux(ind_data, as.character(parameter_combinations$vars[comb]),
                               as.character(parameter_combinations$kernel[comb]),
                               n_cluster, true_labels, ...)
    result_name <- paste("spc", parameter_combinations$kernel[comb],
                         paste(get(as.character(parameter_combinations$vars[comb])), collapse = ""), sep = "_")
    cl_list[[result_name]] <- result$cluster

    metrics_row <- result$time
    if (!tl_null) {
      metrics_row <- c(result$valid, metrics_row)
    }
    metrics_df <- rbind(metrics_df, metrics_row)
    rownames(metrics_df)[nrow(metrics_df)] <- result_name
  }

  if(tl_null) colnames(metrics_df) <- "Time"
  else colnames(metrics_df) <- c("Purity","Fmeasure","RI","Time")

  result_list <- list("cluster" = cl_list, "metrics" = metrics_df)
  return(result_list)
}



#' EHyClus method for clustering. It creates a multivariate dataset containing
#' the epigraph, hypograph and/or modified version on the curves and derivatives
#' and perform hierarchical clustering, kmeans, kernel kmeans, support vector
#' clustering and spectral clustering
#'
#' @param curves Dataset containing the curves to apply a clustering algorithm.
#' The functional dataset can be one dimensional (nxp) where n is the number of
#' curves and p the number of time points, or multidimensional (nxpxk) where k
#' represents the number of dimensions in the data
#' @param t Grid
#' @param nbasis Number of basis for the B-splines
#' @param norder Order of the B-splines
#' @param indices Names ofthe indices that need to be generated. They should be
#' one or more between EI, HI, MEI, MHI. Depending on the dimension on the data
#' they are calculated for one or multiple dimension
#' @param vars_list List of characters representing variable sets for clustering
#' @param l_method_hierarch List of clustering methods for hierarchical clustering
#' @param l_dist_hierarch List of distances for hierarchical clustering
#' @param l_dist_kmeans List of distances for kmeans clustering
#' @param l_kernel List of kernels
#' @param l_method_svc List of clustering methods for support vector clustering
#' @param n_clusters Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param ... Additional arguments (unused)
#'
#' @return A list containing the clustering partition for each method and indexes
#' combination and a data frame containing the time elapsed for obtaining a
#' clustering partition of the indexes dataset for each methodology
#' @export
#'
EHyClus <- function(curves, t, vars_list, nbasis = 30, norder = 4,
                    indices = c("EI", "HI", "MEI", "MHI"),
                    l_method_hierarch = c("single","complete","average","centroid","ward.D2"),
                    l_dist_hierarch = c("euclidean", "manhattan"),
                    l_dist_kmeans =  c("euclidean", "mahalanobis"),
                    l_kernel = c("rbfdot", "polydot"), l_method_svc = c("kmeans", "kernkmeans"),
                    n_clusters = 2, true_labels = NULL, ...){

  # Check all elements in vars_list are characters
  if(!all(is.character(vars_list))){
    stop("Elemens in vars_list should be characters")
  }

  # Check indices names
  if(!all(indices %in% c("EI", "HI", "MEI", "MHI"))){
    stop("Indices should be one or more of the following: EI, HI, MEI, MHI")
  }

  # Generate the dataset with the indexes
  ind_curves <- ind(curves, t, c(min(t), max(t)), nbasis, norder, indices)

    # Hierarchical clustering
  cl_hierarch <- clustInd_hierarch(ind_data = ind_curves, vars_list = vars_list,
                                   method_list = l_method_hierarch, dist_list = l_dist_hierarch,
                                   n_cluster = n_clusters, true_labels = true_labels)

  # kmeans
  cl_kmeans <- clustInd_kmeans(ind_data = ind_curves, vars_list = vars_list,
                               dist_list = l_dist_kmeans, n_cluster = n_clusters,
                               true_labels = true_labels)

  # kernel kmeans
  cl_kkmeans <- clustInd_kkmeans(ind_data = ind_curves, vars_list = vars_list,
                                 kernel_list = l_kernel, n_cluster = n_clusters,
                                 true_labels = true_labels)

  # support vector clustering
  cl_svc <- clustInd_svc(ind_data = ind_curves, vars_list = vars_list,
                         method_list = l_method_svc, n_cluster = n_clusters,
                         true_labels = true_labels)

  # spectral clustering
  cl_spc <- clustInd_spc(ind_data = ind_curves, vars_list = vars_list,
                         kernel_list = l_kernel, n_cluster = n_clusters,
                         true_labels = true_labels)

  cluster <- c(cl_hierarch$cluster, cl_kmeans$cluster, cl_kkmeans$cluster,
               cl_svc$cluster, cl_spc$cluster)

  metrics <- rbind(cl_hierarch$metrics, cl_kmeans$metrics, cl_kkmeans$metrics,
                   cl_svc$metrics, cl_spc$metrics)

  return(list("cluster" = cluster, "metrics" = metrics))
}
