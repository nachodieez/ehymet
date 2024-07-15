library(clusterCrit)
library(tidyverse)

# define the infomax metric for clustering
infomax_metric <- function(cluster_labels) {
  # Calculate the frequency of each cluster
  cluster_freq <- table(cluster_labels)

  # Calculate the total number of data points
  n <- length(cluster_labels)

  # Calculate the probabilities of each cluster
  p <- cluster_freq / n

  # Calculate the entropy
  infomax <- -sum(p * log2(p))

  list(infomax = infomax)
}

n <- 50
curves <- sim_model_ex2(n = n, i_sim = 4)

true_labels <- c(rep(1, n), rep(2, n))
all_metrics <- data.frame()
for (combination in generic_vars_combinations(multidimensional = TRUE)) {
  res <- EHyClus(curves, vars_combinations = list(combination),
                 true_labels = true_labels, only_best = TRUE)

  ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
  all_metrics <- rbind.data.frame(all_metrics, append(
    res$cluster[[1]]$valid,
    intCriteria(ind_curves, res$cluster[[1]]$cluster, "all")
  ))
}

curves <- sim_model_ex2(n = n, i_sim = 3)
for (combination in generic_vars_combinations(multidimensional = TRUE)) {
  res <- EHyClus(curves, vars_combinations = list(combination),
                 true_labels = true_labels, only_best = TRUE)

  ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
  all_metrics <- rbind.data.frame(all_metrics, append(
    res$cluster[[1]]$valid,
    intCriteria(ind_curves, res$cluster[[1]]$cluster, "all")
  ))
}

all_metrics |>
  arrange(desc(RI)) |>
  View()

data.frame(RI = cor(all_metrics)["RI", ])

# In clustering analysis, the trace_WIB can be used to compare different clustering solutions.
# A clustering solution with a lower trace_WIB is generally preferred as it indicates more compact clusters.

## univRIate ##

n <- 50
curves <- sim_model_ex1(n = n)

true_labels <- c(rep(1, n), rep(2, n))

all_metrics2 <- data.frame()
valid_combinations <- list()
for (combination in generic_vars_combinations(multidimensional = FALSE)) {
  tryCatch({
    res <- EHyClus(curves, vars_combinations = list(combination),
                   true_labels = true_labels, only_best = TRUE)

    ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
    all_metrics2 <- rbind.data.frame(all_metrics2, append(
      res$cluster[[1]]$valid,
      intCriteria(ind_curves, res$cluster[[1]]$cluster, "all")
    ))
    valid_combinations <- append(valid_combinations, list(combination))

  },
  error = function(cond) NA
  )
}

curves <- sim_model_ex1(n = n, i_sim = 5)
for (combination in generic_vars_combinations(multidimensional = FALSE)) {
  tryCatch({
    res <- EHyClus(curves, vars_combinations = list(combination),
                   true_labels = true_labels, only_best = TRUE)

    ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
    all_metrics2 <- rbind.data.frame(all_metrics2, append(
      res$cluster[[1]]$valid,
      intCriteria(ind_curves, res$cluster[[1]]$cluster, "all")
    ))
    valid_combinations <- append(valid_combinations, list(combination))

  },
  error = function(cond) NA
  )
}

curves <- sim_model_ex1(n = n, i_sim = 7)
for (combination in generic_vars_combinations(multidimensional = FALSE)) {
  tryCatch({
    res <- EHyClus(curves, vars_combinations = list(combination),
                   true_labels = true_labels, only_best = TRUE)

    ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
    all_metrics2 <- rbind.data.frame(all_metrics2, append(
      res$cluster[[1]]$valid,
      intCriteria(ind_curves, res$cluster[[1]]$cluster, "all")
    ))
    valid_combinations <- append(valid_combinations, list(combination))

  },
  error = function(cond) NA
  )
}

get_metrics_multidimensional <- function(n = 50) {
  metrics <- data.frame()
  true_labels <- c(rep(1, n), rep(2, n))
  for (i_sim in c(3, 4)) {
    curves <- sim_model_ex2(n = n, i_sim = i_sim)
    for (combination in generic_vars_combinations(multidimensional = TRUE)) {
      tryCatch({
        res <- EHyClus(curves, vars_combinations = list(combination),
                       true_labels = true_labels, only_best = TRUE)

        ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
        metrics <- rbind.data.frame(metrics, append(
          res$cluster[[1]]$valid,
          c(intCriteria(ind_curves, res$cluster[[1]]$cluster, "all"), infomax_metric(res$cluster[[1]]$cluster))
        ))
      },
      error = function(cond) NA
      )
    }
  }

  metrics |>
    arrange(desc(RI))
}

## multidimensional
all_metrics2 <- data.frame()
valid_combinations <- list()
curves <- sim_model_ex2(n = n, i_sim = 3)
for (combination in generic_vars_combinations(multidimensional = TRUE)) {
  tryCatch({
    res <- EHyClus(curves, vars_combinations = list(combination),
                   true_labels = true_labels, only_best = TRUE)

    ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
    all_metrics2 <- rbind.data.frame(all_metrics2, append(
      res$cluster[[1]]$valid,
      c(intCriteria(ind_curves, res$cluster[[1]]$cluster, "all"), infomax_metric(res$cluster[[1]]$cluster))
    ))
    valid_combinations <- append(valid_combinations, list(combination))

  },
  error = function(cond) NA
  )
}

curves <- sim_model_ex2(n = n, i_sim = 4)
for (combination in generic_vars_combinations(multidimensional = TRUE)) {
  tryCatch({
    res <- EHyClus(curves, vars_combinations = list(combination),
                   true_labels = true_labels, only_best = TRUE)

    ind_curves <- as.matrix(attr(res, "ind_curves")[combination])
    all_metrics2 <- rbind.data.frame(all_metrics2, append(
      res$cluster[[1]]$valid,
      intCriteria(ind_curves, res$cluster[[1]]$cluster, "all")
    ))
    valid_combinations <- append(valid_combinations, list(combination))

  },
  error = function(cond) NA
  )
}

all_metrics2 |>
  select(RI, davies_bouldin, dunn, silhouette) |>
  arrange(desc(RI)) |>
  View()

all_metrics2 |>
  arrange(desc(RI)) |>
  View()

data.frame(RI = cor(all_metrics)["RI", ])  # mult
data.frame(RI = cor(all_metrics2)["RI", ]) # uni

data.frame(RI = cor(rbind(all_metrics, all_metrics2))["RI", ])

rbind(all_metrics, all_metrics2) |>
  select(RI, trace_wib) |>
  arrange(desc(RI)) |>
  View()
