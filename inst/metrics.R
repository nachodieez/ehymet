library(clusterCrit)
library(tidyverse)

infomax_metric <- function(cluster_labels) {
  # Calculate the frequency of each cluster
  cluster_freq <- table(cluster_labels)

  # Total number of data points
  n <- length(cluster_labels)

  # Calculate the probabilities of each cluster
  p <- cluster_freq / n

  # Calculate the entropy
  infomax <- -sum(p * log2(p))

  list(infomax = infomax)
}

get_metrics_unidimensional <- function(n = 50) {
  metrics <- data.frame()
  true_labels <- c(rep(1, n), rep(2, n))
  for (i_sim in 1:8) {
    curves <- sim_model_ex1(n = n, i_sim = i_sim)
    for (combination in generic_vars_combinations(multidimensional = FALSE)) {
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

metrics_unidimensional   <- get_metrics_unidimensional()
metrics_multidimensional <- get_metrics_multidimensional()

data.frame(RI = cor(metrics_unidimensional)["RI", ])   # uni
data.frame(RI = cor(metrics_multidimensional)["RI", ]) # mult

linear_correlation <- data.frame(RI = cor(rbind(metrics_unidimensional, metrics_multidimensional))["RI", ])

write.csv(metrics_unidimensional,   "inst/metrics_csv/metrics_unidimensional.csv")
write.csv(metrics_multidimensional, "inst/metrics_csv/metrics_multidimensional.csv")
write.csv(linear_correlation,       "inst/metrics_csv/linear_correlation.csv")




