library(clusterCrit)
library(tidyverse)

n <- 20
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
# rownames(all_metrics) <- sapply(generic_vars_combinations(multidimensional = TRUE),
#                                 \(x) paste0(x, collapse = "_"))

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
  arrange(desc(ARI)) |>
  View()

data.frame(ARI = cor(all_metrics)["ARI", ])

# In clustering analysis, the trace_WIB can be used to compare different clustering solutions.
# A clustering solution with a lower trace_WIB is generally preferred as it indicates more compact clusters.

## univariate ##

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


curves <- sim_model_ex2(n = n, i_sim = 1)
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

curves <- sim_model_ex2(n = n, i_sim = 2)
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

all_metrics2 |>
  arrange(desc(ARI)) |>
  View()

data.frame(ARI = cor(all_metrics)["ARI", ])  # mult
data.frame(ARI = cor(all_metrics2)["ARI", ]) # uni

data.frame(ARI = cor(rbind(all_metrics, all_metrics2))["ARI", ])

rbind(all_metrics, all_metrics2) |>
  select(ARI, trace_wib) |>
  arrange(desc(ARI)) |>
  View()
