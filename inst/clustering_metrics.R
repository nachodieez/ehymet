n <- 30
true_labels <- c(rep(1, n), rep(2, n))
global_metrics <- list()

# ex1
for (i_sim in seq_len(8)) {
  for (n_clusters in c(2, 3)) {
    metrics <- list()
    for (iter in seq_len(10)) {
      curves <- sim_model_ex1(n = n, i_sim = i_sim)
      res <- EHyClus(
        curves = curves,
        vars_combinations = generic_vars_combinations(multidimensional = FALSE),
        true_labels = true_labels,
        n_clusters = n_clusters
      )
      metrics[[length(metrics) + 1]] <- res$metrics
    }
    metrics <- do.call(rbind, metrics)
    metrics <- metrics[order(metrics[, "ARI"], decreasing = TRUE), ]

    write.csv(
      x = metrics,
      file = paste0("inst/metrics_csv/ex1_sim", i_sim, "_nclusters", n_clusters, ".csv")
    )

    write.csv(
      x = data.frame(ARI = cor(metrics, use = "pairwise.complete.obs")["ARI", ]),
      file = paste0("inst/metrics_csv/cor_ex1_sim", i_sim, "_nclusters", n_clusters, ".csv")
    )

    global_metrics[[length(global_metrics) + 1]] <- metrics
  }
}

# ex2 --> i_sim 1 and i_sim 2
for (i_sim in seq_len(2)) {
  for (n_clusters in c(2, 3)) {
    metrics <- list()
    for (iter in seq_len(10)) {
      curves <- sim_model_ex2(n = n, i_sim = i_sim)
      res <- EHyClus(
        curves = curves,
        vars_combinations = generic_vars_combinations(multidimensional = FALSE),
        true_labels = true_labels,
        n_clusters = n_clusters
      )
      metrics[[length(metrics) + 1]] <- res$metrics
    }
    metrics <- do.call(rbind, metrics)
    metrics <- metrics[order(metrics[, "ARI"], decreasing = TRUE), ]

    write.csv(
      x = metrics,
      file = paste0("inst/metrics_csv/ex2_sim", i_sim, "_nclusters", n_clusters, ".csv")
    )

    write.csv(
      x = data.frame(ARI = cor(metrics, use = "pairwise.complete.obs")["ARI", ]),
      file = paste0("inst/metrics_csv/cor_ex2_sim", i_sim, "_nclusters", n_clusters, ".csv")
    )

    global_metrics[[length(global_metrics) + 1]] <- metrics
  }
}

# ex2 --> i_sim 3 and i_sim 4
for (i_sim in 3:4) {
  for (n_clusters in c(2, 3)) {
    metrics <- list()
    for (iter in seq_len(10)) {
      curves <- sim_model_ex2(n = n, i_sim = i_sim)
      res <- EHyClus(
        curves = curves,
        vars_combinations = generic_vars_combinations(multidimensional = TRUE),
        true_labels = true_labels,
        n_clusters = n_clusters
      )
      metrics[[length(metrics) + 1]] <- res$metrics
    }
    metrics <- do.call(rbind, metrics)
    metrics <- metrics[order(metrics[, "ARI"], decreasing = TRUE), ]

    write.csv(
      x = metrics,
      file = paste0("inst/metrics_csv/ex2_sim", i_sim, "_nclusters", n_clusters, ".csv")
    )

    write.csv(
      x = data.frame(ARI = cor(metrics, use = "pairwise.complete.obs")["ARI", ]),
      file = paste0("inst/metrics_csv/cor_ex2_sim", i_sim, "_nclusters", n_clusters, ".csv")
    )

    global_metrics[[length(global_metrics) + 1]] <- metrics
  }
}

# writing global metrics
global_metrics <- do.call(rbind, global_metrics)
global_metrics <- global_metrics[order(global_metrics[, "ARI"], decreasing = TRUE), ]

write.csv(
  x = global_metrics,
  file = "inst/metrics_csv/global_metrics.csv"
)

write.csv(
  x = data.frame(ARI = cor(global_metrics, use = "pairwise.complete.obs")["ARI", ]),
  file = "inst/metrics_csv/cor_global_metrics.csv"
)

# separating by clusters
files <- list.files(path = "inst/metrics_csv")
two_clusters_files   <- files[grep("^ex.*nclusters2", files)]
three_clusters_files <- files[grep("^ex.*nclusters3", files)]

two_clusters   <- do.call(rbind, lapply(paste0("inst/metrics_csv/", two_clusters_files), read.csv, row.names = 1))
three_clusters <- do.call(rbind, lapply(paste0("inst/metrics_csv/", three_clusters_files), read.csv, row.names = 1))

write.csv(
  x = data.frame(ARI = cor(two_clusters, use = "pairwise.complete.obs")["ARI", ]),
  file = "inst/metrics_csv/cor_two_clusters.csv"
)

write.csv(
  x = data.frame(ARI = cor(three_clusters, use = "pairwise.complete.obs")["ARI", ]),
  file = "inst/metrics_csv/cor_three_clusters.csv"
)

## combined ##

true_labels <- c(rep(1, 30), rep(2, 30), rep(3, 15), rep(4, 15))
global_metrics <- list()


metrics <- list()
for (iter in seq_len(10)) {
  curves <- rbind(
    sim_model_ex1(n = 30, i_sim = 1),
    sim_model_ex1(n = 15, i_sim = 7)
  )
  res <- EHyClus(
    curves = curves,
    vars_combinations = generic_vars_combinations(multidimensional = FALSE),
    true_labels = true_labels,
    n_clusters = 4
  )
  metrics[[length(metrics) + 1]] <- res$metrics
}
metrics <- do.call(rbind, metrics)
metrics <- metrics[order(metrics[, "ARI"], decreasing = TRUE), ]

write.csv(
  x = metrics,
  file = "inst/metrics_csv/4clusters.csv"
)

write.csv(
  x = data.frame(ARI = cor(metrics, use = "pairwise.complete.obs")["ARI", ]),
  file = "inst/metrics_csv/cor_4clusters.csv"
)
