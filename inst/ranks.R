set.seed(42)

library(tidyverse)

f <- function(n = 30, iters = 5, i_sim = 1, n_clusters = 2) {
  true_labels <- c(rep(1, n), rep(2, n))
  metrics <- list()
  for (iter in seq_len(iters)) {
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
  metrics <- metrics[,-ncol(metrics)]

  ranks <- rbind.data.frame(lapply(metrics, \(x) rank(-x)))
  rownames(ranks) <- rownames(metrics)
  colnames(ranks) <- paste0("rank_", colnames(ranks))
  ranks[["combined_rank"]] <- rank(rowSums(ranks[,-c(1,2,3,4)]))
  ranks[["combined_rank2"]] <- rank(rowSums(ranks[,-c(1,2,3,4,5,6)]))

  df <- merge(metrics, ranks, by = "row.names")
  rownames(df) <- df$Row.names
  df[,1] <- NULL
  as_tibble(df) |>
    arrange(desc(ARI))
}

df <- f(i_sim = 8)
cor(df)[,"ARI"]

# for_comparison <- df[,c("ARI", "rank_ARI", "combined_rank", "combined_rank2")]

compare_with_ranks_1 <- function(df, trials = 10000) {
  successes <- 0
  for (trial in seq_len(trials)) {
    ranks <- df$combined_rank
    ARIs  <- df$ARI

    # Calculate deciles for the ranks
    rank_deciles <- quantile(ranks, probs = seq(0, 1, 0.1))

    # Extract indices for each decile
    second_decile_indices <- which(ranks <= rank_deciles[3])
    fifth_decile_indices  <- which(ranks > rank_deciles[5] & ranks <= rank_deciles[6])
    ninth_decile_indices  <- which(ranks > rank_deciles[9])

    # Randomly sample one value from each decile
    second_decile_value <- sample(ARIs[second_decile_indices], 1)
    fifth_decile_value  <- sample(ARIs[fifth_decile_indices], 1)
    ninth_decile_value  <- sample(ARIs[ninth_decile_indices], 1)

    if (second_decile_value > fifth_decile_value && fifth_decile_value > ninth_decile_value) {
      successes <- successes + 1
    }
  }

  successes / trials
}

compare_with_ranks_2 <- function(df, trials = 10000) {
  successes <- 0
  for (trial in seq_len(trials)) {
    ranks <- df$combined_rank2
    ARIs  <- df$ARI

    # Calculate deciles for the ranks
    rank_deciles <- quantile(ranks, probs = seq(0, 1, 0.1))

    # Extract indices for each decile
    second_decile_indices <- which(ranks <= rank_deciles[3])
    fifth_decile_indices  <- which(ranks > rank_deciles[5] & ranks <= rank_deciles[6])
    ninth_decile_indices  <- which(ranks > rank_deciles[9])

    # Randomly sample one value from each decile
    second_decile_value <- sample(ARIs[second_decile_indices], 1)
    fifth_decile_value  <- sample(ARIs[fifth_decile_indices], 1)
    ninth_decile_value  <- sample(ARIs[ninth_decile_indices], 1)

    if (second_decile_value > fifth_decile_value && fifth_decile_value > ninth_decile_value) {
      successes <- successes + 1
    }
  }

  successes / trials
}


for (i_sim in 1:8) {
  df <- f(i_sim = i_sim, n = 50, iters = 10)
  successes_1 <- df |>
    compare_with_ranks_1()

  successes_2 <- df |>
    compare_with_ranks_2()

  cat(paste0("i_sim = ", i_sim,
             "\nMax ARI: ", df$ARI[1],
             "\nSuccesses (with silhouette and dunn): ", successes_1,
             "\nSuccesses (w/o silhouette and dunn): ", successes_2, "\n"))
}
