valid <- function(true_labels, clusters){

  if (is.integer(true_labels))
    true_labels = as.numeric(true_labels)
  if (is.integer(clusters))
    clusters = as.numeric(clusters)
  if (!is.vector(true_labels) || !is.numeric(true_labels))
    stop("true_labels should be a numeric vector")
  if (!is.vector(clusters) || !is.numeric(clusters))
    stop("clusters should be a numeric vector")
  if (length(true_labels) != length(clusters))
    stop("The length of the true_labels vector should be equal to the length of the clusters vector")

  tbl = table(clusters, true_labels) # contingency table
  conv_df = as.data.frame.matrix(tbl)

  # Purity
  res_purity = sum(apply(conv_df, 1, max))/length(true_labels)

  # True positives(tp), false positives(fp), true negatives(tn), and false negatives(fn)
  # (needed for calculating RI)
  tp_plus_fp = sum(choose(rowSums(conv_df), 2))
  tp_plus_fn = sum(choose(colSums(conv_df), 2))
  tp = sum(choose(as.vector(as.matrix(conv_df)), 2))
  fp = tp_plus_fp - tp
  fn = tp_plus_fn - tp
  tn = choose(sum(as.vector(as.matrix(conv_df))), 2) - tp - fp - fn

  # Precision and recall (needed for calculating Fmeasure)
  prec = tp / (tp + fp) # Precision
  rec = tp / (tp + fn) # Recall

  # (Purity, Fmeasure, RI)
  res <- c(round(res_purity, 4), round(2 * ((prec * rec) / (prec + rec)), 4),
           round((tp + tn) / (tp + fp + fn + tn), 4))
  names(res) <- c("Purity", "Fmeasure", "RI")

  return(as.table(res))
}
