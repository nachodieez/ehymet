#' Epigraph Index (EI) for a univariate functional dataset.
#'
#' The Epigraph Index of a curve x is one minus the proportion of curves in the sample
#' that are above x.
#'
#' @param curves A matrix where each row represents a curve, and each column
#' represents values along the curve.
#'
#' @return Return a numeric vector containing the EI for each curve
#' @export
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7),ncol = 3, nrow = 4)
#' EI(x)
EI <- function(curves){
  n_curves <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  index <- apply(curves,1, function(y) sum(apply(curves,1,function(x) sum(x>=y)==lengthcurves)))/n_curves
  return (1-index)
}
