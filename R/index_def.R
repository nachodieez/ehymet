EI <- function(curves){
  n_curves <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  index <- apply(curves,1, function(y) sum(apply(curves,1,function(x) sum(x>=y)==lengthcurves)))/n_curves
  return (1-index)
}
