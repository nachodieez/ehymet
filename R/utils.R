#' Title
#'
#' @param curves
#' @param t
#' @param rangeval
#' @param nbasis
#' @param norder
#' @param ...
#'
#' @return
#' @noRd
#'
#' @examples
funspline <- function(curves, t, rangeval, nbasis, norder, ...){
  #Create B-spline basis
  basisobj <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis, norder = norder, ...)

  curves_dim <- length(dim(curves))

  if(curves_dim == 2){

    # Smooth data using B-spline basis
    ys <-  fda::smooth.basis(argvals = t, y = t(curves), fdParobj = basisobj)

    # Evaluate smoothed data and derivatives
    smooth <- t(fda::eval.fd(t,ys$fd,0)) # smoothed data
    deriv <- t(fda::eval.fd(t,ys$fd,1)) # first derivatives
    deriv2 <- t(fda::eval.fd(t,ys$fd,2)) # second derivatives
  } else if(curves_dim == 3){
    n_curves <- dim(curves)[1]
    l_curves <- dim(curves)[2]
    d_curves <- dim(curves)[3]

    # Initialize empty dataframes to store the results
    smooth <- array(rep(NaN,n_curves*l_curves),dim=c(n_curves,l_curves,d_curves))
    deriv <- array(rep(NaN,n_curves*l_curves),dim=c(n_curves,l_curves,d_curves))
    deriv2 <- array(rep(NaN,n_curves*l_curves),dim=c(n_curves,l_curves,d_curves))

    for(d in 1:dim(curves)[3]){
      # Smooth data using B-spline basis
      ys <-  fda::smooth.basis(argvals = t, y = t(curves[,,d]), fdParobj = basisobj)

      # Evaluate smoothed data and derivatives
      smooth[,,d] <- t(fda::eval.fd(t,ys$fd,0)) # smoothed data
      deriv[,,d] <- t(fda::eval.fd(t,ys$fd,1)) # first derivatives
      deriv2[,,d] <- t(fda::eval.fd(t,ys$fd,2)) # second derivatives
    }
  } else{
    stop("Invalid number of dimensions")
  }


  # Return a list containing the data and derivatives
  res <- list(
    "smooth" = smooth,
    "deriv" = deriv,
    "deriv2" = deriv2
    )

  return(res)
}
