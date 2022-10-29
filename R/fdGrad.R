#' Finite-difference gradient
#'
#' @param pars Parameters involved in gradient computation
#' @param fun Function for which the finite difference gradient is computed
#' @param ... Other arguments for \code{fun} not involved in gradient computation
#' @param .relStep Amount to shift parameter values when calculating gradient
#' @param minAbsPar Minimum absolute parameter value
#'
#' @return Finite difference gradient
#'
#' @export
#'
fdGrad <- function (pars, fun, ..., .relStep = (.Machine$double.eps)^(1/2), minAbsPar = 0){
  npar <- length(pars)
  incr <- ifelse(abs(pars) <= minAbsPar, .relStep, (abs(pars) -
                                                      minAbsPar) * .relStep)
  sapply(1:npar, function(i) {
    del <- rep(0, npar)
    del[i] <- incr[i]


    # Save two sides of difference below - if both inf (same sign), return 0
    diff1 <- do.call(fun, list(pars + del, ...))
    diff2 <- do.call(fun, list(pars - del, ...))

    if(is.infinite(diff1) & is.infinite(diff2)){
      if(sign(diff1) == sign(diff2)){
        0
      }else{
        Inf
      }
    }else{
      (diff1 - diff2)/incr[i]/2
    }
  })
}
