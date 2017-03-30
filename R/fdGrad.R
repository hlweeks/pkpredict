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
#' @examples
fdGrad <- function (pars, fun, ...,
                    .relStep = (.Machine$double.eps)^(1/2),
                    minAbsPar = 0) {

  npar <- length(pars)
  incr <- ifelse(abs(pars) <= minAbsPar, .relStep,
                 (abs(pars)-minAbsPar) * .relStep)
  ival <- do.call(fun, list(pars, ...))
  diff <- rep(0,npar)
  sapply(1:npar, function(i) {
    del <- rep(0,npar)
    del[i] <- incr[i]
    (do.call(fun, list(pars+del, ...))-ival)/incr[i]
  })
}
