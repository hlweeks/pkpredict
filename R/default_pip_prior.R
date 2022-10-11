#' Default prior specification for piperacillin 
#'
#' The prior values stored here were derived from random effects estimates from a prior pharmacokinetic
#' study using Gaussian kernel density estimation. See vignette for details. This prior specification is 
#' also stored as an option when the package is loaded and can be retreived using 
#' `getOption("pkpredict.pip.default.prior")`. 
#'
#' @format A list with 4 elements specifying the default prior distribution for function arguments: 
#' log_pk_mean, log_pk_vcov, log_err_mean, log_err_sdev
#'
#' @keywords datasets
#'
#' @examples
#' data(default_pip_prior)
"default_pip_prior"
