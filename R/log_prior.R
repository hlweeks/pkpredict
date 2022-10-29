#' Log-prior for pharmacokinetic parameters
#'
#' Evaluates the prior distribution (normal) of the log parameter vector with
#' normally distributed error.
#'
#' Parameter distributions:
#' log(lv_1, lk_10, lk_12, lk_21) ~ N(mu, sig)
#' log(error) ~ N(ler_mean, ler_sdev)
#'
#'
#' @param lpr log-PK parameter vector with error: (lv_1, lk_10, lk_12, lk_21, ler_mean)
#' @param mu log-mean of the PK parameter distribution
#' @param sig log-variance-covariance matrix of the PK parameter distribution
#' @param ler_mean log-mean of the error distribution
#' @param ler_sdev log-standard deviation of the error distribution
#'
#' @import mvtnorm
#'
#' @return The value of the log-prior evaluated at the specified log-parameter vector
#'
#' @export
#'
#' @examples
#' log_prior(lpr = c(getOption("pkpredict.pip.default.prior")$log_pk_mean,
#'                   getOption("pkpredict.pip.default.prior")$log_err_mean))
#'


log_prior <- function(lpr,
                      mu = getOption("pkpredict.pip.default.prior")$log_pk_mean,
                      sig = getOption("pkpredict.pip.default.prior")$log_pk_vcov,
                      ler_mean = getOption("pkpredict.pip.default.prior")$log_err_mean,
                      ler_sdev = getOption("pkpredict.pip.default.prior")$log_err_sd){

  val <- mvtnorm::dmvnorm(lpr[1:4], mu, sig, log = TRUE) +
           dnorm(lpr[5], mean=ler_mean, sd=ler_sdev, log=TRUE)

  names(val) <- NULL
  return(val)
}
