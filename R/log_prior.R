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
#' # Default parameter values based on prior study
#' lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7)
#' lpk_vcov_d <- 300 * matrix(c(   0.00167,  -0.00128,      0,      0,
#'                                -0.00128,   0.00154,      0,      0,
#'                                       0,         0, .00015,      0,
#'                                       0,         0,      0, .00015), 4, 4)
#' ler_mean_d <- 2.33
#' ler_sdev_d <- 0.32
#'
#' lpr_d  <- c(lpk_mean_d, ler_mean_d)
#' log_prior(lpr = lpr_d, mu = lpk_mean_d, sig = lpk_vcov_d,
#'           ler_mean = ler_mean_d, ler_sdev = ler_sdev_d)


log_prior <- function(lpr, mu, sig,
                      ler_mean, ler_sdev){

  val <- mvtnorm::dmvnorm(lpr[1:4], mu, sig, log = TRUE) +
           dnorm(lpr[5], mean=ler_mean, sd=ler_sdev, log=TRUE)

  names(val) <- NULL
  return(val)
}
