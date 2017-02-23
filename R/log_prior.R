#' Log-prior for pharmacokinetic parameters
#'
#' @param lpr log-PK parameter vector
#' @param mu log-mean of the PK parameter distribution
#' @param sig log-variance covariance matrix of the PK parameter distribution
#' @param ler_mean log-mean of the error distribution
#' @param ler_sdev log-standard deviation of the error distribution
#'
#' @return Value of the log-prior for the given \code{lpr} vector
#'
#' @examples
#' lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7)
#' lpk_vcov_d <- 300 * matrix(c(  0.00167, -0.00128,      0,      0,
#'                                -0.00128,  0.00154,      0,      0,
#'                                0,        0, .00015,      0,
#'                                0,        0,      0, .00015), 4, 4)
#' ler_mean_d <- 2.33
#' ler_sdev_d <- 0.32
#' lpr_d  <- c(lpk_mean_d, ler_mean_d)
#' log_prior(lpr_d, lpk_mean_d, lpk_vcov_d, ler_mean_d, ler_sdev_d)
#'
#' @export


log_prior <- function(lpr, mu, sig, ler_mean, ler_sdev)
  dmvnorm(lpr[1:4], mu, sig, log = TRUE) +
  dnorm(lpr[5], mean=ler_mean, sd=ler_sdev, log=TRUE)
