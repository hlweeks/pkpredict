#' Log-posterior for PK parameters
#'
#' Evaluates the prior distribution (normal) of the log parameter vector with
#' normally distributed error.
#'
#' @param lpr log-PK parameter vector with error: (lv_1, lk_10, lk_12, lk_21, ler_mean)
#' @param ivt List with containing start of infusion times (h), end of infusion times (h),
#' and rate of infusion (g/h) at each dose
#' @param dat Concentration data frame of the form: data.frame(time_h, conc_mg_dl)
#'
#' @return The log-posterior distribution evaluated at the specified log-parameter vector
#'
#' @export
#'

log_posterior <- function(lpr, ivt, dat,
                          mu = c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7),
                          sig = 300 * matrix(c(   0.00167,  -0.00128,      0,      0,
                                                  -0.00128,   0.00154,      0,      0,
                                                  0,         0, .00015,      0,
                                                  0,         0,      0, .00015), 4, 4),
                          ler_mean = 2.33, ler_sdev = 0.32) {
  dat <- na.omit(dat)

  if(nrow(dat) < 1) {
    log_prior(lpr)
  } else {
    log_like <- log_likelihood(lpr, ivt, dat)
    res <- log_prior(lpr, mu=mu, sig=sig, ler_mean=ler_mean, ler_sdev=ler_sdev) + log_like
    attr(res, "soln") <- attributes(log_like)$soln
    res
  }
}
