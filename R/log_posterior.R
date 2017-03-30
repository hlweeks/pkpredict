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

log_posterior <- function(lpr, ivt, dat) {
  dat <- na.omit(dat)

  if(nrow(dat) < 1) {
    log_prior(lpr)
  } else {
    log_prior(lpr) + log_likelihood(lpr, ivt, dat)
  }
}
