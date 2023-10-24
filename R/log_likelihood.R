#' Log-likelihood for PK parameters
#'
#' Evaluates the prior distribution (normal) of the log parameter vector with
#' normally distributed error.
#'
#' @param lpr log-PK parameter vector with error: (lv_1, lk_10, lk_12, lk_21, ler_mean)
#' @param ivt List with containing start of infusion times (h), end of infusion times (h),
#' and rate of infusion (g/h) at each dose
#' @param dat Concentration data frame of the form: data.frame(time (h), concentration (milligrams/dl))
#' @param init Initial concentrations in each compartment
#'
#' @return The value of the log-likelihood evaluated at the specified log-parameter vector
#'
#' @export
#' @examples
#' ivt_d <- list(list(begin=0.0, end=0.5, k_R=6),
#'               list(begin=8.0, end=8.5, k_R=6),
#'               list(begin=16.0, end=16.5, k_R=6))
#' dat_d <- data.frame(time_h = c(1,4,40), conc_mcg_ml = c(82.7,80.4,60))
#'
#' log_likelihood(lpr = c(getOption("pkpredict.pip.default.prior")$log_pk_mean,
#'                        getOption("pkpredict.pip.default.prior")$log_err_mean),
#'                ivt = ivt_d, dat = dat_d)
#'

log_likelihood <- function(lpr, ivt, dat, init = c(0,0)) {
  epr <- exp(lpr)

  sol <- pk_solution(v_1 = epr[1], k_10 = epr[2],
                     k_12 = epr[3], k_21 = epr[4], ivt = ivt, init = init)

  colnames(dat) <- c("time_h", "conc_mcg_ml")
  dat$pconc_g_l   <- sol(dat$time_h)[1,]
  dat$pconc_mcg_ml <- 1000*dat$pconc_g_l

  res <- with(dat, sum(dnorm(conc_mcg_ml, pconc_mcg_ml,
                             # proportional error
                             pconc_mcg_ml * epr[5], log=TRUE)))
  attr(res, "soln") <- sol
  return(res)
}
