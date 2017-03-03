#' Log-likelihood for PK parameters
#'
#' @param lpr log-PK parameter vector
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param dat Concentration data frame of the form: data.frame(time_h, conc_mg_dl)
#' @param ini Initial concentrations in each compartment
#'
#' @return The value of the log-likelihood evaluated at the specified log-parameter vector
#'
#' @export
#'
#' @examples

log_likelihood <- function(lpr, ivt, dat, ini=c(0,0)) {
  epr <- exp(lpr)
  sol <- pk_solution(v_1=epr[1], k_10=epr[2], k_12=epr[3], k_21=epr[4], ivt=ivt)
  dat$pconc_g_l   <- sol(dat$time_h)[1,]
  dat$pconc_mg_dl <- 1000*dat$pconc_g_l
  with(dat, sum(dnorm(conc_mg_dl, pconc_mg_dl, epr[5], log=TRUE)))
}
