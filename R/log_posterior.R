#' Log posterior for PK parameters
#'
#' @param lpr log-PK parameter vector
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param dat Concentration data frame of the form: data.frame(time_h, conc_mg_dl)
#'
#' @return The log-posterior distribution evaluated at the specified log-parameter vector
#'
#' @export
#'
#' @examples

log_posterior <- function(lpr, ivt, dat) {
  dat <- na.omit(dat)
  if(nrow(dat) < 1) {
    log_prior(lpr)
  } else {
    log_prior(lpr) + log_likelihood(lpr, ivt, dat)
  }
}
