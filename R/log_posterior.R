#' Log-posterior for PK parameters
#'
#' Evaluates the prior distribution (normal) of the log parameter vector with
#' normally distributed error.
#'
#' @param lpr log-PK parameter vector with error: (lv_1, lk_10, lk_12, lk_21, ler_mean)
#' @param ivt List with containing start of infusion times (h), end of infusion times (h),
#' and rate of infusion (g/h) at each dose
#' @param dat Concentration data frame of the form: data.frame(time_h, conc_mg_dl)
#' @param mu prior pk param mean
#' @param sig prior pk vcov matrix
#' @param ler_mean prior error mean
#' @param ler_sdev prior error sd
#'
#' @return The log-posterior distribution evaluated at the specified log-parameter vector
#'
#' @examples
#' ivt_d <- list(list(begin=0.0, end=0.5, k_R=6),
#'               list(begin=8.0, end=8.5, k_R=6),
#'               list(begin=16.0, end=16.5, k_R=6))
#' dat_d <- data.frame(time_h = c(1,4,40), conc_mcg_ml = c(82.7,80.4,60))
#'
#' log_posterior(lpr = c(getOption("pkpredict.pip.default.prior")$log_pk_mean,
#'                       getOption("pkpredict.pip.default.prior")$log_err_mean),
#'               ivt = ivt_d, dat = dat_d)
#'
#' @export
#'

log_posterior <- function(lpr, ivt, dat,
                          mu = getOption("pkpredict.pip.default.prior")$log_pk_mean,
                          sig = getOption("pkpredict.pip.default.prior")$log_pk_vcov,
                          ler_mean = getOption("pkpredict.pip.default.prior")$log_err_mean,
                          ler_sdev = getOption("pkpredict.pip.default.prior")$log_err_sd){
  if(is.null(mu) | is.null(sig) | is.null(ler_mean) | is.null(ler_sdev)){
    stop("One or more of `mu`, `sig`, `ler_mean`, or `ler_sdev` arguments is NULL. Set prior arguments manually or use `options(pkpredict.pip.default.prior = data(default_pip_prior))`")
  }

  dat <- na.omit(dat)

  if(nrow(dat) < 1) {
    log_prior(lpr, mu=mu, sig=sig, ler_mean=ler_mean, ler_sdev=ler_sdev)
  } else {
    log_like <- log_likelihood(lpr, ivt, dat)
    res <- log_prior(lpr, mu, sig, ler_mean, ler_sdev) + log_like
    attr(res, "soln") <- attributes(log_like)$soln
    res
  }
}
