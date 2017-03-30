#' Confidence intervals for ftMIC statistic
#'
#' @param ftmic "mic" object obtained from \code{mic_stat}
#' @param est Output from \code{optim}
#' @param conf.level Desired confidence level of the interval
#'
#' @return CI
#' @export
#'
#' @examples
#'

confint.mic <- function(ftmic, est, conf.level = .95){
  alp <- 1 - conf.level

  # SE of logit(statistic)
  grd_mic <- fdGrad(est$par, function(pars) {
    mic <- mic_stat(pk_pars = pars, ivt, tms, con, th = thres)
    log(mic/(1-mic)) ## constrain between 0 and 1
  })
  sde_mic <- sqrt(diag(t(grd_mic) %*% solve(-est$hessian) %*% grd_mic))

  # Get CI for logit transformed statistic then backtransform to original scale
  ci_logit_mic <- log(frac_mic/(1-frac_mic)) + c(-1,1)*qnorm(1-alp/2)*sde_mic
  ci_mic <- exp(ci_logit_mic)/(1 + exp(ci_logit_mic))

  res <- list("ftmic" = ftmic, "95% CI lower" = ci_mic[1], "95% CI upper" = ci_mic[2])

  return(res)
}
