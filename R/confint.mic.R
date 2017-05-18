#' Confidence intervals for ftMIC statistic
#'
#' @param ftmic "mic" object obtained from \code{mic_stat}
#' @param conf.level Desired confidence level of the interval
#'
#' @return CI
#' @export
#'
#' @examples
#'

confint.mic <- function(ftmic, conf.level = .95){
  if(abs(conf.level) > 1){stop("conf.level must be between 0 and 1")}
  alp <- 1 - conf.level

  est <- ftmic$est

  # SE of logit(statistic)
  grd_mic <- fdGrad(est$par, function(pars) {
    #try extracting the est object and going direct to the time>mic computation?
    #the optim function might be what's making this take a few seconds to run
    #the optim function is the difference between bayes.R and this
    mic <- mic_stat(pars, ivt = ftmic$ivt, dat = ftmic$dat,
                    times = ftmic$tms, con = ftmic$con, th = ftmic$th)$ftmic
    log(mic/(1-mic)) ## constrain between 0 and 1
  })

  sde_mic <- sqrt(diag(t(grd_mic) %*% solve(-est$hessian) %*% grd_mic))

  # Get CI for logit transformed statistic then backtransform to original scale
  ci_logit_mic <- log(ftmic$ftmic/(1-ftmic$ftmic)) + c(-1,1)*qnorm(1-alp/2)*sde_mic
  ci_mic <- exp(ci_logit_mic)/(1 + exp(ci_logit_mic))

  res <- list("ftmic" = ftmic, "95% CI lower" = ci_mic[1], "95% CI upper" = ci_mic[2])

  return(res)
}
