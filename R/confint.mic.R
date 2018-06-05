#' Confidence intervals for ftMIC statistic
#'
#' @param object "mic" object obtained from \code{mic_stat}
#' @param parm which parameters should be given confidence intervals
#' @param level Desired confidence level of the interval
#' @param ... other arguments
#'
#' @return CI
#' @export
#'

confint.mic <- function(object, parm, level, ...){
  if(abs(level) > 1){stop("conf.level must be between 0 and 1")}
  alp <- 1 - level

  est <- object$ftmic

  # SE of logit(statistic)
  grd_mic <- fdGrad(ftobjectmic$pars, function(pars) {
    #try extracting the est object and going direct to the time>mic computation?
    #the optim function might be what's making this take a few seconds to run
    #the optim function is the difference between bayes.R and this
    mic <- mic_stat(pars, ivt = object$ivt, th = object$th)$ftmic
    log(mic/(1-mic)) ## constrain between 0 and 1
  })

  sde_mic <- sqrt(diag(t(grd_mic) %*% solve(-est$hessian) %*% grd_mic))

  # Get CI for logit transformed statistic then backtransform to original scale
  ci_logit_mic <- log(object$ftmic/(1-object$ftmic)) + c(-1,1)*qnorm(1-alp/2)*sde_mic
  ci_mic <- exp(ci_logit_mic)/(1 + exp(ci_logit_mic))

  res <- list("ftmic" = object, "95% CI lower" = ci_mic[1], "95% CI upper" = ci_mic[2])

  return(res)
}
