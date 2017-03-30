#' Posterior concentration estimates
#'
#' @param est Estimate obtained from \code{optim} function with \code{Hessian = TRUE}
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param dat Concentration data frame of the form: data.frame(time_h, conc_mg_dl)
#' @param alp Value of alpha to use for generating pointwise 1 - alpha confidence bands
#' @param cod Length of time after end of last dose to consider
#' @param thres Threshold for effective treatment
#'
#' @return posterior estimates
#' @export
#'
#' @examples
#'

pkm <- function(pars = c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7, lerr = 2.33),
                ivt, dat, alp=0.05, cod=12, thres=64) {

  est <- optim(pars, log_posterior, ivt = ivt, dat = dat,
               control = list(fnscale=-1), hessian=TRUE)


  ## Compute plotting times
  ## - ensure peak and trough times
  ## - avoid time zero
  tms <- sapply(ivt, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  tms <- pmax(1e-3, tms)

  ## Approximate standard deviation of log concentration-time curve
  grd <- fdGrad(est$par, function(pars) {
    sol <- pk_solution(v_1=exp(pars[1]), k_10=exp(pars[2]),
                       k_12=exp(pars[3]), k_21=exp(pars[4]), ivt=ivt)
    log(sol(tms)[1,]*1000) ## mulitply by 1000: g/l -> ug/ml
  })
  sde <- sqrt(diag(grd %*% solve(-est$hessian) %*% t(grd)))
  sde <- ifelse(is.nan(sde), 0, sde)

  ## Posterior estiamte
  sol <- pk_solution(v_1=exp(est$par[1]), k_10=exp(est$par[2]),
                     k_12=exp(est$par[3]), k_21=exp(est$par[4]), ivt=ivt)
  con <- apply(sol(tms)*1000, 2, function(x) pmax(0,x))

  # MIC statistic information
  ftmic <- mic_stat(pk_pars = est$par, ivt, tms, con[1,], th = thres)


  obj <- list()
  obj$prior <- ?
  obj$doseptrn <- ivt
  obj$data <- dat
  obj$optim <- est
  obj$cod <- cod

  obj$times <- tms
  obj$fitted.values <- con
  obj$se_con <- sde

  obj$thresh <- thres
  obj$ftmic <- conf.int(ftmic)



  class(obj) <- "pkm"

  return(obj)
}

