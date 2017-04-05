#' Posterior concentration estimates
#'
#' Fits a two-compartment model to obtain posterior estimates of concentration of drug
#' over time.
#'
#' @param formula A formula where the left side is the measured concentration of drug
#' and the right side is the times of concentration measurements
#' @param data Data frame with concentration data (time of measurement in hours and concentration in mcg/ml)
#' @param pars Vector of pharmacokinetic parameters of length 5: (v_1, k_10, k_12, k_21, err)
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param alp Value of alpha to use for generating pointwise (1 - \code{alp})% confidence bands
#' @param cod Length of time after end of last dose to consider
#' @param thres Threshold for effective treatment (mcg/ml)
#'
#' @details Measurements must be entered in particular units: mcg/ml for concentrations, g/h in rate of
#' infusion, hours for times.
#'
#' @return posterior estimates
#' @export
#'
#' @examples
#' # Insert example from Bayes.R
#'

pkm <- function(formula, data,
                pars = c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7, lerr = 2.33),
                ivt, alp=0.05, cod=12, thres=64) {

  dat <- model.frame(formula, data)
  # Currently, should only have time and concentration data
  if(ncol(dat) > 2){stop("Too many variables in formula")}

  colnames(dat) <- c("time_h", "conc_mcg_ml")

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
  ftmic <- mic_stat(pk_pars = est$par, ivt, dat,
                    tms, con[1,], th = thres)


  obj <- list()
  obj$prior <- log_prior(pars)
  obj$infsched <- ivt
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

