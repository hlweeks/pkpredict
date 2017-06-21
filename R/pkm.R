#' Posterior concentration estimates
#'
#' Fits a two-compartment model to obtain posterior estimates of concentration of drug
#' over time.
#'
#' @import stats graphics
#'
#' @param formula A formula where the left side is the measured concentration of drug
#' and the right side is the times of concentration measurements
#' @param data Data frame with concentration data (time of measurement in hours and concentration in mcg/ml)
#' @param subset Subset of the \code{dat} data frame to use
#' @param pars Vector of (prior) log-pharmacokinetic parameters of length 5: (lv_1, lk_10, lk_12, lk_21, lerr)
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
#' pkm(concentration ~ time, dat_d, ivt_d) # something like this

pkm <- function(formula, data, subset, ivt,
                pars = c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7, lerr = 2.33),
                alp=0.05, cod=12, thres=64,
                timeint = c(0, max(sapply(ivt, function(x) x$end)) + cod), ...) {

  # Allows formula, data, and subset to be optional (for prior only)
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mc), 0L)

  if(sum(m) > 0){
    mc <- mc[c(1L, m)]
    mc$drop.unused.levels <- TRUE
    mc[[1L]] <- as.name("model.frame")
    dat <- eval(mc, parent.frame())


    # Currently, should only have time and concentration data
    if(ncol(dat) > 2){stop("Too many variables in formula")}

    colnames(dat) <- c("conc_mcg_ml", "time_h")
    dat <- dat[, c("time_h", "conc_mcg_ml")]
  }else{
    dat <- data.frame()
  }

  # Only depends on pars if nrow(dat) == 0
  est <- optim(pars, log_posterior, ivt = ivt, dat = dat,
               control = list(fnscale=-1), hessian=TRUE)

  # Times at which to compute concentration estimates and SE: dose and concentration meas times
  tms <- c(sapply(ivt, function(x) x$begin),
           sapply(ivt, function(x) x$end))
  if(nrow(dat) > 0){tms <- c(tms, dat$time_h)}
  tms <- sort(unique(tms))
  # Prevent log(0)
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
  ftmic <- mic_stat(ivt = ivt, th = thres, dat = dat,
                    pars = pars, cod = cod)

  obj <- list(#"call" = match.call(),
              # Posterior estimate
              "optim" = est,
              # Prior information
              # NEED TO GET PRIOR ESTIMATES??

              # Relevant data
              "infsched" = ivt,
              "data" = dat,
              "alpha" = alp,
              "cod" = cod,

              # At infusion and observed times
              "fitted.values" = data.frame("time" = tms,
                                           "conc" = con[1,],
                                           "se_con" = sde),

              # Needed for predict method
              "sol.eqn" = sol,

              # ft>mic information
              "thresh" = thres,
              "ftmic" = ftmic
              )

  class(obj) <- "pkm"

  return(obj)
}

