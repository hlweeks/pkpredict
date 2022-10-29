#' Posterior concentration estimates
#'
#' Fits a two-compartment model to obtain posterior estimates of concentration of drug
#' over time.
#'
#' @import stats
#'
#' @param formula A formula where the left side is the measured concentration of drug
#' and the right side is the times of concentration measurements
#' @param data Data frame with concentration data (time of measurement in hours and concentration in mcg/ml)
#' @param subset Subset of the \code{dat} data frame to use
#' @param pars Vector of (prior) log-pharmacokinetic parameters of length 5: (lv_1, lk_10, lk_12, lk_21, lerr)
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param alp Value of alpha to use for generating pointwise (1 - \code{alp})\% confidence bands
#' @param cod Length of time after end of last dose to consider
#' @param thres Threshold for effective treatment (mcg/ml)
#' @param mcmc logical: should estimate of time above threshold be computed using MCMC (false = laplace approximation)
#' @param timeint time interval over which to compute estimate
#' @param nreps number of MCMC iterations to perform (including burn in)
#' @param nburnin number of burn in replications to perform
#' @param seed seed for replicating MCMC results
#' @param shiny is this being used within shiny_pkm
#' @param ... additional arguments (e.g., `mu`, `sig`, `ler_mean`, `ler_sdev` for changing the PK parameter prior mean,
#' variance-covariance matrix and error prior mean and standard deviation, respectively)
#'
#' @details Measurements must be entered in particular units: mcg/ml for concentrations, g/h in rate of
#' infusion, hours for times.
#'
#' @return posterior estimates
#' @export
#'
#' @examples
#' ivt_d <- list(list(begin=0.0, end=0.5, k_R=6),
#'               list(begin=8.0, end=8.5, k_R=6),
#'               list(begin=16.0, end=16.5, k_R=6))
#' dat_d <- data.frame(time_h = c(1,4,40), conc_mcg_ml = c(82.7,80.4,60))
#'
#' pkm(conc_mcg_ml ~ time_h, data = dat_d, ivt = ivt_d)

pkm <- function(formula, data, subset, ivt,
                pars = c(getOption("pkpredict.pip.default.prior")$log_pk_mean,
                         getOption("pkpredict.pip.default.prior")$log_err_mean),
                alp=0.05, cod=12, thres=64,
                timeint = c(0, max(sapply(ivt, function(x) x$end)) + cod),
                mcmc = FALSE, nreps = 5000, nburnin = 2000, seed = NULL, shiny = FALSE, ...) {

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
               control = list(fnscale=-1), hessian=TRUE, ...)

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
                    pars = pars, cod = cod, mcmc = mcmc, shiny = shiny, ...)

  obj <- list("call" = match.call(),
              #"units" = unit,
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

