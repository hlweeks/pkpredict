#' Time Above Pharmacokinetic Threshold
#'
#' Compute an estimate of pharmacodynamic target attainment, typically defined
#' to be some multiple of the minimum inhibitory concentration of the infecting microorganism.
#'
#'
#' @param pars Vector of pharmacokinetic parameters of length 5: (log(v_1), log(k_10), log(k_12), log(k_21), log(err))
#' @param ivt List with containing start of infusion times (h), end of infusion times (h),
#' and rate of infusion (g/h) at each dose
#' @param th Threshold value for effective treatment - check units
#' @param timeint Vector with two elements indicating the start and end of the time interval over which
#' the statistic is computed. By default covers the full dosing period + \code{cod}
#' @param cod Length of time after end of last dose to consider
#' @param conf.level Desired confidence level of the interval
#' @param mcmc logical: should estimate of time above threshold be computed using MCMC (false = laplace approximation)
#' @param dat data
#' @param nreps mcmc replications
#' @param nburnin mcmc burn in iterations
#' @param seed seed for replication
#' @param shiny is shiny being used
#' @param ... additional arguments (e.g., `mu`, `sig`, `ler_mean`, `ler_sdev` for changing the PK parameter prior mean,
#' variance-covariance matrix and error prior mean and standard deviation, respectively)
#'
#' @return Fraction of time spent above the specified threshold from time of first dose through
#' \code{cod} hours after end of the last dose
#'
#' @export
#'
#'

mic_stat <- function(ivt, th, dat = data.frame(),
                     pars = c(getOption("pkpredict.pip.default.prior")$log_pk_mean,
                              getOption("pkpredict.pip.default.prior")$log_err_mean),
                     cod = 12, timeint = c(0, max(sapply(ivt, function(x) x$end)) + cod),
                     conf.level = .95, mcmc = FALSE, nreps = 5000, nburnin = 2000,
                     seed = NULL, shiny = FALSE, ...){

  # Times required for computations
  tms_d <- sapply(ivt, function(x) c(x$begin, x$end))
  tms_d <- c(tms_d, max(tms_d)+cod)
  tms_d <- unique(unlist(sapply(1:(length(tms_d)-1), function(i) {
    s1 <- seq(tms_d[i], tms_d[i+1], 0.1)
    if(s1[length(s1)] != tms_d[i+1]){s1 <- c(s1, tms_d[i+1])}
    return(s1)
  })))
  # Restrict to time interval of interest
  tms_d <- subset(tms_d, tms_d >= timeint[1] & tms_d <= timeint[2])
  tms_d <- pmax(1e-3, tms_d)

  # Only changes pars if nrow(dat) > 0 (default is prior: dat is empty)
  # ... = additional arguments to pass to log_posterior
  est <- optim(pars, log_posterior,
               ivt = ivt, dat = dat,
               control = list(fnscale=-1,maxit=1000), hessian=TRUE, ...)

  # Re-compute if hessian does not produce a valid var-covar matrix
  if(any(diag(solve(-est$hessian)) < 0)){
    est0 <- optim(pars, log_posterior,
                  ivt = ivt, dat = dat,
                  method = "CG",
                  control = list(fnscale=-1,maxit=1000), hessian=TRUE, ...)

    if(any(diag(solve(-est0$hessian)) < 0)){
      warning("Estimated variance-covariance matrix is not positive semi-definite")
    }else{
      # If new hessian results in valid vcov matrix, keep it
      est <- est0
    }
  }


  get_stat <- function(pkpars, ivt_d = ivt, tms=tms_d, th_d = th, inherit.soln = FALSE){

    # Get PK solution equation evaluated at parameters
    if(inherit.soln){
      soln <- attributes(pkpars)$soln

    }else{
      soln <- pk_solution(v_1=exp(pkpars[1]), k_10=exp(pkpars[2]),
                          k_12=exp(pkpars[3]), k_21=exp(pkpars[4]), ivt=ivt_d)
    }

    #Values of the posterior concentrations
    con <- apply(soln(tms)*1000, 2, function(x) pmax(0,x))
    conc <- con[1,]

    # Use PK solution to define function that computes
    #   the concentrations centered by threshold
    f_mic <- function(ts = tms){
      val <- apply(soln(ts)*1000, 2, function(x) pmax(0,x)) - th_d
      return(val[1,])
    }

    # Convert start/end dose times to numeric vectors
    ibe <- sapply(ivt_d, `[[`, 'begin')
    ied <- sapply(ivt_d, `[[`, 'end')

    # Initialize time spent above threshold
    t_above <- 0

    # Time between start of dose and end of dose
    for(i in 1:length(ibe)){
      # Concentration values at interval endpoints, centered by threshold
      cb <- conc[tms == ifelse(ibe[i] == 0, .001, ibe[i])] - th_d #ifelse(i > 1, conc[tms == ibe[i]] - th_d, 0 - th_d)
      ce <- unique(conc[tms == ied[i]] - th_d) #Printing two copies for some reason? Inserted unique to fix for now
      if(cb < 0 && ce > 0){
        # Crosses to above threshold during interval
        root <- uniroot(f_mic, lower = ibe[i], upper = ied[i], tol = .01)$root #Must move from below to above
        t_above <- t_above + (ied[i] - root)
      }else if(cb >= 0 && ce >= 0){
        # Above during whole interval
        t_above <- t_above + (ied[i] - ibe[i])
      }
    }

    # Time between end of one dose and start of the next
    for(j in 1:length(ied)){
      ce <- unique(conc[tms == ied[j]] - th_d)
      c_next <- ifelse(j < length(ied), conc[tms == ibe[j+1]] - th_d, conc[tms == max(tms)] - th_d)
      if(ce > 0 && c_next < 0){
        # Crosses to below threshold during interval
        ulim <- ifelse(j < length(ied), ibe[j+1], max(tms))
        root <- uniroot(f_mic, lower = ied[j], upper = ulim, tol = .01)$root #Must move from above to below
        t_above <- t_above + (root - ied[j])
      }else if(ce >= 0 && c_next >= 0){
        # Above during whole interval
        t_add <- ifelse(j < length(ied), ibe[j+1] - ied[j], max(tms) - ied[j])
        t_above <- t_above + t_add
      }
    }
    return(t_above/max(tms))
  }

  stat <- get_stat(pkpars = est$par)

  # Confidence interval
  alp <- 1 - conf.level

  Sigma0 <- solve(-est$hessian)

  if(mcmc){
    # For reproducibility of sampling
    if(!is.null(seed)){
      set.seed(seed)
    }

    post_samp <- metro_iterate(nreps = nreps, theta0 = est$par,
                               ivt = ivt, dat = dat, Sigma = Sigma0,
                               shiny = shiny)

    theta_samples <- post_samp[[1]][nburnin:nreps,]
    acc_rate <- post_samp[[2]]

    mic_samples <- rep(NA, nrow(theta_samples))
    if(shiny){
      shiny::withProgress(message = 'Computing posterior estimates', value = 0, {
                     for(i in 1:nrow(theta_samples)){
                       shiny::incProgress(1/nrow(theta_samples))
                       mic_samples[i] <- get_stat(theta_samples[i,])
                     }})
    }else{
      mic_samples <- apply(theta_samples, MARGIN = 1, get_stat)#, inherit.soln = TRUE)
    }

    ci_mic <- quantile(mic_samples, probs = c(alp/2, 1 - (alp/2)))
  }else{

    # SE of logit(statistic) using laplace approximation
    grd_mic <- fdGrad(est$par, function(p) {
      mic <- get_stat(pkpars = p)
      # log(mic/(1-mic)) ## use identity transformation (better in new sim results)
    })

    sde_mic <- sqrt(diag(t(grd_mic) %*% Sigma0 %*% grd_mic))

    # # Get CI for logit transformed statistic then backtransform to original scale
    # ci_logit_mic <- log(stat/(1-stat)) + c(-1,1)*qnorm(1-alp/2)*sde_mic
    # ci_mic <- exp(ci_logit_mic)/(1 + exp(ci_logit_mic))
    ci_mic <- stat + c(-1,1)*qnorm(1-alp/2)*sde_mic
  }

  ci_mic[1] <- max(ci_mic[1], 0)
  ci_mic[2] <- min(ci_mic[2], 1)


  # Use time spent above threshold to compute proportion
  ftmic <- list("ftmic" = stat,
                "conf.int" = ci_mic,
                "mcmc" = mcmc,
                "acceptance.rate" = if(mcmc){acc_rate}else{NULL})

  class(ftmic) <- c("mic", class(ftmic))

  return(ftmic)
}


