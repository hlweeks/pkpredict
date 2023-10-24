#' Time Above Pharmacokinetic Threshold
#'
#' Compute an estimate of pharmacodynamic target attainment, typically defined
#' to be some multiple of the minimum inhibitory concentration of the infecting microorganism.
#'
#'
#' @param pars Vector of pharmacokinetic parameters of length 5: (log(v_1), log(k_10), log(k_12), log(k_21), log(err))
#' @param ivt List with containing start of infusion times (h), end of infusion times (h),
#' and rate of infusion (g/h) at each dose. First infusion within \code{ivt} must be a true initial infusion,
#' i.e. concentration of the drug in the patient's system must be 0 prior to that time.
#' @param th Threshold value for effective treatment - check units
#' @param timeint Vector with two elements indicating the start and end of the time interval over which
#' the statistic is computed. By default covers the full dosing period + \code{cod}
#' @param cod Length of time after end of last dose to consider
#' @param conf.level Desired confidence level of the interval
#' @param mcmc logical: should estimate of time above threshold be computed using MCMC (false = laplace approximation)
#' @param dat data
#' @param nreps mcmc replications
#' @param nburnin mcmc burn in iterations
#' @param nthin mcmc thinning interval
#' @param seed seed for replication
#' @param shiny is shiny being used
#' @param ... additional arguments (e.g., `mu`, `sig`, `ler_mean`, `ler_sdev` for changing the PK parameter prior mean,
#' variance-covariance matrix and error prior mean and standard deviation, respectively)
#'
#' @return Fraction of time spent above the specified threshold within a certain time interval (by default: time of first dose through
#' \code{cod} hours after end of the last dose).
#'
#' @export
#'
#'

mic_stat <- function(ivt, th, dat = data.frame(),
                     pars = c(getOption("pkpredict.pip.default.prior")$log_pk_mean,
                              getOption("pkpredict.pip.default.prior")$log_err_mean),
                     cod = 12, timeint = NULL,
                     conf.level = .95, mcmc = FALSE, nreps = 5000, nburnin = 2000, nthin = 10,
                     seed = NULL, shiny = FALSE, ...){

  # Times required for computations
  tms0 <- sapply(ivt, function(x) c(x$begin, x$end))
  tms0 <- sort(c(tms0, max(tms0)+cod))
  if(!is.null(timeint)){tms0 <- c(tms0, timeint)}
  tms0 <- sort(tms0)
  tms0 <- unique(unlist(sapply(1:(length(tms0)-1), function(i) {
    s1 <- seq(tms0[i], tms0[i+1], 0.1)
    if(s1[length(s1)] != tms0[i+1]){s1 <- c(s1, tms0[i+1])}
    return(s1)
  })))
  # Restrict to time interval of interest
  if(!is.null(timeint)){
    if(timeint[1] < tms0[1] | timeint[2] > tms0[length(tms0)]){
      stop("'timeint' out of bounds")
    }

    if(timeint[2] <= timeint[1]){
      stop("'timeint' has interval end <= interval start")
    }

    tms <- tms0[tms0 >= timeint[1] & tms0 <= timeint[2]]
  }else{
    # No time interval provided - do full dosing history
    tms <- tms0
  }
  tms <- pmax(1e-5, tms)

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


  get_stat <- function(pkpars, ivt_d = ivt, tms_d=tms, th_d = th){

    # Get PK solution equation evaluated at parameters
    soln <- pk_solution(v_1=exp(pkpars[1]), k_10=exp(pkpars[2]),
                        k_12=exp(pkpars[3]), k_21=exp(pkpars[4]),
                        ivt=ivt_d)

    #Values of the posterior concentrations
    con <- apply(soln(tms_d)*1000, 2, function(x) pmax(0,x))
    conc <- con[1,]

    # Use PK solution to define function that computes
    #   the concentrations centered by threshold
    f_mic <- function(ts = tms_d){
      val <- apply(soln(ts)*1000, 2, function(x) pmax(0,x)) - th_d
      return(val[1,])
    }

    # Helper function to identify root within sub-intervals
    time_add_root <- function(conc_start, conc_end,
                              t_start, t_end){

      # If no conditions are met, add 0
      t_add <- 0

      # T/F is concentration increasing during interval
      conc_incr <- conc_end > conc_start

      # Relabel higher/lower concentration
      conc_lower <- min(conc_start, conc_end)
      conc_upper <- max(conc_start, conc_end)


      # Concentration crosses threshold during interval
      if(conc_lower < 0 & conc_upper >= 0){

        root <- uniroot(f_mic,
                        lower = t_start,
                        upper = t_end,
                        tol = .01)$root

        if(conc_incr){
          # If conc increasing, add later part of interval
          t_add <- t_end - root
        }else{
          # If conc decreasing, add earlier part of interval
          t_add <- root - t_start
        }

      }else if(conc_lower >= 0 & conc_upper >= 0){
        # Above during whole interval
        t_add <- t_end - t_start
      }

      # Return time to be added to t_above
      return(t_add)

    }

    # Convert start/end dose times to numeric vectors
    ibe <- pmax(1e-5, sapply(ivt_d, `[[`, 'begin'))
    ied <- sapply(ivt_d, `[[`, 'end')

    # Determine endpoints needed for root checking
    root_times <- unique(c(ibe, ied, max(tms_d), timeint))
    if(!is.null(timeint)){
      root_times <- sort(subset(root_times,
                                root_times >= timeint[1],
                                root_times <= timeint[2]))
    }else{
      root_times <- sort(root_times)
    }

    t_above_interval <- sapply(1:(length(root_times) - 1), function(ix){

      time_add_root(conc_start = conc[tms_d == root_times[ix]] - th_d,
                    conc_end = conc[tms_d == root_times[ix+1]] - th_d,
                    t_start = root_times[ix],
                    t_end = root_times[ix+1])

    }, USE.NAMES = FALSE)

    # ftmic = proportion of timeint with concentration above threshold
    t_above <- sum(t_above_interval)/(max(tms_d) - min(tms_d))

    return(t_above)
  }

  # ftmic stat based on posterior pkpar estimate
  stat <- get_stat(pkpars = est$par)

  # Confidence interval
  alp <- 1 - conf.level
  # Posterior-estimated var-covar matrix
  Sigma0 <- solve(-est$hessian)

  if(mcmc){
    if(nburnin >= nreps){
      stop("nburnin must be less than nreps")
    }

    # For reproducibility of sampling
    if(!is.null(seed)){
      set.seed(seed)
    }

    post_samp <- metro_iterate(nreps = nreps, theta0 = est$par,
                               ivt = ivt, dat = dat, Sigma = Sigma0,
                               shiny = shiny)

    theta_samples <- post_samp$theta_df[((nburnin+1):nreps)%%nthin == 0,]
    acc_rate <- mean(post_samp$acceptance)

    mic_samples <- rep(NA, nrow(theta_samples))
    if(shiny){
      shiny::withProgress(message = 'Computing posterior estimates', value = 0, {
                     for(i in 1:nrow(theta_samples)){
                       shiny::incProgress(1/nrow(theta_samples))
                       mic_samples[i] <- get_stat(theta_samples[i,])
                     }})
    }else{
      mic_samples <- apply(theta_samples, MARGIN = 1, get_stat)
    }

    ci_mic <- quantile(mic_samples, probs = c(alp/2, 1 - (alp/2)))
  }else{
    # Non-MCMC

    # SE of logit(statistic) using laplace approximation
    grd_mic <- fdGrad(est$par, function(p) {
      mic <- get_stat(pkpars = p)
      # log(mic/(1-mic)) ## use "logit" transformation
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
                "est" = est,
                "mcmc" = mcmc,
                "acceptance.rate" = if(mcmc){acc_rate}else{NULL})

  class(ftmic) <- c("mic", class(ftmic))

  return(ftmic)
}


