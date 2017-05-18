#' Time Above Pharmacokinetic Threshold
#'
#' Compute an estimate of pharmacodynamic target attainment, typically defined
#' to be some multiple of the minimum inhibitory concentration of the infecting microorganism
#'
#' @param pars Vector of pharmacokinetic parameters of length 5: (log(v_1), log(k_10), log(k_12), log(k_21), log(err))
#' @param ivt List with containing start of infusion times (h), end of infusion times (h),
#' and rate of infusion (g/h) at each dose
#' @param th Threshold value for effective treatment - check units
#' @param timeint Vector with two elements indicating the start and end of the time interval over which
#' the statistic is computed
#' @param cod Length of time after end of last dose to consider
#'
#' @return Fraction of time spent above the specified threshold from time of first dose through
#' \code{cod} hours after end of the last dose
#'
#' @export
#'
#' @examples
#'
#'

mic_stat <- function(ivt, th,
                     pars = c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7, lerr = 2.33),
                     cod = 12, timeint = c(0, max(sapply(ivt, function(x) x$end)) + cod)){

  # Times required for computations
  tms <- sapply(ivt, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  # Restrict to time interval of interest
  tms <- subset(tms, tms >= timeint[1] & tms <= timeint[2])
  tms <- pmax(1e-3, tms)


  # Get PK solution equation evaluated at parameters
  soln <- pk_solution(v_1=exp(pars[1]), k_10=exp(pars[2]),
                      k_12=exp(pars[3]), k_21=exp(pars[4]), ivt=ivt)
  #Values of the posterior concentrations
  con <- apply(soln(tms)*1000, 2, function(x) pmax(0,x))
  conc <- con[1,]

  # Use PK solution to define function that computes
  #   the concentrations centered by threshold
  f_mic <- function(ts = tms){
    val <- apply(soln(ts)*1000, 2, function(x) pmax(0,x)) - th
    return(val[1,])
  }

  # Convert start/end dose times to numeric vectors
  ibe <- sapply(ivt, `[[`, 'begin')
  ied <- sapply(ivt, `[[`, 'end')

  # Initialize time spent above threshold
  t_above <- 0

  # Time between start of dose and end of dose
  for(i in 1:length(ibe)){
    # Concentration values at interval endpoints, centered by threshold
    cb <- ifelse(i > 1, conc[tms == ibe[i]] - th, 0 - th)
    ce <- unique(conc[tms == ied[i]] - th) #Printing two copies for some reason? Inserted unique to fix for now
    if(cb < 0 && ce > 0){
      # Crosses to above threshold during interval
      root <- uniroot(f_mic, lower = ibe[i], upper = ied[i])$root #Must move from below to above
      t_above <- t_above + (ied[i] - root)
    }else if(cb >= 0 && ce >= 0){
      # Above during whole interval
      t_above <- t_above + (ied[i] - ibe[i])
    }
  }

  # Time between end of one dose and start of the next
  for(j in 1:length(ied)){
    ce <- unique(conc[tms == ied[j]] - th)
    c_next <- ifelse(j < length(ied), conc[tms == ibe[j+1]] - th, conc[tms == max(tms)] - th)
    if(ce > 0 && c_next < 0){
      # Crosses to below threshold during interval
      ulim <- ifelse(j < length(ied), ibe[j+1], max(tms))
      root <- uniroot(f_mic, lower = ied[j], upper = ulim)$root #Must move from above to below
      t_above <- t_above + (root - ied[j])
    }else if(ce >= 0 && c_next >= 0){
      # Above during whole interval
      t_add <- ifelse(j < length(ied), ibe[j+1] - ied[j], max(tms) - ied[j])
      t_above <- t_above + t_add
    }
  }


  # Use time spent above threshold to compute proportion
  ftmic <- list("ftmic" = t_above/max(tms),
                "pars" = pars,
                "ivt" = ivt, "th" = th)

  class(ftmic) <- c("mic", class(ftmic))

  return(ftmic)
}
