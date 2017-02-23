#' Time Above Pharmacokinetic Threshold
#'
#' @param pk_pars Vector of pharmacokinetic parameters of length 4: (v_1, k_10, k_12, k_21)
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param tms Times to use in function evaluation
#' @param con Concentration of drug at each of the specified times
#' @param th Threshold value for effective treatment
#'
#' @return Fraction of time spent above the specified threshold from time of first dose through
#' 12 hours after end of the last dose
#'
#' @export

mic_stat <- function(pk_pars, ivt, tms, con, th){
  #Values of the plotted posterior concentrations
  conc <- con[1,]

  # Get PK solution equation evaluated at parameters
  soln <- pk_solution(v_1=exp(pk_pars[1]), k_10=exp(pk_pars[2]),
                      k_12=exp(pk_pars[3]), k_21=exp(pk_pars[4]), ivt=ivt)
  # Use PK solution to define function that computes
  #   the concentrations centered by threshold
  f_mic <- function(times = tms){
    val <- apply(soln(times)*1000, 2, function(x) pmax(0,x)) - th
    return(val[1,])
  }

  # Convert start/end dose times to numeric vectors
  ibe <- sapply(ivt, `[[`, 'begin')
  ied <- sapply(ivt, `[[`, 'end')

  # Initialize time spent above 4*mic
  t_above <- 0

  # Time between start of dose and end of dose
  for(i in 1:length(ibe)){
    # Concentration values at interval endpoints, centered by threshold
    cb <- ifelse(i > 1, conc[tms == ibe[i]] - th, 0-64)
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

  # Use time spent above 4*mic to compute proportion
  frac_time <- t_above/max(tms)

  return(frac_time)
}
