#### PK functions

# Slightly modified function within mic_stat - extracted to compute truth and to compute ft>k x mic with different transformations
# Had to add back in a couple arguments intended for mic_stat, and move a couple key steps within this function
get_stat <- function(pkpars, ivt, th){



  # Full time interval
  tms <- sapply(ivt_i, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+12)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  tms <- unique(pmax(1e-3, tms))
  # Time interval over which to compute ft>k x mic
  ftmic_int <- fs_time_std %>% filter(record_id == opp_id_i) %>%
    ungroup() %>%
    select(ftmic_begin, pk_ftmic_end)
  tms <- tms[tms >= ftmic_int$ftmic_begin & tms <= ftmic_int$pk_ftmic_end]
  # Append start/end times if needed - NEW
  if(tms[1] != ftmic_int$ftmic_begin){tms <- c(ftmic_int$ftmic_begin, tms)}
  if(tms[length(tms)] != ftmic_int$pk_ftmic_end){tms <- c(tms, ftmic_int$pk_ftmic_end)}



  # Get PK solution equation evaluated at parameters
  soln <- pk_solution(v_1=exp(pkpars[1]), k_10=exp(pkpars[2]),
                      k_12=exp(pkpars[3]), k_21=exp(pkpars[4]), ivt=ivt)

  #Concentration values
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
  # Restrict to only times within the tms start/stop - NEW
  ibe <- ibe[ibe >= min(tms) & ibe <= max(tms)]
  ied <- ied[ied >= min(tms) & ied <= max(tms)]

  # Initialize time spent above threshold
  t_above <- 0

  # Concentration values at the interval start and stop points
  cstart <- conc[1]
  cstop <- conc[length(conc)]


  # Time between start of dose and end of dose
  for(i in 1:length(ibe)){
    # Concentration values at interval endpoints, centered by threshold
    cb <- ifelse(i > 1, conc[which.min(abs(tms - ibe[i]))] - th, cstart - th) # NEW LOGIC
    ce <- unique(conc[which.min(abs(tms - ied[i]))] - th) #Printing two copies for some reason? Inserted unique to fix for now
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
    ce <- unique(conc[which.min(abs(tms - ied[j]))] - th)
    c_next <- ifelse(j < length(ied), conc[which.min(abs(tms - ibe[j+1]))] - th,
                     cstop - th)
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
  return(t_above/(max(tms)-min(tms)))
}

# MH functions - currently used within pkpredict package but not exported
metropolis <- function(theta, ivt, dat, Sigma){
  thetastar <- rmvnorm(1, theta, Sigma)

  post_tstar <- log_posterior(lpr = thetastar, ivt, dat = dat)
  post_theta <- log_posterior(lpr = theta, ivt, dat = dat)

  lratio = post_tstar - post_theta

  if(log(runif(1)) < lratio){
    theta <- thetastar
    attr(theta, "soln") <- attributes(post_tstar)$soln
  }else{
    attr(theta, "soln") <- attributes(post_theta)$soln
  }

  return(theta)
}

metro_iterate <- function(nreps = 1000,
                          theta0,
                          return.AR = T,
                          n.reeva.sigma = nreps,
                          ivt, dat, Sigma,
                          shiny = FALSE)
{
  theta <- matrix(NA, nreps, length(theta0))
  soln_list <- vector(mode = "list", length = nreps)
  colnames(theta) <- c("lv_1", "lk_10", "lk_12", "lk_21", "lerr")
  theta[1,] <- theta0
  accept.count = 1


  for(n in 2:nreps){
    metro <- metropolis(theta[n-1,], ivt, dat, Sigma)
    theta[n,] <- metro
    soln_list[[n]] <- attributes(metro)$soln

    if(sum(theta[n,] != theta[n-1,])>0) accept.count = accept.count + 1

    if(n == n.reeva.sigma){
      theta.bar = apply(theta[1:n,], 2, mean)
      theta.dif = sweep(theta[1:n,], 2, theta.bar)
      Sigma = (2.4^2/3)*Reduce("+", lapply(1:n, function(i){
        theta.dif[i,] %*% t(theta.dif[i,])})) / n
    }
  }



  accept.rate = accept.count / nreps

  theta_df <- data.frame(theta = theta)
  attr(theta_df, "soln_list") <- soln_list
  if(return.AR){
    return(list(theta = theta_df, acceptance.rate = accept.rate))
  }else{return(theta_df)}

}

