#'
#' Functions needed for MH and MCMC, currently not exported.
#'
#' @param theta pk parameters
#' @param ivt infusion schedule
#' @param dat data
#' @param Sigma covariance matrix
#' @param ... additional arguments (e.g., `mu`, `sig`, `ler_mean`, `ler_sdev` for changing the PK parameter prior mean,
#' variance-covariance matrix and error prior mean and standard deviation, respectively)




metropolis <- function(theta, ivt, dat, Sigma, ...){
  # Sample point from proposal distribution
  thetastar <- rmvnorm(1, theta, Sigma)

  post_tstar <- log_posterior(lpr = thetastar, ivt, dat = dat, ...)
  post_theta <- log_posterior(lpr = theta, ivt, dat = dat, ...)

  # log-LR
  lratio = post_tstar - post_theta

  # Accept point with some probability
  if(log(runif(1)) < lratio){
    theta_keep <- thetastar
    # attr(theta, "soln") <- attributes(post_tstar)$soln
  }else{
    theta_keep <- theta
    # attr(theta, "soln") <- attributes(post_theta)$soln
  }

  return(theta_keep)
}

metro_iterate <- function(nreps = 1000,
                          theta0,
                          ivt, dat, Sigma,
                          shiny = FALSE)
{
  # Initialize objects
  theta <- matrix(NA, nreps, length(theta0))
  colnames(theta) <- c("lv_1", "lk_10", "lk_12", "lk_21", "lerr")
  # soln_list <- vector(mode = "list", length = nreps)
  accept = rep(NA, nreps)

  # Starting theta guess
  theta[1,] <- theta0

  if(shiny){
    shiny::withProgress(message = 'Sampling from the posterior distribution', value = 0, min = 0, max = 1,
                 {for(n in 2:nreps){
                     if(n %% 250 == 0){
                       shiny::incProgress(250/nreps)
                     }

                     theta[n,] <- metropolis(theta[n-1,], ivt, dat, Sigma)
                     accept[n] <- ifelse(all.equal(theta[n,], theta[n-1,]), 0, 1)

                   }
                 })
  }else{
    for(i in 2:nreps){
      metro <- metropolis(theta[i-1,], ivt, dat, Sigma)
      theta[i,] <- metro
      # soln_list[[i]] <- attributes(metro)$soln

      accept[i] <- ifelse(all.equal(theta[i,], theta[i-1,]) == TRUE,
                          0, 1)

    }
  }

  # Remove first element (NA)
  accept <- accept[-1]

  return(list(theta_df = data.frame(theta),
              acceptance = accept))

}

