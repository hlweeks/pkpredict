#'
#' Functions needed for MH and MCMC, currently NOT exported.
#'
#' @param theta pk parameters
#' @param ivt infusion schedule
#' @param dat data
#' @param Sigma covariance matrix




metropolis <- function(theta, ivt, dat, Sigma){
  thetastar <- rmvnorm(1, theta, Sigma)

  post_tstar <- log_posterior(lpr = thetastar, ivt = ivt_d, dat = dat)
  post_theta <- log_posterior(lpr = theta, ivt = ivt_d, dat = dat)

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

  if(shiny){
    shiny::withProgress(message = 'Sampling from the posterior distribution', value = 0, min = 0, max = 1,
                 {for(n in 2:nreps){
                     if(n %% 250 == 0){
                       shiny::incProgress(250/nreps)
                     }

                     theta[n,] <- metropolis(theta[n-1,], ivt, dat, Sigma)
                     if(sum(theta[n,] != theta[n-1,])>0) accept.count = accept.count + 1

                     if(n == n.reeva.sigma){
                       theta.bar = apply(theta[1:n,], 2, mean)
                       theta.dif = sweep(theta[1:n,], 2, theta.bar)
                       Sigma = (2.4^2/3)*Reduce("+", lapply(1:n, function(i){
                         theta.dif[i,] %*% t(theta.dif[i,])})) / n
                     }
                   }
                 })
  }else{
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
  }




  accept.rate = accept.count / nreps
  if(return.AR) return(list(theta = theta, acceptance.rate = accept.rate))

  theta_mat <- data.frame(theta = theta)
  attr(theta_mat, "soln_list") <- soln_list
  if(!return.AR) return(theta_mat)
}

