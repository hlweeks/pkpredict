#'
#' Functions needed for MH and MCMC, currently NOT exported.
#'
#'
#'
#'
#' @examples
#'
#'


metropolis <- function(theta, ivt, dat, Sigma){
  thetastar <- rmvnorm(1, theta, Sigma)
  lratio = log_posterior(lpr = thetastar, ivt = ivt_d, dat = dat) - log_posterior(lpr = theta, ivt = ivt_d, dat = dat)
  if(log(runif(1)) < lratio){theta <- thetastar}
  return(theta)
}

metro_iterate <- function(nreps = 1000,
                          theta0,
                          return.AR = T,
                          n.reeva.sigma = nreps,
                          ivt, dat, Sigma)
{
  theta <- matrix(NA, nreps, length(theta0))
  colnames(theta) <- c("lv_1", "lk_10", "lk_12", "lk_21", "lerr")
  theta[1,] <- theta0
  accept.count = 1

  for(n in 2:nreps){
    if(n %% 1000 == 0) print(n)
    theta[n,] <- metropolis(theta[n-1,], ivt, dat, Sigma)
    if(sum(theta[n,] != theta[n-1,])>0) accept.count = accept.count + 1

    if(n == n.reeva.sigma){
      theta.bar = apply(theta[1:n,], 2, mean)
      theta.dif = sweep(theta[1:n,], 2, theta.bar)
      Sigma = (2.4^2/3)*Reduce("+", lapply(1:n, function(i){
        theta.dif[i,] %*% t(theta.dif[i,])})) / n
    }
  }

  accept.rate = accept.count / nreps
  if(return.AR) return(list(theta = theta, acceptance.rate = accept.rate))
  if(!return.AR) return(data.frame(theta = theta))
}

