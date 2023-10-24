##### Estimate what proportion of the posterior samples result in MIC statistics that fall within the bounds of the laplace approximation
# select a single observation

library(pkpredict)
library(parallel)

### Generate sample data
# sample observations from the prior distribution
set.seed(25)
n.samples = 1000
# Prior distribution for PK parameters and error
lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7)
lpk_vcov_d <- 300 * matrix(c(  0.00167, -0.00128,      0,      0,
                               -0.00128,  0.00154,      0,      0,
                               0,        0, .00015,      0,
                               0,        0,      0, .00015), 4, 4)
ler_mean_d <- 2.33
ler_sdev_d <- 0.32

sample.params <- as.data.frame(exp(rmvnorm(n.samples, mean = lpk_mean_d, sigma = lpk_vcov_d)))
sample.params$err <- exp(rnorm(n.samples, ler_mean_d, ler_sdev_d))
names(sample.params) <- c("v_1", "k_10", "k_12", "k_21", "err")

# Set times for blood draws
tm_draws = c(32,32.5,33,34,36,38)

# Calculate predicted concentration for each set of prior parameters
mean.pconc_g_l <- t(apply(sample.params, MARGIN = 1, function(row){
    sol = pk_solution(v_1 = row[1], k_10 = row[2], k_12 = row[3], k_21 = row[4])
    sol(tm_draws)[1,]
  }))


# Convert concentration to mg
mean.pconc_mg_l <- 1000*mean.pconc_g_l

# Draw sample observations from a normal distribution with mean equal to the predicted concentration and the prior variance
sample.obs <- matrix(rnorm(length(tm_draws)*n.samples, mean.pconc_mg_l, sample.params$err),
                     n.samples, length(tm_draws), byrow = F)
#make matrix/vectors of errors and add to concentrations

set.seed(25)
num.clus <- 4
cl <- makeCluster(num.clus)

dat <- data.frame(time_h = tm_draws, conc_mg_dl = NA)


# coverage_logit = coverage_probit = coverage_frac = rep(NA, n.samples)
# ci_mic_probit = ci_mic_frac = ci_mic_logit = matrix(NA, nrow = n.samples, 2)
# ncov_logit = ncov_probit = ncov_frac = data.frame("lower" = rep(NA, n.samples),
#                                                   "upper" = rep(NA, n.samples))


for(i in 1:n.samples){
  print(paste("Calculating coverage for sample",i))

  obs = sample.obs[i,]
  lpr <- log(sample.params[i,])
  dat$conc_mg_dl <- sample.obs[i,]

  # get laplace approximation
  est <- optim(lpr_mean_d, # which parameters to optimize - starting values
               log_posterior, # function to minimize -- this includes the lpr values as an argument
               ivt=ivt_d, # dosage schedule - argument to log_posterior, loaded from "model.R"
               dat=dat,
               control = list(fnscale=-1),
               hessian=TRUE)

  # # estimate theta using metropolis algorithm
  Sigma0 = solve(-est$hessian)
  res = metro_iterate(nreps = nreps, theta0 = lpr_mean_d, ivt = ivt_d, dat = dat)
  theta = res[[1]][nburnin:nreps,]

  # calculate laplace MIC bounds
  # Estimate time spent above threshold
  tms <- sapply(ivt_d, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  tms <- pmax(1e-3, tms)

  # estimate percent of time above threshold / se for Laplace approx
  sol_lap <- pk_solution(v_1=exp(est$par[1]), k_10=exp(est$par[2]),
                         k_12=exp(est$par[3]), k_21=exp(est$par[4]), ivt=ivt_d)
  con_lap <- apply(sol_lap(tms)*1000, 2, function(x) pmax(0,x))

  frac_mic_lap <- mic_stat(pk_pars = est$par, ivt = ivt_d, tms, con_lap, th = thres)

  # calculate gradients for various transformations: logit, probit, identity
  grd_mic_logit <- fdGrad(est$par, function(pars) {
    mic <- mic_stat(pk_pars = pars, ivt_d, tms, con_lap, th = thres)
    log(mic/(1-mic)) ## constrain between 0 and 1
  })

  grd_mic_probit <- fdGrad(est$par, function(pars) {
    mic <- mic_stat(pk_pars = pars, ivt_d, tms, con_lap, th = thres)
    qnorm(mic)
  })

  grd_mic_frac <- fdGrad(est$par, function(pars) {
    mic_stat(pk_pars = pars, ivt_d, tms, con_lap, th = thres)
  })

  # estimate standard errors using delta method
  sde_mic_logit <- sqrt(diag(t(grd_mic_logit) %*% solve(-est$hessian) %*% grd_mic_logit))
  sde_mic_probit <- sqrt(diag(t(grd_mic_probit) %*% solve(-est$hessian) %*% grd_mic_probit))
  sde_mic_frac <- sqrt(diag(t(grd_mic_frac) %*% solve(-est$hessian) %*% grd_mic_frac))

  # Get CI for transformed statistic then backtransform to original scale
  ci_logit_mic_lap <- log(frac_mic_lap/(1-frac_mic_lap)) + c(-1,1)*qnorm(1-alp/2)*sde_mic_logit
  ci_probit_mic_lap <- qnorm(frac_mic_lap) + c(-1,1)*qnorm(1-alp/2)*sde_mic_probit
  ci_frac_mic_lap <- frac_mic_lap + c(-1,1)*qnorm(1-alp/2)*sde_mic_frac

  ci_mic_logit[i,] <- exp(ci_logit_mic_lap)/(1 + exp(ci_logit_mic_lap))
  ci_mic_probit[i,] <- pnorm(ci_probit_mic_lap)
  ci_mic_frac[i,] <- ci_frac_mic_lap

  # draw samples from posterior distribution and calculate MIC statistic for each
  theta.samples = theta

  mic.samples = parApply(cl, theta.samples, 1, function(vals){
    source("pk_functions.R")
    thres = 64
    alp=0.05
    cod=12

    tms <- sapply(ivt_d, function(x) c(x$begin, x$end))
    tms <- c(tms, max(tms)+cod)
    tms <- unlist(sapply(1:(length(tms)-1), function(i) {
      s1 <- seq(tms[i], tms[i+1], 1/10)
      if(tms[i+1] %% 1/10)
        s1 <- c(s1, tms[i+1])
      return(s1)
    }))
    tms <- pmax(1e-3, tms)

    sol_sample <- pk_solution(v_1=exp(vals[1]), k_10=exp(vals[2]),
                              k_12=exp(vals[3]), k_21=exp(vals[4]), ivt=ivt_d)
    con_sample <- apply(sol_sample(tms)*1000, 2, function(x) pmax(0,x))
    frac_mic_sample <- mic_stat(pk_pars = vals, ivt = ivt_d, tms, con_sample, th = thres)
    return(frac_mic_sample)
  })

  coverage_logit[i] = mean(mic.samples >= ci_mic_logit[i,1] & mic.samples <= ci_mic_logit[i,2])
  coverage_probit[i] = mean(mic.samples >= ci_mic_probit[i,1] & mic.samples <= ci_mic_probit[i,2])
  coverage_frac[i] = mean(mic.samples >= ci_mic_frac[i,1] & mic.samples <= ci_mic_frac[i,2])

  ncov_logit$lower[i] = mean(mic.samples < ci_mic_logit[i,1])
  ncov_logit$upper[i] = mean(mic.samples > ci_mic_logit[i,2])
  ncov_probit$lower[i] = mean(mic.samples < ci_mic_probit[i,1])
  ncov_probit$upper[i] = mean(mic.samples > ci_mic_probit[i,2])
  ncov_frac$lower[i] = mean(mic.samples < ci_mic_frac[i,1])
  ncov_frac$upper[i] = mean(mic.samples > ci_mic_frac[i,2])
}

# FINISHED

# coverage_res = data.frame(coverage_logit, coverage_probit, coverage_frac)
# # save(coverage_res, file = "cov.RData")
# # save(ncov_logit, ncov_probit, ncov_frac, file = "all_ncov.rdta")
#
# save(coverage_res, ncov_logit, ncov_probit, ncov_frac,
#      file = "all_cov.rdta")
# load("coverage_results.RData")

stopCluster(cl)

#########################################################
#########################################################
#########################################################
#########################################################
#########################################################


# plot coverage probabilities
plot(1:n.sample.coverage, coverage_res[,1], xlab = "sample", ylab = "posterior coverage", ylim = c(0.5,1))
points(1:n.sample.coverage, coverage_res[,2], col = 2)
points(1:n.sample.coverage, coverage_res[,3], col = 4)

# summary of each method
apply(coverage_res,2,summary,na.rm = T)

par(mfrow = c(3,1))
hist(coverage_res[,1], breaks = 20)
hist(coverage_res[,2], breaks = 20)
hist(coverage_res[,3], breaks = 20)
par(mfrow = c(1,1))

# estimate monte carlo error
mc.error = function(p,N) sqrt(p*(1-p)/(N*p))
mean(mc.error(coverage_res[,1], nreps - nburnin), na.rm = T)
mean(mc.error(coverage_res[,2], nreps - nburnin), na.rm = T)
mean(mc.error(coverage_res[,3], nreps - nburnin), na.rm = T)


# identify and graph outliers
outlier.ids = which(apply(coverage_res < 0.8, 1, sum) >0)
outlier.obs = sample.obs[outlier.ids,]
outlier.params = sample.params[outlier.ids,]

par(mfrow = c(length(outlier.ids),1))
for(obs in 1:length(outlier.ids)){
  outlier.sol = pk_solution(v_1 = outlier.params[obs,1],
                            k_10 = outlier.params[obs,2],
                            k_12 = outlier.params[obs,3],
                            k_21 = outlier.params[obs,4],
                            ivt = ivt_d)
  outlier.con <- apply(outlier.sol(tms)*1000, 2, function(x) pmax(0,x))
  plot(tms, outlier.con[1,], xlab="Time (h)", ylab="Central Concentration (g/L)", type='l')
  points(dat$time_h, outlier.obs[obs,])

}
par(mfrow = c(1,1))
