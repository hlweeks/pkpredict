##### Estimate what proportion of the posterior samples result in MIC statistics that fall within the bounds of the laplace approximation
# select a single observation

library(pkpredict)
library(parallel)

#source("sim_help_funs.R")
par_names <- c("lv_1", "lk_10", "lk_12", "lk_21", "lerr")

set.seed(314)

# Number of "patients" to sample from the prior
n_samples = 1000
# Number of MCMC reps to perform and to consider burnin
nreps = 5000
nburnin = 2000


# Alpha used for laplace-approximated credible intervals
alp <- .05
# Prior distribution for PK parameters and error - from kde_emp_prior.RData
# written explicitly so sims can be done in statcomp2
emp_prior_kde <- list("log_pk_mean" = c("lv_1"=2.7760840, "lk10"=-0.8777951,
                                        "lk12"=-4.7726861, "lk21"=-3.9570281),
                      "log_pk_vcov" = matrix(c(0.203362856, -0.009494482, -0.0025150202,  0.0011462905,
                                               -0.009494482,  0.153998510, -0.0027744820, -0.0064576674,
                                               -0.0025150202, -0.0027744820,  0.0013094598,  0.0003188164,
                                               0.0011462905, -0.0064576674,  0.0003188164,  0.0012682863), nrow=4,ncol=4,byrow=T),
                      "log_err_mean"=-0.2741756, "log_err_sd"=0.02111675)


# rmvnorm from mvtnorm package
log_sample_params <- as.data.frame(rmvnorm(n_samples,
                                       mean = emp_prior_kde$log_pk_mean,
                                       sigma = emp_prior_kde$log_pk_vcov))
log_sample_params$err <- rnorm(n_samples, emp_prior_kde$log_err_mean, emp_prior_kde$log_err_sd)
names(log_sample_params) <- par_names


# patient data
source('explore_opppk_data.R')

# which ivt/cod/sampletime set to use
opp_id <- sample(1:length(ivt_master), size = n_samples, replace=T)


# Concentration measurements at these times
# "Observed" concentration measurements with sigma error
sim_conc <- lapply(1:n_samples, function(i){

  opp_id_i <- names(ivt_master)[opp_id[i]]

  # Solution for this patient
  sol <- pk_solution(k_10 = exp(log_sample_params[i,'lk_10']),
                     k_12 = exp(log_sample_params[i,'lk_12']),
                     k_21 = exp(log_sample_params[i,'lk_21']),
                     v_1 = exp(log_sample_params[i,'lv_1']),
                     ivt=ivt_master[[opp_id_i]])

  # Apply to all draw times with measurement error, convert to mg/L
  true_conc <- sol(samptime_master[samptime_master$record_id==opp_id_i,
                                   'reltime_h'])[1,]*1000
  # Add measurement error (sd on mg/L scale)
  sol_err <- true_conc + rnorm(sum(samptime_master$record_id==opp_id_i),
                               mean = 0, sd = exp(log_sample_params[i,'lerr']))
  sol_err <- pmax(sol_err, 1e-3)

  names(sol_err) <- NULL
  # Convert to mg/L
  return(sol_err)
})


setup <- list("log_sample_params" = log_sample_params,
              "opp_id" = opp_id,
              "sim_conc" = sim_conc)
save(setup, file = "../Data/setup.RData")

# Initialize objects that will hold results

bias_res <- data.frame(true_ftmic = rep(NA, n_samples),
                       model_ftmic = rep(NA, n_samples))

grd_list <- list("ident" = matrix(ncol=5,nrow=n_samples),
                 "logit" = matrix(ncol=5,nrow=n_samples),
                 "probit" = matrix(ncol=5,nrow=n_samples))
colnames(grd_list$ident) <- colnames(grd_list$logit) <- colnames(grd_list$probit) <- par_names

est_list <- vector("list", n_samples)

post_pkpar_list <- vector("list", n_samples)

post_ftmic_list <- vector("list", n_samples)

bcov_res <- data.frame(ident = rep(NA, n_samples),
                       logit = rep(NA, n_samples),
                       probit = rep(NA, n_samples))

nbcov_res <- data.frame(ident_lower = rep(NA, n_samples),
                        ident_upper = rep(NA, n_samples),
                        logit_lower = rep(NA, n_samples),
                        logit_upper = rep(NA, n_samples),
                        probit_lower = rep(NA, n_samples),
                        probit_upper = rep(NA, n_samples))




# Function to write status file for monitoring progress
t0 <- Sys.time()
write_status <- function(j, n = n_samples){
  status <- paste0("On sample ", j, " of ", n, ".\n",
                   (j-1)/n*100,"% complete.\n",
                   "Elapsed time: ", round(difftime(Sys.time(), t0, units = "mins"), 3),
                   " minutes")

  write(status, file="status.txt")
}


#### Temporary: using partial to ensure prior parameters are correct ####
log_post <- purrr::partial(log_posterior,
                           mu = emp_prior_kde$log_pk_mean, sig = emp_prior_kde$log_pk_vcov,
                           ler_mean = emp_prior_kde$log_err_mean, ler_sdev = emp_prior_kde$log_err_sd)
# Save repetitive arguments
get_mic_stat <- purrr::partial(get_stat, ivt = ivt_d, th = 64)

# time after last dose to use in ftmic computation
cod <- 8

# Compute ft>k x mic (truth and estimates) as well as 'Bayesian coverage' for each sample drawn from the prior
for(i in 1:n_samples){
  write_status(j=i)

  # log-pk parameter vector and data frame for sample i
  lpr <- log_sample_params[i,]

  # Get opp_id used to define ivt and dat
  opp_id_i <- names(ivt_master)[opp_id[i]]

  ivt_d <- ivt_master[[opp_id_i]]
  dat <- data.frame(time_h = samptime_master[samptime_master$record_id==opp_id_i,
                                             'reltime_h_'],
                    conc_mgL = NA)


  # Time interval over which to compute ft>k x mic
  tms <- sapply(ivt_d, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    print(i)
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  tms <- pmax(1e-3, tms)

  # Compute truth based on sampled parameters
  true_ftmic <- get_mic_stat(pkpars = lpr)

  # Get MAP estimate of pk parameters
  est <- optim(c(emp_prior_kde$log_pk_mean, emp_prior_kde$log_err_mean), # which parameters to optimize - starting values
               log_post, # function to minimize -- this includes the lpr values as an argument
               ivt=ivt_d, # dosage schedule - argument to log_posterior
               dat=dat, # observed concentration measurements - argument to log_posterior
               control = list(fnscale=-1),
               hessian=TRUE)

  # Get model-based estimate of ft>k x mic
  model_ftmic <- get_mic_stat(pkpars = est$par)


  # Calculate gradients for various transformations: logit, probit, identity
  grd_ident <- fdGrad(est$par, function(pars) {
    get_mic_stat(pkpars = pars)
  })

  grd_logit <- fdGrad(est$par, function(pars) {
    mic <- get_mic_stat(pkpars = pars)
    log(mic/(1-mic))
  })

  grd_probit <- fdGrad(est$par, function(pars) {
    mic <- get_mic_stat(pkpars = pars)
    qnorm(mic)
  })

  # Estimate standard errors on transformed scales using delta method
  sde_ident <- sqrt(diag(t(grd_ident) %*% solve(-est$hessian) %*% grd_ident))
  sde_logit <- sqrt(diag(t(grd_logit) %*% solve(-est$hessian) %*% grd_logit))
  sde_probit <- sqrt(diag(t(grd_probit) %*% solve(-est$hessian) %*% grd_probit))

  # Get laplace-approximated CI for transformed statistic then backtransform to original scale
  ci_ident_trans <- model_ftmic + c(-1,1)*qnorm(1-alp/2)*sde_ident
  ci_logit_trans <- log(model_ftmic/(1-model_ftmic)) + c(-1,1)*qnorm(1-alp/2)*sde_logit
  ci_probit_trans <- qnorm(model_ftmic) + c(-1,1)*qnorm(1-alp/2)*sde_probit

  # Transform approximate CI back to original scale
  ci_ident <- ci_ident_trans
  ci_logit <- exp(ci_logit_trans)/(1 + exp(ci_logit_trans))
  ci_probit <- pnorm(ci_probit_trans)


  # Now, need to determine the true proportion of the posterior contained within each of these approcimate intervals
  # estimate theta using MH algorithm
  Sigma0 <- solve(-est$hessian)
  # nreps x 5 (dim pkpars plus err par) matrix
  post_samples <- metro_iterate(nreps = nreps, theta0 = c(emp_prior_kde$log_pk_mean, emp_prior_kde$log_err_mean),
                                ivt = ivt_d, dat = dat, Sigma = Sigma0)
  theta <- post_samples$theta[(nburnin+1):nreps,]


  # Calculate ft>k x mic for each of the posterior samples
  post_ftmic <- as.numeric(apply(theta, MARGIN = 1, get_mic_stat))


  # Determine proportion of posterior samples within each model-based interval type
  bcov_res$ident[i] = mean(post_ftmic >= ci_ident[1] & post_ftmic <= ci_ident[2])
  bcov_res$logit[i] = mean(post_ftmic >= ci_logit[1] & post_ftmic <= ci_logit[2])
  bcov_res$probit[i] = mean(post_ftmic >= ci_probit[1] & post_ftmic <= ci_probit[2])

  # Split non-coverage into upper vs lower to assess symmetry
  nbcov_res$ident_lower[i] = mean(post_ftmic < ci_ident[1])
  nbcov_res$ident_upper[i] = mean(post_ftmic > ci_ident[2])
  nbcov_res$logit_lower[i] = mean(post_ftmic < ci_logit[1])
  nbcov_res$logit_upper[i] = mean(post_ftmic > ci_logit[2])
  nbcov_res$probit_lower[i] = mean(post_ftmic < ci_probit[1])
  nbcov_res$probit_upper[i] = mean(post_ftmic > ci_probit[2])

  bcov_list <- list("bcov_res" = bcov_res, "nbcov_res" = nbcov_res)

  # Add other results for sample i to overall results
  bias_res$true_ftmic[i] = true_ftmic
  bias_res$model_ftmic[i] = model_ftmic

  grd_list$ident[i,] = grd_ident
  grd_list$logit[i,] = grd_logit
  grd_list$probit[i,] = grd_probit

  est_list[[i]] <- est

  post_ftmic_list[[i]] <- post_ftmic

  post_list <- list("est_list" = est_list,
                    "post_ftmic_list" = post_ftmic_list)

  # Write results to file with each loop iteration
  save(bias_res, file = "results/bias_res.RData")
  save(grd_list, file = "results/grd_list.RData")
  save(post_list, file = "results/sample_post.RData")
  save(bcov_list, file = "results/bcoverage_res.RData")
}

write_status(j=n_samples)
