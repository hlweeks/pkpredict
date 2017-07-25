#' Plot method for PK model
#'
#' @param object An object of class \code{pkm}
#' @param alp Value used to produce (1 - \code{alp})% credible intervals
#' @param ... Additional arguments
#'
#' @return Concentration-time curve for the fitted pharmacokinetic model
#' @export
#'
#' @examples

plot.pkm <- function(object, alp = 0.05, ...){
  col_bg <- "#AAAAB5"
  col_fg <- "#3E465A"

  #############
  est <- object$optim
  ivt <- object$infsched
  cod <- object$cod

  ## Compute plotting times
  ## - ensure peak and trough times
  ## - avoid time zero
  tms <- sapply(ivt, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  tms <- pmax(1e-3, tms)
  # Compute concentration at each time
  sol <- pk_solution(v_1=exp(est$par[1]), k_10=exp(est$par[2]),
                     k_12=exp(est$par[3]), k_21=exp(est$par[4]), ivt=ivt)
  con <- apply(sol(tms)*1000, 2, function(x) pmax(0,x))
  con <- con[1,]

  ## Approximate standard deviation of log concentration-time curve
  grd <- fdGrad(est$par, function(pars) {
    sol <- pk_solution(v_1=exp(pars[1]), k_10=exp(pars[2]),
                       k_12=exp(pars[3]), k_21=exp(pars[4]), ivt=ivt)
    log(sol(tms)[1,]*1000) ## mulitply by 1000: g/l -> ug/ml
  })
  sde <- sqrt(diag(grd %*% solve(-est$hessian) %*% t(grd)))
  sde <- ifelse(is.nan(sde), 0, sde)



  #############

  par(mfrow=c(1,1))
  plot(tms, con, xlab="Time (h)", ylab="Concentration (ug/mL)",
       ylim=c(0, max(exp(log(con)+qnorm(1-alp/2)*sde), na.rm=TRUE)),
       type='n', main="Concentration vs. Time")

  ## Plot 95% credible bands
  polygon(c(tms,rev(tms)),
          c(exp(log(con) + qnorm(1-alp/2)*sde),
            rev(exp(log(con) - qnorm(1-alp/2)*sde))),
          col=col_bg, border=NA)

  ## Plot posterior estimate
  lines(tms, con, lwd=2, col=col_fg)

  ## Plot measured points
  dat <- object$data
  if(nrow(dat) > 0){
    points(dat$time_h, dat$conc_mcg_ml, pch=16)
  }

  ## Create legend
  legend('topleft', c("Predicted", "95% Credible Band", "Measured"),
         lwd=c(2,4,NA), pch=c(NA,NA,16), col=c(col_fg,col_bg,'black'),
         border=NA, bty='n')

  #Plotting elements for mic statistic
  ftmic <- object$ftmic
  thres <- object$thresh
  abline(h = thres, lty = 2)
  legend("topright", bty = 'n',
         legend = c(paste("fT > threshold:", round(ftmic$ftmic, 3)),
                    paste(ifelse(ftmic$mcmc, "Exact", "Approximate"), "95% CI: (", round(ftmic$conf.int[1], 3), ",",
                          round(ftmic$conf.int[2], 3), ")")))
}
