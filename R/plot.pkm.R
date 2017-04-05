#' Plot method for PK model
#'
#' @param object An object of class \code{pkm}
#' @param alp Value used to produce (1 - \code{alp})% credible intervals
#'
#' @return Concentration-time curve for the fitted pharmacokinetic model
#' @export
#'
#' @examples

plot.pkm <- function(object, alp = 0.05){
  col_bg <- "#AAAAB5"
  col_fg <- "#3E465A"

  tms <- object$times
  con <- object$fitted.values
  sde <- object$se_con

  par(mfrow=c(1,1))
  plot(tms, con[1,], xlab="Time (h)", ylab="Concentration (ug/mL)",
       ylim=c(0, max(exp(log(con[1,])+qnorm(1-alp/2)*sde), na.rm=TRUE)),
       type='n', main="Concentration vs. Time")

  ## Plot 95% credible bands
  polygon(c(tms,rev(tms)),
          c(exp(log(con[1,]) + qnorm(1-alp/2)*sde),
            rev(exp(log(con[1,]) - qnorm(1-alp/2)*sde))),
          col=col_bg, border=NA)

  ## Plot posterior estimate
  lines(tms, con[1,], lwd=2, col=col_fg)

  ## Plot measured points
  dat <- object$data
  points(dat$time_h, dat$conc_mcg_ml, pch=16)

  ## Create legend
  legend('topleft', c("Predicted", "95% Credible Band", "Measured"),
         lwd=c(2,4,NA), pch=c(NA,NA,16), col=c(col_fg,col_bg,'black'),
         border=NA, bty='n')

  #Plotting elements for mic statistic
  ftmic <- object$ftmic
  abline(h = thres, lty = 2)
  legend("topright", bty = 'n',
         legend = c(paste("fT > threshold:", round(ftmic[1], 3)),
                    paste("95% CI: (", round(ftmic[2], 3), ",", round(ftmic[3], 3), ")")))
}
