#' Predict method for PK model
#'
#' @param object Object of class "pkm"
#' @param newivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose for which predictions will be obtained
#'
#' @return
#' @export
#'
#' @examples
#'

predict.pkm <- function(object, newivt){

  # Parameter estimates from fitted pkm model
  est <- object$est

  tms <- sapply(newivt, function(x) c(x$begin, x$end))
  tms <- c(tms, max(tms)+cod)
  tms <- unlist(sapply(1:(length(tms)-1), function(i) {
    s1 <- seq(tms[i], tms[i+1], 1/10)
    if(tms[i+1] %% 1/10)
      s1 <- c(s1, tms[i+1])
    return(s1)
  }))
  tms <- pmax(1e-3, tms)

  grd <- fdGrad(est$par, function(pars) {
    sol <- pk_solution(v_1=exp(pars[1]), k_10=exp(pars[2]),
                       k_12=exp(pars[3]), k_21=exp(pars[4]), ivt=newivt)
    log(sol(tms)[1,]*1000) ## mulitply by 1000: g/l -> ug/ml
  })
  sde <- sqrt(diag(grd %*% solve(-est$hessian) %*% t(grd)))
  sde <- ifelse(is.nan(sde), 0, sde)

  ## Posterior estiamte
  sol <- pk_solution(v_1=exp(est$par[1]), k_10=exp(est$par[2]),
                     k_12=exp(est$par[3]), k_21=exp(est$par[4]), ivt=newivt)
  con <- apply(sol(tms)*1000, 2, function(x) pmax(0,x))

  # Predicted ft>mic estimate
  ftmic <- mic_stat(pars = est$par, newivt, object$data,
                    tms, con[1,], th = object$thresh)

  res <- list("predict" = data.frame("conc" = con, "tms" = tms),
              "infsched" = newivt,
              "pftmic" = confint(ftmic))
  return(res)
}
