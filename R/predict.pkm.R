#' Predict method for PK model
#'
#' @param object Object of class "pkm"
#' @param times Times at which predictions of blood concentration measurements are desired
#'
#' @return Predicted blood concentration measurements
#' @export
#'
#' @examples
#'

predict.pkm <- function(object, times){

  # Parameter estimates from fitted pkm model
  est <- object$est
  sde <- object$fitted.values$se_con


  ## Posterior estiamte
  con <- apply(object$sol.eqn(times)*1000, 2, function(x) pmax(0,x))

  # # Predicted ft>mic estimate
  # ftmic <- mic_stat(pars = est$par, newivt, object$data,
  #                   tms, con[1,], th = object$thresh)


  return(data.frame("conc" = con, "tms" = times))
}
