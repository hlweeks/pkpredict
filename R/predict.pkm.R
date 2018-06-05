#' Predict method for PK model
#'
#' @param object Object of class "pkm"
#' @param times Times at which predictions of blood concentration measurements are desired
#'
#' @return Predicted blood concentration measurements
#' @export
#'
#'

predict.pkm <- function(object, times = NULL, ...){

  # Parameter estimates from fitted pkm model
  est <- object$est
  sde <- object$fitted.values$se_con


  ## Posterior estiamte
  if(is.null(times)){
    return(data.frame("Time" = object$fitted.values$time,
                      "Concentration" = object$fitted.values$conc))
  }else{
    con <- apply(object$sol.eqn(times)*1000, 2, function(x) pmax(0,x))
    return(data.frame("Time" = times,
                      "Concentration" = con))
  }

  # # Predicted ft>mic estimate
  # ftmic <- mic_stat(pars = est$par, newivt, object$data,
  #                   tms, con[1,], th = object$thresh)
}
