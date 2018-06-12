#' Update method for PK model
#'
#' @param object An object of class "pkm
#' @param formula A formula where the left side is the measured concentration of drug
#' and the right side is the times of concentration measurements
#' @param newdat Data frame with concentration data (time of measurement in hours and concentration in mcg/ml)
#' @param newivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose for which predictions will be obtained
#' @param ... other params
#'
#' @details Update a fitted PK model with either new infusion information, new obseved concentration data,
#' or both.
#'
#' @return Updated PK model
#' @export
#'
#'

update.pkm <- function(object, formula = NULL, newdat = NULL, newivt = NULL, ...){
  if(is.null(newivt) & is.null(newdat)){stop("At least one of newivt, newdat must be specified")}

  dta <- object$data
  ivt <- object$infsched

  if(!is.null(newdat)){
    if(is.null(formula)){stop("concentration ~ time formula not specified")}
    dat_add <- model.frame(formula, dta)

    dta <- cbind(dta, dat_add)
  }

  if(!is.null(newivt)){
    #Need some check of same formatting - see .Rmd for comment on this
    ivt <- cbind(ivt, newivt)
  }

  pkm(formula = conc_mcg_ml ~ time_h,
      data = dta, ivt = ivt, pars = object$est,
      alp = object$alpha, cod = object$cod, thres = object$thresh)
}
