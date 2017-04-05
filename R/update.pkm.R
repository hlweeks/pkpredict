#' Update method for PK model
#'
#' @param object An object of class "pkm
#' @param newivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose for which predictions will be obtained
#' @param formula A formula where the left side is the measured concentration of drug
#' and the right side is the times of concentration measurements
#' @param data Data frame with concentration data (time of measurement in hours and concentration in mcg/ml)
#'
#' @return Updated PK model
#' @export
#'
#' @examples
#'

update.pkm <- function(object, formula, newdat = NULL, newivt = NULL){
  if(is.null(newivt) & is.null(newdat)){stop("At least one of newivt, newdat must be specified")}

  ivt <- ifelse(is.null(newivt), object$infsched, newivt)
  dat <- ifelse(is.null(newdat), object$data, newdat)

}
