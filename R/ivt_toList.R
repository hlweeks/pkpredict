#' Convert infusion schedule characteristics to list format
#'
#' Helper function: Given vectors that characterize the infusion schedule, this helper function converts the vectors
#' to list format for compatibility with other package functions
#'
#' @param begin Times at which each infusion began (relative to start of first infusion)
#' @param end Times at which each infusion ended (relative to start of first infusion)
#' @param dur Duration of each dose (h)
#' @param rate Rate of infusion for each dose (g/h)
#'
#' @return A list characterizing the infusion schedule
#'
#' @export
#'
#' @examples
#' b <- c(0, 8, 16, 24)
#' e <- c(0.5, 8.5, 16.5, 24.5) # alternatively, specify durations: d <- rep(0.5, 4)
#' r <- rep(6, 4)
#' ivt_toList(begin = b, end = e, rate = r)

ivt_toList <- function(begin, end = NULL, dur = NULL, rate){
  if(is.null(end) & is.null(duration)){
    stop('At least one of end or dur must be specified')
  }
  if(!is.null(end) & !is.null(duration)){
    stop('Only one of end or dur can be specified')
  }

  val <- ifelse(is.null(end), begin + dur, end)

  df <- data.frame(begin, val, rate) # Will produce error is args are not of the same length

  ivt <- apply(df, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))

  return(ivt)
}
