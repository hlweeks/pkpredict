#' Convert infusion schedule characteristics to list format
#'
#' Helper function: Given vectors that characterize the infusion schedule, this helper function converts the vectors
#' to list format for compatibility with other package functions
#'
#' @param begin Times at which each infusion began (relative to start of first infusion)
#' @param end Times at which each infusion ended (relative to start of first infusion)
#' @param dur Duration of each dose
#' @param rate Rate of infusion for each dose
#'
#' @return A list characterizing the infusion schedule
#'
#' @examples
#' # Using begin/end times
#' b <- c(0, 8, 16, 24)
#' e <- c(0.5, 8.5, 16.5, 24.5)
#' r <- rep(6, 4)
#' ivt_toList(begin = b, end = e, rate = r)
#'
#' Using begin time with infusion duration
#' b <- c(0, 8, 16, 24)
#' d <- 0.5 # Or d <- rep(0.5, 4)
#' r <- 6
#' ivt_toList(begin = b, dur = d, rate = r)
#'

ivt_toList <- function(begin, end = NULL, dur = NULL, rate){

  if(is.null(end) & is.null(dur)){
    stop('At least one of end or dur must be specified')
  }

  if(!is.null(end) & !is.null(dur)){
    stop('Only one of end or dur can be specified')
  }

  if(!is.null(dur)){
    if(length(dur) == 1){
      dur <- rep(dur, length(begin))
    }else if(length(dur != length(begin))){
      stop('`dur` must have length = 1 or length = length(begin)')
    }
  }


  if(length(rate) == 1){
    rate <- rep(rate, length(begin))
  }else if(length(rate != length(begin))){
    stop('`rate` must have length = 1 or length = length(begin)')
  }


  if(!is.null(end) & (length(end) != length(begin))){
    stop('`begin` and `end` must have the same length')
  }


  end <- if(is.null(end)){begin + dur}else{end}

  df <- data.frame(begin, end, rate)
  names(df) <- NULL

  ivt <- apply(df, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))

  return(ivt)
}
