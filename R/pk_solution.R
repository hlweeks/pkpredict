#' PK solution over dosing cycle
#'
#' The 'pk_solution' function impelments a piecewise solution
#' to the two-compartment model, taking into account the sequence
#' of dosing events:
#'
#' @param k_10 Elimination rate from central compartment
#' @param k_12 Elimination rate from central to peripheral compartment
#' @param k_21 Elimination rate from peripheral to central compartment
#' @param v_1 Volume of central compartment
#' @param ivt List with containing start of infusion times, end of infusion times,
#' and rate of infusion at each dose
#' @param init Initial concentration of drug in each compartment
#'
#' @return A function that will compute concentration of drug in the body given time since
#' first infusion
#'
#' @export

pk_solution <-
  function(k_10, k_12, k_21, v_1, ivt, init=c(0,0)) {
    ## create a list of event times
    ibe <- sapply(ivt, `[[`, 'begin')
    ied <- sapply(ivt, `[[`, 'end')
    prd <- sort(unique(c(0, c(ibe,ied), Inf)))
    rits <- list()
    ## compute basic solution in each interval
    for(i in 1:(length(prd)-1)) {
      civt <- sapply(ivt, function(iv) {
        if(prd[i] >= iv$begin && prd[i] < iv$end) {
          iv$k_R
        } else { 0 }
      })

      rit <- list(begin=prd[i], end=prd[i+1],
                  idose=sum(civt))
      if(i == 1) {
        rit$init <- init
      } else {
        rit$init <- c(rits[[i-1]]$c_1(rits[[i-1]]$end-rits[[i-1]]$begin),
                      rits[[i-1]]$c_2(rits[[i-1]]$end-rits[[i-1]]$begin))
      }

      sol <- pk_basic_solution(k_10, k_12, k_21, v_1,
                               k_R=rit$idose, c_0=rit$init)

      rits[[i]] <- c(rit, sol)
    }

    # both c_1 and c_2
    ret <- function(tms) {
      sapply(tms, function(t) {
        val <- NA
        for(rit in rits) {
          if(t >= rit$begin && t <= rit$end) {
            val <- c(rit$c_1(t-rit$begin),
                     rit$c_2(t-rit$begin))
            break
          }
        }
        val
      })
    }
    return(ret)
  }
