#' Unit conversion
#'
#' Converts units for compatibility with package functions.
#'
#'
#' @param x The value of the number to be converted, or to convert 2 (depending on `forward` argument)
#' @param mvunits Units of the mass and volume of the data. Must be entered as "mass/volume"
#' @param forward Boolean. Tells whether or not the conversion is forward (into g/L) or backward (from g/L to original)
#'
#' @details
#'
#' @return posterior estimates
#' @export
#'
#' @examples
#' # Insert example from Bayes.R
#' pkm(concentration ~ time, dat_d, ivt_d) # something like this

unit_conversion <- function(x, mvunits, forward = TRUE){
  if(1 - grepl("/", mvunits)){stop("Units must be entered as 'mass/volume'")}

  slash_loc <- gregexpr(pattern = '/', mvunits)[[1]]
  mass <- substr(mvunits, 1, slash_loc-1)
  vol <- substr(mvunits, slash_loc+1, nchar(mvunits))

  if(!(tolower(mass) %in% c("mcg", "mg", "kg", "g"))){stop("Mass must be one of 'mcg', 'mg', 'kg', 'g'")}
  if(!(tolower(vol) %in% c("ml", "l"))){stop("Mass must be one of 'mL', 'L'")}

  if(forward){ # convert to g/L units for compatibility with functions
    if(tolower(mass) == "mcg") x <- x*1e-6
    if(tolower(mass) == "mg") x <- x*.001
    if(tolower(mass) == "kg") x <- x*1000

    if(tolower(vol) == "ml") x <- x*1000
  }
  if(1-forward){ # convert from g/L units to original
    if(mass == "mcg") x <- x/1e-6
    if(mass == "mg") x <- x/.001
    if(mass == "kg") x <- x/1000

    if(vol == "mL") x <- x/1000
  }

  return_unit <- ifelse(forward, "g/L", mvunits)

  attr(x, "units") <- return_unit

  return(x)
}

