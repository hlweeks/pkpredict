#' Unit conversion
#'
#' Converts units for compatibility with package functions.
#'
#' @param x The value of the number to be converted, or to convert to (depending on `to_gL` argument)
#' @param mvunits Units of the mass and volume of the data. Must be entered as "mass/volume"
#' @param to_gL Boolean. Tells whether or not the conversion is into g/L or bfrom g/L to original
#'
#' @return Value equivalent to input with the appropriate units.
#'
#' @export
#'
#' @examples
#' unit_conversion(64, "mcg/mL")

unit_conversion <- function(x, mvunits, to_gL = TRUE){
  if(1 - grepl("/", mvunits)){stop("Units must be entered as 'mass/volume'")}

  slash_loc <- gregexpr(pattern = '/', mvunits)[[1]]
  if(length(slash_loc) > 1){stop("Units must be entered as 'mass/volume'")}

  mass <- substr(mvunits, 1, slash_loc-1)
  vol <- substr(mvunits, slash_loc+1, nchar(mvunits))

  if(!(tolower(mass) %in% c("ng", "mcg", "mg", "kg", "g"))){
    stop("Mass must be one of 'ng', 'mcg', 'mg', 'g', 'kg'")
  }
  if(!(tolower(vol) %in% c("ml", "l"))){stop("Mass must be one of 'mL', 'L'")}


  if(class(to_gL) != "logical"){stop("`forward` must be T/F")}

  if(to_gL){ # convert to g/L units for compatibility with functions
    if(tolower(mass) == "ng") x <- x*1e-9
    if(tolower(mass) == "mcg") x <- x*1e-6
    if(tolower(mass) == "mg") x <- x*.001
    if(tolower(mass) == "kg") x <- x*1000

    if(tolower(vol) == "ml") x <- x*1000
  }
  if(1-to_gL){ # convert from g/L units to original
    if(tolower(mass) == "ng") x <- x/1e-9
    if(tolower(mass) == "mcg") x <- x/1e-6
    if(tolower(mass) == "mg") x <- x/.001
    if(tolower(mass) == "kg") x <- x/1000

    if(tolower(vol) == "mL") x <- x/1000
  }

  return_unit <- ifelse(to_gL, "g/L", mvunits)

  attr(x, "units") <- return_unit

  return(x)
}

