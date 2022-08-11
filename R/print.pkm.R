#' Print method for pkm object
#'
#' @param x pkm model
#' @param digits number of digits to round to
#' @param ... other parameters
#'
#' @return Printed model
#' @export
#'

print.pkm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")


  ivt_print <- as.data.frame(t(sapply(x$infsched, unlist)))
  colnames(ivt_print) <- c("Start (h)", "End (h)", paste0("Inf Rate (", x$units, ")"))

  cat("Infusion Schedule:\n")
  print.data.frame(ivt_print, row.names=F)

  cat("\nft > k x MIC estimate:\n",
      paste0(round(x$ftmic$ftmic, digits), " (",
             round(x$ftmic$conf.int[1], digits), ", ",
             round(x$ftmic$conf.int[2], digits), ")"))

  cat("\n")
  invisible(x)
}
