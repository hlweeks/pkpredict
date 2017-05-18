#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

print.pkm <- function(object, ...){
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(object), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}