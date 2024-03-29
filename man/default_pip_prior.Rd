% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_pip_prior.R
\docType{data}
\name{default_pip_prior}
\alias{default_pip_prior}
\title{Default prior specification for piperacillin}
\format{
A list with 4 elements specifying the default prior distribution for function arguments: 
log_pk_mean, log_pk_vcov, log_err_mean, log_err_sdev
}
\usage{
default_pip_prior
}
\description{
The prior values stored here were derived from random effects estimates from a prior pharmacokinetic
study using Gaussian kernel density estimation. See vignette for details. This prior specification is 
also stored as an option when the package is loaded and can be retreived using 
`getOption("pkpredict.pip.default.prior")`.
}
\examples{
data(default_pip_prior)
}
\keyword{datasets}
