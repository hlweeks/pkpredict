% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_prior.R
\name{log_prior}
\alias{log_prior}
\title{Log-prior for pharmacokinetic parameters}
\usage{
log_prior(lpr, mu, sig, ler_mean, ler_sdev)
}
\arguments{
\item{lpr}{log-PK parameter vector with error: (lv_1, lk_10, lk_12, lk_21, ler_mean)}

\item{mu}{log-mean of the PK parameter distribution}

\item{sig}{log-variance-covariance matrix of the PK parameter distribution}

\item{ler_mean}{log-mean of the error distribution}

\item{ler_sdev}{log-standard deviation of the error distribution}
}
\value{
The value of the log-prior evaluated at the specified log-parameter vector
}
\description{
Evaluates the prior distribution (normal) of the log parameter vector with
normally distributed error.
}
\details{
Parameter distributions:
log(lv_1, lk_10, lk_12, lk_21) ~ N(mu, sig)
log(error) ~ N(ler_mean, ler_sdev)
}
\examples{
# Default parameter values based on prior study
lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7)
lpk_vcov_d <- 300 * matrix(c(   0.00167,  -0.00128,      0,      0,
                               -0.00128,   0.00154,      0,      0,
                                      0,         0, .00015,      0,
                                      0,         0,      0, .00015), 4, 4)
ler_mean_d <- 2.33
ler_sdev_d <- 0.32

lpr_d  <- c(lpk_mean_d, ler_mean_d)
log_prior(lpr = lpr_d, mu = lpk_mean_d, sig = lpk_vcov_d,
          ler_mean = ler_mean_d, ler_sdev = ler_sdev_d)
}
