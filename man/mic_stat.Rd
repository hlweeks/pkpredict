% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mic_stat.R
\name{mic_stat}
\alias{mic_stat}
\title{Time Above Pharmacokinetic Threshold}
\usage{
mic_stat(ivt, th, dat = data.frame(), pars = c(lv_1 = 3.223, lk_10 = -1.65,
  lk_12 = -7, lk_21 = -7, lerr = 2.33), mu = c(lv_1 = 3.223, lk_10 = -1.65,
  lk_12 = -7, lk_21 = -7), sig = 300 * matrix(c(0.00167, -0.00128, 0, 0,
  -0.00128, 0.00154, 0, 0, 0, 0, 0.00015, 0, 0, 0, 0, 0.00015), 4, 4),
  ler_mean = 2.33, ler_sdev = 0.32, cod = 12, timeint = c(0,
  max(sapply(ivt, function(x) x$end)) + cod), conf.level = 0.95,
  mcmc = FALSE, nreps = 5000, nburnin = 2000, seed = NULL,
  shiny = FALSE)
}
\arguments{
\item{ivt}{List with containing start of infusion times (h), end of infusion times (h),
and rate of infusion (g/h) at each dose}

\item{th}{Threshold value for effective treatment - check units}

\item{dat}{data}

\item{pars}{Vector of pharmacokinetic parameters of length 5: (log(v_1), log(k_10), log(k_12), log(k_21), log(err))}

\item{mu}{prior pk param mean}

\item{sig}{prior pk param vcov matrix}

\item{ler_mean}{prior error dist mean}

\item{ler_sdev}{prior error dist sd}

\item{cod}{Length of time after end of last dose to consider}

\item{timeint}{Vector with two elements indicating the start and end of the time interval over which
the statistic is computed. By default covers the full dosing period + \code{cod}}

\item{conf.level}{Desired confidence level of the interval}

\item{mcmc}{logical: should estimate of time above threshold be computed using MCMC (false = laplace approximation)}

\item{nreps}{mcmc replications}

\item{nburnin}{mcmc burn in iterations}

\item{seed}{seed for replication}

\item{shiny}{is shiny being used}
}
\value{
Fraction of time spent above the specified threshold from time of first dose through
\code{cod} hours after end of the last dose
}
\description{
Compute an estimate of pharmacodynamic target attainment, typically defined
to be some multiple of the minimum inhibitory concentration of the infecting microorganism.
}
