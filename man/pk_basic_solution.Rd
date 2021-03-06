% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_basic_solution.R
\name{pk_basic_solution}
\alias{pk_basic_solution}
\title{Two-compartment model solution}
\usage{
pk_basic_solution(k_10, k_12, k_21, v_1, k_R, c_0 = c(0, 0))
}
\arguments{
\item{k_10}{Elimination rate from central compartment}

\item{k_12}{Elimination rate from central to peripheral compartment}

\item{k_21}{Elimination rate from peripheral to central compartment}

\item{v_1}{Volume of central compartment}

\item{k_R}{Rate of infusion (g/h)}

\item{c_0}{Initial concentrations of c(central, peripheral) compartments}
}
\value{
List containing concentration in central and peripheral compartments, respectively.
}
\description{
Calculates the solution of a two-comparment model.
}
\details{
The 'pk_basic_solution' function implements the
two-compartment model solution for fixed parameters k_10,
k_12, k_21, v_1, k_R, and c_0 (initial concentrations).
This function is based on the solution to a nonhomogeneous
linear system of ODEs. The basic solution could also be
computed using numerical integration, which is commonly used
for more complicated models.
}
