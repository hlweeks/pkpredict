#' Two-compartment model solution
#'
#' Calculates the solution of a two-comparment model.
#'
#' @param k_10 Elimination rate from central compartment
#' @param k_12 Elimination rate from central to peripheral compartment
#' @param k_21 Elimination rate from peripheral to central compartment
#' @param v_1 Volume of central compartment
#' @param k_R Rate of infusion (g/h)
#' @param c_0 Initial concentrations of c(central, peripheral) compartments
#'
#' @details
#' The 'pk_basic_solution' function implements the
#' two-compartment model solution for fixed parameters k_10,
#' k_12, k_21, v_1, k_R, and c_0 (initial concentrations).
#' This function is based on the solution to a nonhomogeneous
#' linear system of ODEs. The basic solution could also be
#' computed using numerical integration, which is commonly used
#' for more complicated models.
#'
#' @return List containing concentration in central and peripheral compartments, respectively.
#'
#' @export

pk_basic_solution <-
  function(k_10, k_12, k_21, v_1, k_R, c_0=c(0,0)) {

    K <- matrix(c(-(k_10+k_12), k_21, k_12, -k_21),
                2,2, byrow=TRUE)
    T_K <- -(k_10 + k_12 + k_21)
    D_K <- (k_10+k_12)*k_21 - k_12*k_21
    lambda <- T_K/2 + c(-1,1)*(T_K^2/4 - D_K)^(1/2)
    N_21 <- (k_21+lambda[2])/k_21
    N_22 <- (k_21+lambda[1])/k_21
    N <- matrix(c(-1,N_21,-1,N_22),2,2)
    gamma <- c(k_R/k_10/v_1, k_12*k_R/k_10/k_21/v_1)
    r <- c(( N_22*(c_0[1]-(k_R/v_1/k_10)) +
               (c_0[2]-(k_R/v_1/k_10)*k_12/k_21)),
           (-N_21*(c_0[1]-(k_R/v_1/k_10)) -
              (c_0[2]-(k_R/v_1/k_10)*k_12/k_21)))
    r <- r * k_21/(diff(lambda))
    c_1 <- function(t)
      gamma[1] + r[1]*N[1,1]*exp(lambda[1]*t) +
      r[2]*N[1,2]*exp(lambda[2]*t)
    c_2 <- function(t)
      gamma[2] + r[1]*N[2,1]*exp(lambda[1]*t) +
      r[2]*N[2,2]*exp(lambda[2]*t)

    return(list(c_1 = c_1, c_2 = c_2))
  }
