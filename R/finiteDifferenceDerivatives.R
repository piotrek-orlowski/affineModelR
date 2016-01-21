#' @title Finite difference coefficients.
#' @description Calculation of coefficients for finite-difference approximation of derivative of order p with error of order O(h^(p+k)), given a vector of function values.
#' @param p Derivative order.
#' @param k Desired accuracy.
#' @param h Grid step size.
#' @param f0 Value of function at differentiation point.
#' @param fVal Values of function at points that will be used for approximation.
#' @export
#' @return \code{fdCoefficients} Vector of coefficients for linear combination of function values that will give the desired finite-difference approximation. \cr \cr \code{fdDerivative} numerical derivative based on vector of function values.
#' @details If p+k is odd, the finite difference scheme will be asymmetric, with one more point on the negative side.

fdCoefficients <- function(p,k){
  # number of equations
  M <- p+k
  
  m <- c(seq(floor(-M/2),-1),seq(1,floor(M/2)))
  A <- t(apply(matrix(1:M),1,function(pp) m^pp ))
  b <- rep(0,M)
  b[p] <- 1
  
  DCoef <- solve(A,b)
  DCoef <- DCoef
  
  return(DCoef)
}

#' @export fdDerivative

fdDerivative <- function(f0,fVal,p,k,h){
  DCoef <- fdCoefficients(p,k)
  deriv <- sum(DCoef * fVal) - sum(DCoef)*f0
  deriv <- deriv / h^p * factorial(p)
  return(deriv)
}