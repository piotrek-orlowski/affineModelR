#' @title Characteristic function, moments, cumulant-generating function.
#' @name affineCFandDerivs
#' @description Functions to evaluate the P- or Q-measure characteristic function, evaluate CF derivatives with respect to the first argument (log-asset price).
#' @param u \code{U x (N.factors+1)} matrix of points at which the CF and its derivative should be evaluated.
#' @param params.Q Q measure (pricing) parameters. See \code{\link{jumpDiffusionODEs}}.
#' @param params.P P measure (statistical) parameters, optional. See \code{\link{jumpDiffusionODEs}}.
#' @param t.vec \code{T}-length numeric vector of maturities
#' @param v.0 \code{S x N.factors} matrix of volatility factor values
#' @param jumpTransform string, name of function to evaluate the jumpTransform in the model
#' @param N.factors number of stochastic volatility factors, 3 by default
#' @param CGF return Cumulant-Generating Function or Characteristic/Moment Generating Function? CGF if \code{TRUE}. Log of CF if \code{u} is complex.
#' @param ... Further arguments passed to \code{\link{solveODE}}
#' @export affineCF
#' @return \code{affineCF} evaluates the CF/CGF of an affine model under P or Q measures, at matrix \code{u} of size  \code{U x (N.factors+1)}, maturity vector \code{t.vec} of length \code{T}, and variance factor matrix of size \code{S x N.factors}. The result is a \code{U x T x S} matrix. \cr 
#' \code{affineCFderivs} evaluates derivatives of the characteristic function with respect to its first argument via ODE solutions of an extended system. A list of length 4 is returned, each holding an \code{U x T x S} matrix. This is useful for calculating moments of log-returns.

affineCF <- function(u, params.Q, params.P = NULL, t.vec, v.0, jumpTransform = getPointerToJumpTransform(fstr = 'expNormJumpTransform')$TF, N.factors = 3, CGF= FALSE, mod.type = "standard", ...){
  
  # define mkt
  mkt <- data.frame(p=1,q=0,r=0,t=t.vec)
  
  # solve ODEs
  if(is.null(params.P)){
    ode.sol <- jumpDiffusionODEs(u = u, params = params.Q, mkt = mkt, jumpTransform = jumpTransform, N.factors = N.factors, mod.type = mod.type, ...)  
  } else {
    ode.sol <- jumpDiffusionODEsP(u = u, params.P = params.P, params.Q = params.Q, mkt = mkt, jumpTransform = jumpTransform, N.factors = N.factors, mod.type = mod.type, ...)
  }
  
  
  # evaluate CF
  log.cf.val <- apply(v.0,1,function(v.vec){
    v.vec <- c(v.vec,1)
    v.vec <- array(v.vec, dim = dim(ode.sol)[c(3,1,2)])
    v.vec <- aperm(v.vec, c(2,3,1))
    res <- v.vec * ode.sol
    res <- apply(res,c(1,2),sum)
  })
  
  log.cf.val <- array(as.complex(log.cf.val), dim = c(nrow(u),length(t.vec),nrow(v.0)))
  
  # return
  if(CGF){
    return(log.cf.val)
  } else {
    return(exp(log.cf.val))
  }
}

#' @export affineCFderivs
#' @describeIn affineCFandDerivs

affineCFderivs <- function(u, params.Q, params.P = NULL, t.vec, v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 3, mod.type = 'standard', ...){
  
  # define mkt
  mkt <- data.frame(p=1,q=0,r=0,t=t.vec)
  
  # solve ODEs
  ode.sol <- odeExtSolveWrap(u = u, params.Q = params.Q, params.P = params.P, mkt = mkt, N.factors = N.factors, jumpTransform = jumpTransform, mod.type = mod.type, ...)
 
  ode.sol.cf <- ode.sol[,,c(paste0("b",1:N.factors),"a"),drop=FALSE]
  
  # evaluate CF
  log.cf.val <- apply(v.0,1,function(v.vec){
    v.vec <- c(v.vec,1)
    v.vec <- array(v.vec, dim = dim(ode.sol.cf)[c(3,1,2)])
    v.vec <- aperm(v.vec, c(2,3,1))
    res <- v.vec * ode.sol.cf
    res <- apply(res,c(1,2),sum)
  })
  
  log.cf.val <- array(as.numeric(log.cf.val), dim = c(nrow(u),length(t.vec),nrow(v.0)))
  
  # module with first derivative
  ode.sol.d1 <- ode.sol[,,c(paste0("bp1",1:N.factors),"ap1"),drop=FALSE]
  d1.val <- apply(v.0,1,function(v.vec){
    v.vec <- c(v.vec,1)
    v.vec <- array(v.vec, dim = dim(ode.sol.d1)[c(3,1,2)])
    v.vec <- aperm(v.vec, c(2,3,1))
    res <- v.vec * ode.sol.d1
    res <- apply(res,c(1,2),sum)
  })
  
  d1.val <- array(as.numeric(d1.val), dim = c(nrow(u),length(t.vec),nrow(v.0)))
  
  cf.d1 <- exp(log.cf.val)*d1.val
  
  # module with second derivative
  ode.sol.d2 <- ode.sol[,,c(paste0("bp2",1:N.factors),"ap2"),drop=FALSE]
  d2.val <- apply(v.0,1,function(v.vec){
    v.vec <- c(v.vec,1)
    v.vec <- array(v.vec, dim = dim(ode.sol.d2)[c(3,1,2)])
    v.vec <- aperm(v.vec, c(2,3,1))
    res <- v.vec * ode.sol.d2
    res <- apply(res,c(1,2),sum)
  })
  
  
  d2.val <- array(as.numeric(d2.val), dim = c(nrow(u),length(t.vec),nrow(v.0)))
  
  cf.d2 <- exp(log.cf.val)*(d1.val^2 + d2.val)

  # module with third derivative
  ode.sol.d3 <- ode.sol[,,c(paste0("bp3",1:N.factors),"ap3"),drop=FALSE]
  d3.val <- apply(v.0,1,function(v.vec){
    v.vec <- c(v.vec,1)
    v.vec <- array(v.vec, dim = dim(ode.sol.d3)[c(3,1,2)])
    v.vec <- aperm(v.vec, c(2,3,1))
    res <- v.vec * ode.sol.d3
    res <- apply(res,c(1,2),sum)
  })
  
  d3.val <- array(as.numeric(d3.val), dim = c(nrow(u),length(t.vec),nrow(v.0)))
  
  cf.d3 <- exp(log.cf.val)*(d1.val^3 + 3*d1.val*d2.val + d3.val)
  
  dimnames(log.cf.val) <- dimnames(cf.d3) <- dimnames(cf.d2) <- dimnames(cf.d1) <- list(paste0("u.",1:nrow(u)),t.vec,paste0("v.",1:nrow(v.0)))
  
  return(list(cf = exp(log.cf.val), cf.d1 = cf.d1, cf.d2 = cf.d2, cf.d3 = cf.d3))
}

#' @export affineCFderivsNumerical
#' @describeIn affineCFandDerivs

affineCFderivsNumerical <- function(u, params.Q, params.P = NULL, t.vec, v.0, N.factors = 3, hh = 1e-4, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type, ...){
  
  uu <- seq(-5*hh,5*hh,by=hh)
  u <- apply(u,1,function(u.row){
    u.mat <- matrix(data = u.row, nrow = length(uu), ncol = length(u.row), byrow = TRUE)
    u.mat[,1] <- u.mat[,1]+uu
    return(u.mat)
  })
  
  u <- matrix(u, nrow = length(u)/(N.factors+1), ncol = N.factors+1, byrow= FALSE)
  cf.for.diff <- affineCF(u = u, params.Q = params.Q, params.P = params.P, t.vec = t.vec, v.0 = v.0, jumpTransform = jumpTransform, N.factors = dim(v.0)[2], CGF = FALSE, mod.type = mod.type)
  cf.d1 <- apply(X = cf.for.diff, MARGIN = c(2,3), FUN = function(cff){
    f0 <- cff[median(order(cff))]
    fVal <- cff[-median(order(cff))]
    fVal <- fVal[-c(1,10)]
    deriv <- fdDerivative(f0 = f0, fVal = fVal, p = 1, k = 7, h = hh)
    return(deriv)
  })
  
  cf.d2 <- apply(X = cf.for.diff, MARGIN = c(2,3), FUN = function(cff){
    f0 <- cff[median(order(cff))]
    fVal <- cff[-median(order(cff))]
    fVal <- fVal[-c(1,10)]
    deriv <- fdDerivative(f0 = f0, fVal = fVal, p = 2, k = 6, h = hh)
    return(deriv)
  })
  
  cf.d3 <- apply(X = cf.for.diff, MARGIN = c(2,3), FUN = function(cff){
    f0 <- cff[median(order(cff))]
    fVal <- cff[-median(order(cff))]
    fVal <- fVal[-c(1,10)]
    deriv <- fdDerivative(f0 = f0, fVal = fVal, p = 3, k = 5, h = hh)
    return(deriv)
  })
  
  return(list(cf.d1 = cf.d1, cf.d2 = cf.d2, cf.d3= cf.d3))
}
