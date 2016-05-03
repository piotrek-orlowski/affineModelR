#' @title ODE solver for jump diffusion models.
#' @rdname jumpDiffODEs
#' @description Return the solution of the ODE describing the conditional Laplace transform of the affine states (log stock + volatility factors). Calculations performed under the pricing measure Q, i.e. with value of Laplace transform constrained to 1 for \code{u = c(1,rep(0,N.factors))}, where \code{N.factors} is the number of stochastic volatility factors, or under the statistical measure P, with Q parameters required to set risk premia and dynamics.
#' @param u \code{Ux(N.factors+1)} matrix of complex numbers.
#' @param params A list, containing the structures defining the stochastic volatitilites and jumps. See Details.
#' @param params.P P (statistical measure) parameters, see Details.
#' @param params.Q Q (pricing measure) parameters, see Details.
#' @param mkt data.frame describing the market structure (times to maturity, interest rate, dividend yield, current stock price). See Details.
#' @param jumpTransform string indicating which jump transform to contain in the specification. Available values: \code{'expNormJumpTransform'} for the co-jump specification from DPS (2000) ``concrete examples'' section, \code{'kouExpJumpTransform'} for a similar co-jump specification with double exponential jumps in the log-asset price.
#' @param rtol relative tolerance for the ODE solver functions, \code{\link[deSolve]{vode}} and \code{\link[deSolve]{zvode}}.
#' @param atol absolute tolerance for the ODE solver functions, \code{\link[deSolve]{vode}} and \code{\link[deSolve]{zvode}}.
#' @param mf Integration method to use in the ODE solver functions, \code{\link[deSolve]{vode}} and \code{\link[deSolve]{zvode}}.
#' @param N.factors The number of volatility factors, of which the first one can co-jump with the stock.
#' @param mod.type string indicating whether the \code{'standard'} DPS (2000) specification should be used, or a variation as seen in Duffie, Pedersen and Singleton (2003), where the volatility of volatility of the first factor can be driven by other factors.
#' @details \code{jumpDiffusionODEs} solves the Riccati equations associated with a jump-diffusion model under the restrictions of the pricing measure: the stock price process is a martingale. \cr
#' \code{jumpDiffusionODEsP} uses \code{params.Q} and \code{params.P} to construct a model with well-defined risk premia and solves the Riccati equation under the statistical measure; this construction is required because the parameters of the stock equation under the statistical measure depend, in some specifications, on the parameters of the stock and volatility equations under the pricing measure.
#' \cr \cr
#' \code{params}, \code{params.P} and \code{params.Q} are \code{list}s with \code{N.factors+1} fields. \cr \cr
#' The fields must be named \code{'1'} through \code{'N.factors'} and \code{'jmp'}. The numbered fields correspond to volatility factor specifications with parameters: \code{kpp} (mean reversion speed), \code{lmb} (volatility-of-volatility), \code{rho} (correlation with log-asset price, leverage effect), \code{eta} (long-run mean), and \code{phi}: volatility scaling in the log-asset price equation. \code{kpp}, \code{lmb} greater than 0, \code{phi} greater than or equal to 0. If \code{phi = 0}, the volatility factor does not drive the log-asset price, but can potentially drive the jump intensity; either \code{eta} or \code{phi} should be normalised to \code{1}. The parameter \code{params.P$[["1"]]$erp} sets how the risk premium is driven by the total level of volatility. \cr \cr
#' The list \code{$jmp} has fields: \code{lvec} (double): jump intensity, constant part, \code{lprop} jump intensity loadings on volatility factors, \code{length(lprop) = N.factors}, \code{muSc} (positive, mean vol jump for 'expNormJumpTransform' or one over mean vol jump for 'kouExpJumpTransform'), \code{rhoc} (jump leverage), \code{muYc} (jump distribution location parameter), \code{sigmaYc} (jump distribution scale parameter).
#' \cr \cr
#' \code{mkt} data.frame describing maturities and corresponding interest rates, dividend yields. Fields: \code{p}: initial stock price, normalized to 1 without loss of generality, \code{q} dividend yield per annum, \code{r} interest rate, per annum, \code{t} maturity for which the ODE solutions are to be calculated.
#' @useDynLib affineModelR
#' @export jumpDiffusionODEs
#' @return From \code{jumpDiffusionODEs} and \code{jumpDiffusionODEsP}: An array of size \code{UxTx(N.factors + 1)} where \code{U = nrow(u)}, \code{T = length(mkt$t)} (number of maturities), code{N.factors+1} is the number of coefficients in the exponentially affine characteristic function. \cr 
#' From \code{odeEstSolveWrap}: an array of UxTx(4x(N.factors+1)): affine coefficients and their derivatives with respect to \code{u[,1]}. This allows to for highly accurate evaluation of the derivatives of the characteristic function of the log-asset price.

jumpDiffusionODEs <- function(u,params,mkt,jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, rtol=1e-12, atol=1e-30, mf = 22, N.factors = 3, mod.type = "standard") {
  
  # sanity checks. Make sure we conform with "new" setup
  stopifnot(ncol(u)== (N.factors+1))
  
  ode.structs <- ODEstructs(params,jumpTransform,mkt,N.factors,mod.type = mod.type)
  
  solMat <- solveODE(u, mkt, ode.structs$K0, ode.structs$K1, ode.structs$l0, ode.structs$l1, ode.structs$H1, jmp = params$jmp,jumpTransform = jumpTransform, mf = mf, rtol=rtol, atol=atol, N.factors = N.factors)
  return(solMat)  
}

#' @describeIn jumpDiffODEs
#' @export jumpDiffusionODEsP

jumpDiffusionODEsP <- function(u,params.P,params.Q,mkt,jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, rtol=1e-13, atol=1e-30, mf = 22, N.factors = 3, mod.type = "standard") {
  
  # sanity checks. Make sure we conform with "new" setup
  stopifnot(ncol(u)== (N.factors+1))
  
  # check measure consistency
  for (nn in 1:N.factors) {
    stopifnot(params.P[[as.character(nn)]]$rho == params.Q[[as.character(nn)]]$rho & params.P[[as.character(nn)]]$lmb == params.Q[[as.character(nn)]]$lmb & params.P[[as.character(nn)]]$phi == params.Q[[as.character(nn)]]$phi)
  }
  
  # if orthogonal erp not defined, then let's assume it is 0
  if (is.null(params.P[[as.character(1)]]$erp)) {
    params.P[[as.character(1)]]$erp <- 0
  }
  if (is.null(params.P[[as.character(1)]]$erp.0)) {
    params.P[[as.character(1)]]$erp.0 <- 0
  }
  
  ode.structs.P <- ODEstructs(params.P,jumpTransform,mkt,N.factors,mod.type)
  ode.structs.Q <- ODEstructs(params.Q,jumpTransform,mkt,N.factors,mod.type)
  
  # now add Q drift to P structs (only the first values, since we need only to have the stock equation)
  ode.structs.P$K1[1,] <- ode.structs.Q$K1[1,]
  ode.structs.P$K0[1] <- ode.structs.Q$K0[1]
  
  # correct stock part
  for (nn in 1:N.factors) {
    phi <- params.P[[as.character(nn)]]$phi
    rho <- params.P[[as.character(nn)]]$rho
    lmb <- params.P[[as.character(nn)]]$lmb
    
    kpp.P <- params.P[[as.character(nn)]]$kpp
    kpp.Q <- params.Q[[as.character(nn)]]$kpp
    eta.P <- params.P[[as.character(nn)]]$eta
    eta.Q <- params.Q[[as.character(nn)]]$eta
    
    ode.structs.P$K1[1,1+nn] <- ode.structs.P$K1[1,1+nn] +  phi * rho * (kpp.Q - kpp.P) / lmb[1] 
    ode.structs.P$K1[1,1+nn] <- ode.structs.P$K1[1,1+nn] + phi * sqrt(1-rho^2) * params.P[[as.character(1)]]$erp
    # add the constant part of the erp that is due to the vrp
    ode.structs.P$K0[1] <- ode.structs.P$K0[1] + phi * rho * (kpp.P * eta.P - kpp.Q * eta.Q) / lmb[1]
  }
  # add identically constant part of the erp
  ode.structs.P$K0[1] <- ode.structs.P$K0[1] + params.P[[as.character(1)]]$erp.0
  
  # correct coefficients in K1 for drift from other factors in the cascade model
  if(mod.type == "cascade.vol"){
    kpp.P <- params.P[["1"]]$kpp
    kpp.Q <- params.Q[["1"]]$kpp
    ode.structs.P$K1[2,3:(3+length(params.Q[["1"]]$lmb)-2)] <- - (kpp.P - kpp.Q)/params.Q[["1"]]$lmb[1] * params.Q[["1"]]$lmb[-1]
  }
  
  solMat <- solveODE(u, mkt, ode.structs.P$K0, ode.structs.P$K1, ode.structs.P$l0, ode.structs.P$l1, ode.structs.P$H1, jmp = params.P$jmp, jumpTransform = jumpTransform, mf = mf, rtol=rtol, atol=atol, N.factors=N.factors)
  return(solMat)  
}

#' @describeIn jumpDiffODEs
#' @export odeExtSolveWrap

odeExtSolveWrap <- function(u, params.Q, params.P = NULL, mkt, rtol = 1e-12, atol = 1e-30, mf = 12, N.factors = 3, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), mod.type = 'standard', ...){
  # sanity checks. Make sure we conform with "new" setup
  stopifnot(ncol(u)== (N.factors+1))
  
  if(!is.null(params.P)){
    # check measure consistency
    for (nn in 1:N.factors) {
      stopifnot(params.P[[as.character(nn)]]$rho == params.Q[[as.character(nn)]]$rho & params.P[[as.character(nn)]]$lmb == params.Q[[as.character(nn)]]$lmb & params.P[[as.character(nn)]]$phi == params.Q[[as.character(nn)]]$phi)
    }
    
    # if orthogonal erp not defined, then let's assume it is 0
    if (is.null(params.P[[as.character(1)]]$erp)) {
      params.P[[as.character(1)]]$erp <- 0
    }
    if (is.null(params.P[[as.character(1)]]$erp.0)) {
      params.P[[as.character(1)]]$erp.0 <- 0
    }
    
    ode.structs.P <- ODEstructs(params.P,jumpTransform,mkt,N.factors,mod.type)
    ode.structs.Q <- ODEstructs(params.Q,jumpTransform,mkt,N.factors,mod.type)
    
    # now add Q drift to P structs (only the first values, since we need only to have the stock equation)
    ode.structs.P$K1[1,] <- ode.structs.Q$K1[1,]
    ode.structs.P$K0[1] <- ode.structs.Q$K0[1]
    
    # correct stock part
    for (nn in 1:N.factors) {
      phi <- params.P[[as.character(nn)]]$phi
      rho <- params.P[[as.character(nn)]]$rho
      lmb <- params.P[[as.character(nn)]]$lmb
      
      kpp.P <- params.P[[as.character(nn)]]$kpp
      kpp.Q <- params.Q[[as.character(nn)]]$kpp
      eta.P <- params.P[[as.character(nn)]]$eta
      eta.Q <- params.Q[[as.character(nn)]]$eta
      
      ode.structs.P$K1[1,1+nn] <- ode.structs.P$K1[1,1+nn] +  phi * rho * (kpp.Q - kpp.P) / lmb[1] 
      ode.structs.P$K1[1,1+nn] <- ode.structs.P$K1[1,1+nn] + phi * sqrt(1-rho^2) * params.P[[as.character(1)]]$erp
      # add the constant part of the erp that is due to the vrp
      ode.structs.P$K0[1] <- ode.structs.P$K0[1] + phi * rho * (kpp.P * eta.P - kpp.Q * eta.Q) / lmb[1]
    }
    # add identically constant part of the erp
    ode.structs.P$K0[1] <- ode.structs.P$K0[1] + params.P[[as.character(1)]]$erp.0
    
    # correct coefficients in K1 for drift from other factors in the cascade model
    if(mod.type == "cascade.vol"){
      kpp.P <- params.P[["1"]]$kpp
      kpp.Q <- params.Q[["1"]]$kpp
      ode.structs.P$K1[2,3:(3+length(params.Q[["1"]]$lmb)-2)] <- - (kpp.P - kpp.Q)/params.Q[["1"]]$lmb[1] * params.Q[["1"]]$lmb[-1]
    }  
  } else {
    ode.structs.P <- ODEstructs(params = params.Q,jumpTransform = jumpTransform,mkt = mkt,N.factors = N.factors,mod.type = mod.type)
    params.P <- params.Q
  }
  
  solMat <- solveExtendedODE(u, mkt, ode.structs.P$K0, ode.structs.P$K1, ode.structs.P$l0, ode.structs.P$l1, ode.structs.P$H1, jmp = params.P$jmp, mf = mf, rtol=rtol, atol=atol, N.factors=N.factors, jumpTransform = jumpTransform, ...)
  
  return(solMat)
}