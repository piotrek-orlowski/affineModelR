#' @title ODE specification structures
#' @description This function translates the parameters lists (see \code{jumpDiffusionODEs}) to model specification matrices and vectors for ODE solution. Expand this function if you want to add other specifications.
#' @param params A structure containing the volatility factor + jump descriptions.
#' @param jumpTransform String, the function that calculates the jump transform.
#' @param N.factors The number of stochastic volatility factors (stock is not a factor).
#' @param mod.type string. \code{"standard"} forces the basic DPS 'concrete example' parametrisation, \code{"cascade.vol"} allows for richer vol-of-vol structure
#' @details using \code{mod.type = "cascade.vol"} allows that the vol of vol of the first factor is an affine function of other factors. The constraint is that this factor can't have a leverage effect (i.e. its rho parameter must be 0), otherwise you're not in the affine setting.
#' @export
#' @return Return a list with elements K0, l0, K1, l1, H1, see DPS (2000) for their meaning

ODEstructs <- function(params,jumpTransform,mkt,N.factors,mod.type = "standard") {
  
  # Error control
  
  stopifnot(length(params$jmp$lvec) == 1)
  
  stopifnot(length(params$jmp$lprop) == N.factors)
  
  if(mod.type == "cascade.vol" & !(params[["1"]]$rho==0)){
    stop("The first factor rho is not 0")
  }
  
  if(jumpTransform == "expNormJumpTransform"){
    jumpTransform <- expNormJumpTransform
  } else if(jumpTransform == "kouExpJumpTransform"){
    jumpTransform <- kouExpJumpTransform
  }
  
  #drift that depends on state variables
  K1 <- matrix(0,N.factors+1,N.factors+1)
  
  for (nn in 1:N.factors) {
    K1[1,1+nn] = -0.5 * params[[as.character(nn)]]$phi^2 - params$jmp$lprop[nn]*(jumpTransform(c(1,0,0),params$jmp))
    K1[1+nn,1+nn] = -params[[as.character(nn)]]$kpp
  }
  
  l0 = params$jmp$lvec
  
  K0 <- rep(0, N.factors + 1)
  # Here we explicitly IGNORE mkt$r and mkt$q in the construction of the drift matrix for the ODE solver. The adjustments for r and q are made inside solvODE, after the solution has been obtained.
  K0[1] <- -l0 * (jumpTransform(c(1,rep(0,N.factors)),params$jmp))
  
  for (nn in 1:N.factors) {
    K0[1+nn] <- params[[as.character(nn)]]$eta*params[[as.character(nn)]]$kpp
  }
  
  H1 = array(0,rep(N.factors+1,3))
  
  if(mod.type == "standard"){
    for (nn in 1:N.factors) {
      H1[1,1,1+nn]  <- params[[as.character(nn)]]$phi^2
      H1[1,1+nn,1+nn] <- params[[as.character(nn)]]$phi * params[[as.character(nn)]]$rho * params[[as.character(nn)]]$lmb
      H1[1+nn,1,1+nn] <- params[[as.character(nn)]]$phi * params[[as.character(nn)]]$rho * params[[as.character(nn)]]$lmb
      H1[1+nn,1+nn,1+nn] <- params[[as.character(nn)]]$lmb^2
    }
  } else if(mod.type == "cascade.vol"){
    for (nn in 1:N.factors) {
      H1[1,1,1+nn]  <- params[[as.character(nn)]]$phi^2
      H1[1,1+nn,1+nn] <- params[[as.character(nn)]]$phi * params[[as.character(nn)]]$rho * params[[as.character(nn)]]$lmb[1]
      H1[1+nn,1,1+nn] <- params[[as.character(nn)]]$phi * params[[as.character(nn)]]$rho * params[[as.character(nn)]]$lmb[1]
      H1[1+nn,1+nn,1+nn] <- params[[as.character(nn)]]$lmb[1]^2
    }
    H1[2,2,(2:(1+length(params[["1"]]$lmb)))] <- params[["1"]]$lmb
  }
  
  l1 = c(0, params$jmp$lprop)
  return(list(K0 = K0, l0 = l0, K1 = K1, l1 = l1, H1 = H1))
}