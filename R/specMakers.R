#' Calculate vectors and matrices for propagation equations in a N-factor model with co-jumps in the underlying and the first volatility factor
#' @param params.P A structure containing the volatility factor + jump descriptions under measure P
#' @param params.Q A structure containing the volatility factor + jump descriptions under measure Q
#' @param jumpTransform Pointer to the basic jumpTransform c++ function. Via this mechanism it is easy for the user to provide their external jump transform functions for both simulation and ODE calculations.
#' @param N.factors The number of stochastic volatility factors (stock is not a factor)
#' @param rf.rate risk-free rate
#' @export
#' @return Return a list with elements \code{terms.dt}, \code{terms.vdt}, \code{terms.vdW}, \code{terms.vdWort}, \code{vterms.dt}, \code{vterms.vdt}, \code{vterms.vdW}, \code{jmpPar}, \code{intensity.terms}. These are all elements for simulating the N-factor model.

ODEstructsForSim <- function(params.P = NULL, params.Q, jumpTransformPointer = getPointerToJumpTransform(fstr = 'expNormJumpTransform'), N.factors, rf.rate = 0.0,mod.type = "standard") {
  
  ### if no params.P structure is given, the simulation is done under Q!
  doP <- !is.null(params.P)
  if(doP){
    stopifnot(length(params.P$jmp$lvec) == 1)
    stopifnot(length(params.P$jmp$lprop) == N.factors)
  }
  stopifnot(length(params.Q$jmp$lvec) == 1)
  stopifnot(length(params.Q$jmp$lprop) == N.factors)
  
  qJmpCompensator <- Re(evaluateTransform(genPtr_ = jumpTransformPointer, beta = c(1,rep(0,N.factors)), jmpPar = params.Q$jmp))
  
  ### Prepare the Q terms.dt
  
  ## arrays and vectors of parameters for the stock equation propagation
  # terms.dt : constants, i.e. risk-free rate plus jump compensatro times constant intensity
  terms.dt <- rf.rate - params.Q$jmp$lvec * qJmpCompensator
  
  ### Prepare the P terms.dt
  if(doP){
    for(nn in 1:N.factors){
      loc.g1 <- -(params.Q[[as.character(nn)]]$kpp * params.Q[[as.character(nn)]]$eta - params.P[[as.character(nn)]]$kpp * params.P[[as.character(nn)]]$eta) / params.P[[as.character(nn)]]$lmb[1]
      #     print(loc.g1)
      terms.dt <- terms.dt + params.P[[as.character(nn)]]$phi * params.P[[as.character(nn)]]$rho * loc.g1
    }
  }
  
  # terms.vdt: what in the stock price equation is multiplied by vt*dt
  ### Q terms.vdt
  terms.vdt <- rep(0,N.factors)
  for(nn in 1:N.factors){
    terms.vdt[[nn]] <- -0.5*params.Q[[as.character(nn)]]$phi^2 - params.Q$jmp$lprop[nn] * qJmpCompensator
  }
  
  ### P terms.vdt: two risk-premium components, one from the erp parameter -- price of risk on the stock-specific BM, the other from the price of variance risk through the leverage effect, price of risk on the variance-specific BM
  if(doP){
    for(nn in 1:N.factors){
      loc.g2 <- -(params.P[[as.character(nn)]]$kpp - params.Q[[as.character(nn)]]$kpp)/params.P[[as.character(nn)]]$lmb[1]
      #     print(loc.g2)
      terms.vdt[[nn]] <- terms.vdt[[nn]] + params.P[[as.character(nn)]]$phi * params.P[[as.character(nn)]]$rho * loc.g2 + params.P[[as.character(nn)]]$phi * sqrt(1 - params.P[[as.character(nn)]]$rho^2) * params.P[[as.character(1)]]$erp
    }
  }
  
  ### terms.vdW: what gets multiplied by the sqrt of vol state and by BM driver. Identical under P and Q
  terms.vdW <- terms.vdWort <- rep(0,N.factors)
  for(nn in 1:N.factors){
    terms.vdW[nn] <- params.Q[[as.character(nn)]]$phi * params.Q[[as.character(nn)]]$rho
    terms.vdWort[nn] <- params.Q[[as.character(nn)]]$phi * sqrt(1 - params.Q[[as.character(nn)]]$rho^2)
  }
  
  ### Deal with P and Q measure jumps
  if(doP){
    terms.intensity <- rep(0,N.factors + 1)
    terms.intensity[1] <- params.P$jmp$lvec
    terms.intensity[2:(N.factors+1)] <- params.P$jmp$lprop
    jmpPar <- unlist(params.P$jmp[c("muYc","sigmaYc","muSc","rhoc")])
  } else {
    terms.intensity <- rep(0,N.factors + 1)
    terms.intensity[1] <- params.Q$jmp$lvec
    terms.intensity[2:(N.factors+1)] <- params.Q$jmp$lprop
    jmpPar <- unlist(params.Q$jmp[c("muYc","sigmaYc","muSc","rhoc")])
  }
  
  ## volatility part
  vterms.vdW <- vterms.dt <- vterms.vdt <- matrix(0,N.factors,N.factors)
  if(mod.type == "standard"){
    if(doP){
      for(nn in 1:N.factors){
        vterms.dt[nn,nn] <- params.P[[as.character(nn)]]$kpp * params.P[[as.character(nn)]]$eta
        vterms.vdt[nn,nn] <- - params.P[[as.character(nn)]]$kpp
        # lmb below has to be squared in the 'standard' setup because vterms.dW enters under the square root for the needs of cascade sim
        vterms.vdW[nn,nn] <- params.P[[as.character(nn)]]$lmb[1]^2
      }
    } else {
      for(nn in 1:N.factors){
        vterms.dt[nn,nn] <- params.Q[[as.character(nn)]]$kpp * params.Q[[as.character(nn)]]$eta
        vterms.vdt[nn,nn] <- - params.Q[[as.character(nn)]]$kpp
        # lmb below has to be squared in the 'standard' setup because vterms.dW enters under the square root for the needs of cascade sim
        vterms.vdW[nn,nn] <- params.Q[[as.character(nn)]]$lmb[1]^2
      }
    }
  } else if(mod.type == "cascade.vol"){
    if(doP){
      for(nn in 1:N.factors){
        vterms.dt[nn,nn] <- params.P[[as.character(nn)]]$kpp * params.P[[as.character(nn)]]$eta
        vterms.vdt[nn,nn] <- - params.P[[as.character(nn)]]$kpp
        vterms.vdW[nn,nn] <- params.P[[as.character(nn)]]$lmb[1]^2
      }
      vterms.vdt[2:(2+length(params.P[["1"]]$lmb)-2),1] <- - (params.P[["1"]]$kpp - params.Q[["1"]]$kpp)/params.P[["1"]]$lmb[1] * params.P[["1"]]$lmb[-1]
      vterms.vdW[1:length(params.P[["1"]]$lmb),1] <- params.P[["1"]]$lmb
    } else {
      for(nn in 1:N.factors){
        vterms.dt[nn,nn] <- params.Q[[as.character(nn)]]$kpp * params.Q[[as.character(nn)]]$eta
        vterms.vdt[nn,nn] <- - params.Q[[as.character(nn)]]$kpp
        vterms.vdW[nn,nn] <- params.Q[[as.character(nn)]]$lmb[1]^2
      }
      vterms.vdW[1:length(params.Q[["1"]]$lmb),1] <- params.Q[["1"]]$lmb
    }
  }
  
  #  intensity.terms <- c(params.P$jmp$lvec,params.P$jmp$lprop)
  
  return(list(terms.dt = terms.dt, terms.vdt = terms.vdt, terms.vdW = terms.vdW, terms.vdWort = terms.vdWort, vterms.dt = vterms.dt, vterms.vdt = vterms.vdt, vterms.vdW = vterms.vdW, jmpPar = jmpPar, intensity.terms = terms.intensity))
}