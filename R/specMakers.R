#' Calculate vectors and matrices for propagation equations in a N-factor model with co-jumps in the underlying and the first volatility factor
#' @param params.P A structure containing the volatility factor + jump descriptions under measure P
#' @param params.Q A structure containing the volatility factor + jump descriptions under measure Q
#' @param jumpTransform Pointer to the basic jumpTransform c++ function. Via this mechanism it is easy for the user to provide their external jump transform functions for both simulation and ODE calculations.
#' @param N.factors The number of stochastic volatility factors (stock is not a factor)
#' @param rf.rate risk-free rate
#' @export
#' @return Return a list with elements \code{terms.dt}, \code{terms.vdt}, \code{terms.vdW}, \code{terms.vdWort}, \code{vterms.dt}, \code{vterms.vdt}, \code{vterms.vdW}, \code{jmpPar}, \code{intensity.terms}. These are all elements for simulating the N-factor model.

ODEstructsForSim <- function(params.P = NULL, params.Q, jumpTransformPointer = getPointerToJumpTransform(fstr = 'expNormJumpTransform'), N.factors, rf.rate = 0.0) {
    
  
  if(!is.null(jumpTransformPointer)){
    
    ### if no params.P structure is given, the simulation is done under Q!
    doP <- !is.null(params.P)
    if(doP){
      stopifnot(length(params.P$jmp$lvec) == 1)
      stopifnot(length(params.P$jmp$lprop) == N.factors)
    }
    stopifnot(length(params.Q$jmp$lvec) == 1)
    stopifnot(length(params.Q$jmp$lprop) == N.factors)
    
    if(inherits(jumpTransformPointer,'list')){
      jumpTrPtr <- list(jumpTransformPointer$TF)
    } else if(inherits(jumpTransform,'externalptr')){
      jumpTrPtr <- list(jumpTransformPointer)
    } 
    
    if(doP){
      jmp_lvec <- matrix(params.P$jmp$lvec)
      jmp_lprop <- matrix(params.P$jmp$lprop, ncol = 1L)
    } else {
      jmp_lvec <- matrix(params.Q$jmp$lvec)
      jmp_lprop <- matrix(params.Q$jmp$lprop, ncol = 1L) 
    }
    
  } else {
    par_list_names <- names(params)
    
    if(doP){
      params <- params.P
    } else {
      params <- params.Q
    }
    
    # constant part of each jump intensity
    jmp_lvec <- do.call(what = cbind
                        , args = lapply(X = params[grepl("jmp", par_list_names)]
                                        , FUN = function(jmp_list){
                                          jmp_list$lvec
                                        }))
    # time-varying part of each jump intensity
    jmp_lprop <- do.call(what = cbind
                         , args = lapply(X = params[grepl("jmp", par_list_names)]
                                         , FUN = function(jmp_list){
                                           jmp_list$lprop
                                         }))
    # transform pointers -- only TF element (no derivatives)
    jumpTrPtr <- lapply(X = params[grepl("jmp", par_list_names)]
                        , FUN = function(jmp_list){
                          jmp_list$jumpTransform$TF
                        })
    
    rm(params)
  }
  # How many different jump transforms?
  N.jumps <- length(jumpTrPtr)
   
  # Q-expectations of all jumps
  # qJmpCompensator <- Re(evaluateTransform(genPtr_ = jumpTransformPointer, beta = c(1,rep(0,N.factors)), jmpPar = params.Q$jmp))
  qJmpCompensator <- sapply(seq_along(jumpTrPtr)
                            , function(jmp_index){
                              evaluateTransform(genPtr_ = jumpTrPtr[[jmp_index]]
                                                , beta = beta = c(1,rep(0,N.factors))
                                                , params.Q[[sprintf("jmp%d", jj)]])
                            })
  
  ### Prepare the Q terms.dt
  
  ## arrays and vectors of parameters for the stock equation propagation
  # terms.dt : constants, i.e. risk-free rate plus jump compensatro times constant intensity
  terms.dt <- rf.rate - crossprod(jmp_lvec, qJmpCompensator)
  
  ### Prepare the P terms.dt
  if(doP){
    # if orthogonal erp not defined, then let's assume it is 0
    if (is.null(params.P[[as.character(1)]]$erp)) {
      params.P[[as.character(1)]]$erp <- 0
    }
    if (is.null(params.P[[as.character(1)]]$erp0)) {
      params.P[[as.character(1)]]$erp0 <- 0
    }
    
    if(N.factors > 1){
      for(nn in 2:N.factors){
        if(is.null(params.P[[as.character(nn)]]$erp)){
          params.P[[as.character(nn)]]$erp <- 0
        }
      }
    }
    
    # for(nn in 1:N.factors){
      # loc.g1 <- -(params.Q[[as.character(nn)]]$kpp * params.Q[[as.character(nn)]]$eta - params.P[[as.character(nn)]]$kpp * params.P[[as.character(nn)]]$eta) / params.P[[as.character(nn)]]$lmb[1]
      #     print(loc.g1)
      # terms.dt <- terms.dt + params.P[[as.character(nn)]]$phi * params.P[[as.character(nn)]]$rho * loc.g1
    # }
    terms.dt <- terms.dt + params.P[[as.character(1)]]$erp0
  }
  
  # terms.vdt: what in the stock price equation is multiplied by vt*dt
  ### Q terms.vdt
  terms.vdt <- rep(0,N.factors)
  for(nn in 1:N.factors){
    terms.vdt[[nn]] <- -0.5*params.Q[[as.character(nn)]]$phi^2 - jmp_lprop[nn,,drop=FALSE] %*% qJmpCompensator
  }
  
  ### P terms.vdt: two risk-premium components, one from the erp parameter -- price of risk on the stock-specific BM, the other from the price of variance risk through the leverage effect, price of risk on the variance-specific BM
  if(doP){
    for(nn in 1:N.factors){
      # loc.g2 <- -(params.P[[as.character(nn)]]$kpp - params.Q[[as.character(nn)]]$kpp)/params.P[[as.character(nn)]]$lmb[1]
      #     print(loc.g2)
      # terms.vdt[[nn]] <- terms.vdt[[nn]] + params.P[[as.character(nn)]]$phi * params.P[[as.character(nn)]]$rho * loc.g2 +  params.P[[as.character(nn)]]$erp
      terms.vdt[[nn]] <- terms.vdt[[nn]] +  params.P[[as.character(nn)]]$erp
    }
    # terms.vdt[[1]] <- terms.vdt[[1]] + params.P[[as.character(1)]]$erp0
  }
  
  ### terms.vdW: what gets multiplied by the sqrt of vol state and by BM driver. Identical under P and Q
  terms.vdW <- terms.vdWort <- rep(0,N.factors)
  for(nn in 1:N.factors){
    terms.vdW[nn] <- params.Q[[as.character(nn)]]$phi * params.Q[[as.character(nn)]]$rho
    terms.vdWort[nn] <- params.Q[[as.character(nn)]]$phi * sqrt(1 - params.Q[[as.character(nn)]]$rho^2)
  }
  
  ### Deal with P and Q measure jumps
  if(doP){
    terms.intensity <- matrix(0, N.factors + 1, N.jumps)
    terms.intensity[1,] <- jmp_lvec
    terms.intensity[2:(N.factors+1),] <- jmp_lprop
    
    jmpPar <- params.P[grepl("jmp", names(params.P))]
  } else {
    terms.intensity <- rep(0,N.factors + 1)
    terms.intensity <- matrix(0, N.factors + 1, N.jumps)
    terms.intensity[1,] <- jmp_lvec
    terms.intensity[2:(N.factors+1),] <- jmp_lprop
    
    jmpPar <- params.Q[grepl("jmp", names(params.Q))]
  }
  
  ## volatility part
  vterms.vdW <- vterms.dt <- vterms.vdt <- matrix(0,N.factors,N.factors)
  
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
 
  
  #  intensity.terms <- c(params.P$jmp$lvec,params.P$jmp$lprop)
  
  return(list(terms.dt = terms.dt, terms.vdt = terms.vdt, terms.vdW = terms.vdW, terms.vdWort = terms.vdWort, vterms.dt = vterms.dt, vterms.vdt = vterms.vdt, vterms.vdW = vterms.vdW, jmpPar = jmpPar, intensity.terms = terms.intensity))
}