#' @title ODE specification structures
#' @description This function translates the parameters lists (see \code{jumpDiffusionODEs}) to model specification matrices and vectors for ODE solution. Expand this function if you want to add other specifications.
#' @param params A structure containing the volatility factor + jump descriptions.
#' @param jumpTransform String, the function that calculates the jump transform.
#' @param N.factors The number of stochastic volatility factors (stock is not a factor).
#' @export
#' @return Return a list with elements K0, l0, K1, l1, H1, see DPS (2000) for their meaning

ODEstructs <- function(params, jumpTransform = NULL, mkt, N.factors) {
  
  # for compatibility, leave this mechanism whenever jumpTransform is provided
  if(!is.null(jumpTransform)){
    
    # Error control
    stopifnot(length(params$jmp$lvec) == 1)
    
    stopifnot(length(params$jmp$lprop) == N.factors)
    
    if(inherits(jumpTransform,'list')){
      jumpTrPtr <- list(jumpTransform$TF)
    } else if(inherits(jumpTransform,'externalptr')){
      jumpTrPtr <- list(jumpTransform)
    } 
    
    jmp_lvec <- matrix(params$jmp$lvec)
    jmp_lprop <- matrix(params$jmp$lprop, ncol = 1L)
    
  } else {
    par_list_names <- names(params)
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
  }
  # How many different jump transforms?
  N.jumps <- length(jumpTrPtr)
  
  #drift that depends on state variables
  K1 <- matrix(0,N.factors+1,N.factors+1)
  
  for (nn in 1:N.factors) {
    K1[1,1+nn] <- -0.5 * params[[as.character(nn)]]$phi^2
    for(jj in 1:N.jumps){
      # stock drift components under Q -- jump compensation (P correction is done outside)
      K1[1,1+nn] <- K1[1,1+nn] - jmp_lprop[nn, jj]*Re(evaluateTransform(genPtr_ = jumpTrPtr[[jj]], beta = c(1,rep(0,N.factors)),jmpPar = params[[sprintf("jmp%d", jj)]])) 
    }
    # vol components
    K1[1+nn,1+nn] <- -params[[as.character(nn)]]$kpp
    # use the tv_lr_mean flag to work with factors that have a long-run time-varying mean
    if(!is.null(params[[as.character(nn)]]$tv_lr_mean)){
      if(params[[as.character(nn)]]$tv_lr_mean){
        K1[1+nn,1+nn+1] <- params[[as.character(nn)]]$kpp
      }
    }
  }
  
  K0 <- rep(0, N.factors + 1)
  # Here we explicitly IGNORE mkt$r and mkt$q in the construction of the drift matrix for the ODE solver. The adjustments for r and q are made inside solvODE, after the solution has been obtained.
  for(jj in 1:N.jumps){
    K0[1] <- K0[1] - jmp_lvec[jj] * Re(evaluateTransform(genPtr_ = jumpTrPtr[[jj]], beta = c(1,rep(0,N.factors)),jmpPar = params[[sprintf("jmp%d", jj)]]))
  }
  
  for (nn in 1:N.factors) {
    K0[1+nn] <- params[[as.character(nn)]]$eta*params[[as.character(nn)]]$kpp
    # use the tv_lr_mean flag to work with factors that have a long-run time-varying mean
    if(!is.null(params[[as.character(nn)]]$tv_lr_mean)){
      if(params[[as.character(nn)]]$tv_lr_mean){
      K0[1+nn] <- 0.0   
      }
    }
  }
  
  H1 = array(0,rep(N.factors+1,3))
  
  
  for (nn in 1:N.factors) {
    H1[1,1,1+nn]  <- params[[as.character(nn)]]$phi^2
    H1[1,1+nn,1+nn] <- params[[as.character(nn)]]$phi * params[[as.character(nn)]]$rho * params[[as.character(nn)]]$lmb
    H1[1+nn,1,1+nn] <- params[[as.character(nn)]]$phi * params[[as.character(nn)]]$rho * params[[as.character(nn)]]$lmb
    H1[1+nn,1+nn,1+nn] <- params[[as.character(nn)]]$lmb^2
  }
  
  return(list(K0 = K0, l0 = jmp_lvec, K1 = K1, l1 = jmp_lprop, H1 = H1))
}