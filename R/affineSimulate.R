#' This function simulates paths from a stochastic volatility model with three volatility factors. One of the factors and the underlying asset prices can co-jump. The jump specification follows Duffie, Pan and Singleton (2000), p. 1361.
#' @param paramsList A list containing fields \code{P} and \code{Q}, containing, respectively, the P and Q parameters.
#' @param N.factors number of stochastic volatility factors (stock is not a factor)
#' @param t.days Number of days to simulate
#' @param t.freq Number of intraday prices to simulate
#' @param freq.subdiv Number of periods each segment should be broken into when simulating
#' @param rng.seed The random seed used for generations
#' @param init.vals Initial values for the of the stock/volatility to start from (in case we want to continue a given path, if NULL stock will be started from 1 and the volatilities from stationary distribution.)
#' @param rf.rate Assumed constant risk-free rate (annualized, continuously compounded)
#' @param jumpTransform The function that calculates the jump Laplace transform
#' @param specMaker function that takes arguments \code{params.P}, \code{params.Q}, \code{jumpTransform}, \code{N.factors}, \code{rf.rate}, \code{mod.type} and potentially more via the ellipsis operator (\code{...}) and returns the DPS matrices. See documentation of \link{ODEstructs}.
#' @export
#' @useDynLib affineModelR
#' @return A list containing the simulated stock & volatility paths + stock price jump times
#' 

affineSimulate <- function(paramsList, N.factors = 3, t.days = 1, t.freq = 1/78, freq.subdiv = 24, rng.seed = 42, init.vals = NULL, rf.rate = 0, jumpGeneratorPtr = getPointerToGenerator(fstr = 'expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform(fstr = 'expNormJumpTransform')$TF, mod.type = "standard", specMaker = ODEstructsForSim, nrepl = 1,...){
  
  # Set random seed
  set.seed(rng.seed)
  
  # sanity checks: are volatility inits of the right size?
  if(!is.null(init.vals)){
    stopifnot(length(init.vals$V.array) == N.factors)
  }
  
  # Make parameter structures which are to be passed to the C propagation code.
  paramsListCpp <- specMaker(paramsList$P, paramsList$Q, jumpTransformPtr, N.factors, rf.rate, mod.type = mod.type,...)
  
  # Obtain pointer to jump generator: you can provide your own c++ jump generating function
  jmpPtr = jumpGeneratorPtr
  
  # Set up time grid & random numbers for the BMs
  #   time.grid <- simulator2fMakeTimeGrid(t.days,t.freq,freq.subdiv,N.factors)
  h <- t.freq / freq.subdiv
  time.grid.length <- floor(t.days/h)
  dt <- h/252
  TT <- time.grid.length
  
  #   random.number.grids <- simulator2fGenerateRandomDraws(rng.seed,time.grid,N.factors)
  
  init.vals.cpp <- list("S.array" = 0, "V.array" = rep(1,N.factors))
  day.offset <- 0
  # Initialise stock array to 0 (we will be simulating log prices)
  if(!is.null(init.vals)){
    init.vals.cpp$S.array <- init.vals$S.array
    init.vals.cpp$V.array <- init.vals$V.array
    day.offset <- init.vals$day.offset
  }
  
  # Initialise V.array values from the stationary distribution of the factors
  if(is.null(init.vals)){
    if(!is.null(paramsList$P)){
      for(nn in 1:N.factors){
        der.mat <- matrix(0,nrow=2,ncol = N.factors+1)
        der.mat[2,nn+1] <- 1e-4
        cf <- affineCF(u = der.mat, params.Q = paramsList$Q, params.P = paramsList$P, t.vec = 20, v.0 = matrix(1,nrow=1,ncol=N.factors), jumpTransform = jumpTransformPtr, N.factors = N.factors, CGF = F, mod.type = mod.type, rtol = 1e-6, atol = 1e-12)
        cf <- drop(cf)
        cf.mom <- 1e4*diff(cf)
        init.vals.cpp$V.array[nn] <- cf.mom
      }
    } else {
      for(nn in 1:N.factors){
        der.mat <- matrix(0,nrow=2,ncol = N.factors+1)
        der.mat[2,nn+1] <- 1e-4
        cf <- affineCF(u = der.mat, params.Q = paramsList$Q, params.P = NULL, t.vec = 20, v.0 = matrix(1,nrow=1,ncol=N.factors), jumpTransform = jumpTransformPtr, N.factors = N.factors, CGF = F, mod.type = mod.type, rtol = 1e-6, atol = 1e-12)
        cf <- drop(cf)
        cf.mom <- 1e4*diff(cf)
        init.vals.cpp$V.array[nn] <- cf.mom
      }
    }
  } else {
    init.vals.cpp$V.array <- init.vals$V.array
    init.vals.cpp$S.array <- init.vals$S.array
  }
  
  # Pass the stuff to cpp and let it do the hard work.
  # simArrays <- .Call("affineOption_propagate3f", TT, 2*N.factors, paramsListCpp, dt, init.vals.cpp, rng.seed, PACKAGE = "affineOption")
  simArraysList <- vector(length = nrepl, mode = "list")
  V.array.list <- vector(length = nrepl, mode = "list")
  S.array.list <- vector(length = nrepl, mode = "list")
  numJumpsList <- numeric(nrepl)
  jumpTimesList <- vector(length = nrepl, mode = "list")
  for(kk in 1:nrepl){
    simArraysList[[kk]] <- affineSimulateCpp(TT, 2*N.factors, paramsListCpp, dt, init.vals.cpp, jmpPtr)
    
    # Thin the results to the desired frequency, add timestamps
    dt.loc <- simArraysList[[kk]]$dt*252*freq.subdiv
    
    pick.day.freq <- seq(from = 1, to = t.days/t.freq*freq.subdiv, by = freq.subdiv)
    
    V.array.list[[kk]] <- cbind((seq(0+day.offset,t.days+day.offset-dt.loc,by = dt.loc)),simArraysList[[kk]]$V.array[pick.day.freq,])
    S.array.list[[kk]] <- cbind((seq(0+day.offset,t.days+day.offset-dt.loc,by = dt.loc)),simArraysList[[kk]]$S.array[pick.day.freq])
    
    numJumpsList[kk] <- simArraysList[[kk]]$num.jumps
    
    loc.jumpTimes <- simArraysList[[kk]]$jump.times[which(simArraysList[[kk]]$jump.times!=0)]*252 # given in DAYS
    loc.jumpSizes <- simArraysList[[kk]]$jump.sizes[,which(simArraysList[[kk]]$jump.times!=0),drop=FALSE]
    loc.jumpDF <- as.data.frame(cbind(loc.jumpTimes, t(loc.jumpSizes)))
    colnames(loc.jumpDF) <- c("day","J_logS",sapply(1:(nrow(loc.jumpSizes)-1), function(ss) paste0("J_v",ss)))
    jumpTimesList[[kk]] <- as.matrix(loc.jumpDF)
    
    colnames(V.array.list[[kk]]) <- c("day",paste0("v",1:N.factors))
    colnames(S.array.list[[kk]]) <- c("day","F")
  }
  require(abind)
  V.array <- do.call(what = abind, args = V.array.list)
  S.array <- do.call(what = abind, args = S.array.list)
  numJumps <- do.call(what = c, args = as.list(numJumpsList))
  
  return(list(S.array = S.array, V.array = V.array, dt=dt/252, numJumps = numJumps, jumpSizes = jumpTimesList))
}