#' @title C++ extended ODE solve wraper
#' @description This function passes the ODE parameters, in a specification defined via \code{\link{ODEstructs}}, to the C++ solution routines.
#' @param u u \code{Ux(N.factors+1)} matrix of complex numbers
#' @param mkt data.frame describing the market structure (times to maturity, interest rate, dividend yield, current stock price). See Details in \code{\link{jumpDiffusionODEs}}.
#' @param K0,K1,l0,l1,H1 output from \code{\link{ODEstructs}}, for details see DPS (2000).
#' @param jmp jump parameters, element of \code{params} list from \code{\link{jumpDiffusionODEs}}.
#' @param jumpTransform pointer to jump transform function.
#' @param N.factors number of stochastic volatility factors.
#' @param mf,rtol,atol ODE solution precision parameters, see \code{\link{jumpDiffusionODEs}}.
#' @details This is an internal function not intended for the user.
#' 

solveODE <- function(u, mkt, K0, K1, l0, l1, H1, jmp, jumpTransform = getPointerToJumpTransform(fstr = 'expNormJumpTransform')$TF, mf = 22, rtol=1e-14, atol=1e-30, N.factors = 3) {  
  # prepare structures where we save
  N <- nrow(u)
  TT <- dim(mkt$t)[1]
  if(is.null(TT)){TT <- length(mkt$t)}
  
  solMat <- array(0,c(N,TT,N.factors+1))
  # set column names for the solutions, b1 is the first element of the vector beta
  dimnames(solMat) <- list(paste("u",seq(1:N),sep="."),mkt$t,c(paste("b",1:N.factors,sep=""),"a"))
  
  # load C parameter passing structures
  odeList <- list()
  odeList$N.factors <- N.factors
  odeList$H1r <- Re(as.vector(H1))
  odeList$H1i <- Im(as.vector(H1))
  odeList$K1r <- Re(K1)
  odeList$K1i <- Im(K1)
  
  odeList$K0 <- K0
  odeList$l0 <- l0
  odeList$l1 <- l1
  odeList$jmpPar <- jmp
  odeList$jumpTransformPtr <- jumpTransform

  # now solve for all frequencies
  for (uu in 1:nrow(u)) {
    u.vec <- u[uu,]
    # initialize ODE
    betaalpha0 <- c(u.vec,0)
    names(betaalpha0) <- c("c1",c(paste("b",1:N.factors,sep=""),"a"))
    
    sol <- zvode(betaalpha0, times=c(0,mkt$t), func = "derivs2", parms = odeList, dllname = "affineModelR", initfunc = "initmod",nout=0, mf=mf, rtol=rtol,atol=atol,maxsteps=1500000)
    solMat[uu,,] <- sol[-1,c(paste("b",1:N.factors,sep=""),"a")]
  }
  
  # make correction for the part of drift implied by (r-q) (interest rate less dividend yield)
  drift.mat <- array(0, dim= dim(solMat))
  drift.mat[,,N.factors+1] <- u[,1,drop= FALSE] %*% matrix(mkt$t*(mkt$r - mkt$q), nrow= 1)
  solMat <- solMat + drift.mat
  
  return(solMat)
}