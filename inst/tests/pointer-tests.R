### test all the pointer passing

library(affineModelR)
library(testthat)
N.factors <- 1

params <- list()
params[["1"]] <- list(phi=1,eta=.04,lmb=0.15,kpp=2,rho=0)
params$jmp <- list(muSc = 0.05, muYc = 0, sigmaYc=1e-2,lprop=1,lvec=0,rhoc=-0.2)

params.P <- params
params.P[["1"]]$eta <- 0.03
params.P[["1"]]$kpp <- 3
params.P$jmp <- list(muSc = 0.07, muYc = -0.03, sigmaYc=3e-2,lprop=1,lvec=0,rhoc=-0.2)
params.P$`1`$erp <- 1

uVec <- 1:10

mkt <- list(p=1,q=-.5,r=.2,t=.15)

jumpTrPtrs <- getPointerToJumpTransform(fstr = 'expNormJumpTransform')

ode.structs.check <- ODEstructs(params = params, jumpTransform = jumpTrPtrs, mkt = mkt, N.factors = 1, mod.type = 'standard')

evaluateTransform(genPtr_ = jumpTrPtrs$TF, beta = c(1,0), jmpPar = params$jmp)

sol.check <- jumpDiffusionODEs(u = cbind(0,uVec), params = params, mkt = mkt, jumpTransform = jumpTrPtrs$TF, N.factors = 1, mod.type = 'standard', rtol = 1e-6, atol = 1e-12)

sol.check.P <- jumpDiffusionODEsP(u = cbind(0, uVec), params.P = params.P, params.Q = params, mkt = mkt, jumpTransform = jumpTrPtrs$TF, N.factors = 1, mod.type = 'standard', rtol = 1e-6, atol = 1e-12)

sol.check.ext <- odeExtSolveWrap(u = cbind(0,uVec), params.Q = params, params.P = NULL, mkt = mkt, rtol = 1e-8, atol = 1e-12, N.factors = 1, jumpTransform = jumpTrPtrs, mod.type = 'standard')

cf.check <- affineCF(u = cbind(0,uVec), params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), jumpTransform = jumpTrPtrs$TF, N.factors = 1, CGF = F, mod.type = 'standard')

cf.derivs.check <- affineCFderivs(u = cbind(0,uVec), params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), jumpTransform = jumpTrPtrs, N.factors = 1, CGF = F, mod.type = 'standard')

### parallel ? ###
library(parallel)
cl <- makeCluster(2)
clusterEvalQ(cl,library(affineModelR))

clusterExport(cl, c("uVec","mkt","params","params.P"))
clusterEvalQ(cl, ls())

affineCFpar <- function(u, params.Q, params.P = NULL, t.vec, v.0, N.factors = 3, CGF = FALSE, mod.type = "standard", ptrFoo, ptrArg){
  return(affineCF(t(u),params.Q,params.P,t.vec,v.0,N.factors,CGF,mod.type,jumpTransform = eval(expr = quote(ptrFoo(ptrArg)$TF))))
}

clusterExport(cl,'affineCFpar')

clusterEvalQ(cl, affineCFpar(u = cbind(0,uVec), params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), N.factors = 1, CGF = F, mod.type = 'standard', ptrFoo = getPointerToJumpTransform, ptrArg = 'kouExpJumpTransform'))

# the following does not work
clusterEvalQ(cl, some.pointer <- getPointerToJumpTransform('expNormJumpTransform'))

# in general can't get pointers on the cluster
call.check <- apply(X = cbind(0,uVec), MARGIN = 1, FUN = function(u) affineCF(u = t(u), params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), jumpTransform = getPointerToJumpTransform(fstr = 'expNormJumpTransform')$TF, N.factors = 1))

parallel.check <- parApply(cl = cl, X = cbind(0,uVec), MARGIN = 1, FUN = function(u) affineCF(u = t(u), params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), jumpTransform = getPointerToJumpTransform(fstr = 'kouExpJumpTransform')$TF, N.factors = 1))

parallelWrapper <- function(foo, ptrFoo, ptrArg, argName, ptrElem = NULL, ...){
  if(is.null(ptrElem)){
      eval(parse(text = paste0('resFoo <- function(...){
        res <- foo(..., ' , argName ,' = ptrFoo(ptrArg))
        return(res)
      }
      return(resFoo)')))
  } else {
    eval(parse(text = paste0('resFoo <- function(...){
        res <- foo(..., ' , argName ,' = ptrFoo(ptrArg)$',ptrElem,')
                             return(res)
        }
        return(resFoo)')))
    }
}

cfPar <- parallelWrapper(foo = affineCFderivs, ptrFoo = getPointerToJumpTransform, ptrArg = 'expNormJumpTransform', argName = 'jumpTransform')

call.check.2 <- lapply(X = uVec.list, FUN = cfPar, params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), N.factors = 1, mod.type ='standard')


parallel.check.2 <- parLapply(cl = cl, X = uVec.list, fun = cfPar, params.Q = params, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), N.factors = 1, mod.type ='standard')


uVec.list <- list(cbind(0,uVec),cbind(uVec,0))


### simulate some ###

paths <- affineSimulate(paramsList = list(P = params.P, Q = params), N.factors = 1, t.days = 10, t.freq = 1/78, freq.subdiv = 10, jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = jumpTrPtrs$TF, specMaker = ODEstructsForSim, mod.type = 'standard', nrepl = 1)
